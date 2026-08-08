using GadgetSearch
using Graphs
using Random
using Test
function cross_graph()
    graph = SimpleGraph(4)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 2, 4)
    return graph
end
function graph_from_edges(order, edge_list)
    graph = SimpleGraph(order)
    foreach(edge -> add_edge!(graph, edge...), edge_list)
    return graph
end
@testset "Unweighted search" begin
    @testset "fixed verifier" begin
        reduced = calculate_reduced_alpha_tensor(cross_graph(), [1, 2, 3, 4])
        @test GadgetSearch.inf_mask(reduced) == BigInt(60576)
        @test is_diff_by_constant(reduced .+ 3, reduced) == (true, 3.0)
    end
    @testset "crossing frame" begin
        square = [(0, 0), (-1, 0), (0, 1), (1, 0), (0, -1)]
        pins = square[2:5]
        rays = [(-1, 0), (0, 1), (1, 0), (0, -1)]
        @test all(check_crossing_frame(Square(), square, pins, rays))
        @test !check_crossing_frame(Square(), [square; (-2, 1)], pins, rays).G4
        triangular = [(0, 0), (-1, 0), (0, 1), (1, 0), (-1, -1)]
        triangular_pins = triangular[2:5]
        @test all(check_crossing_frame(Triangular(), triangular, triangular_pins, rays))
    end

    @testset "logical states are dynamic and planar" begin
        target = cross_graph()
        target_reduced = vec(calculate_reduced_alpha_tensor(target, [1, 2, 3, 4]))
        skeletons, evaluated, _, best, trace = GadgetSearch._search_logical_skeletons(
            target, [1, 2, 3, 4], target_reduced, Triangular();
            min_vertices=9, max_vertices=15, max_evaluations=120,
            beam_width=8, mutations_per_candidate=4,
            random_candidates_per_generation=3, exploration_fraction=0.25,
            rng=MersenneTwister(4),
        )
        @test evaluated == 120
        @test length(trace) == evaluated
        @test best[1] >= 0
        @test all(record -> record.stage == :logical, trace)
        @test all(record -> GadgetSearch._has_alternating_planar_frame(
            graph_from_edges(record.vertices, record.graph_edges),
            record.boundary_vertices,
        ), trace)
        @test all(state -> is_gadget_replacement(
            target, state.graph, [1, 2, 3, 4], state.boundary,
        )[1], skeletons)
    end

    @testset "exact rewrites preserve the tensor" begin
        logical = graph_from_edges(13, [(1,11),(2,5),(2,6),(2,7),(3,7),(3,10),
            (4,10),(4,12),(4,13),(5,6),(5,8),(5,11),(5,13),(6,7),(6,8),
            (7,8),(8,9),(9,10),(9,12),(9,13),(10,12),(11,13),(12,13)])
        @test is_gadget_replacement(cross_graph(), logical, [1,2,3,4], [1,2,3,4]) == (true, 3.0)

        split = GadgetSearch._split_vertex(logical, 5, [2, 6])
        @test is_gadget_replacement(logical, split, [1,2,3,4], [1,2,3,4]) == (true, 1.0)
        subdivided = GadgetSearch._even_subdivide_edges(logical, [(1, 11)])
        @test is_gadget_replacement(logical, subdivided, [1,2,3,4], [1,2,3,4]) == (true, 1.0)
    end

    @testset "triangular positive control" begin
        axial = [(0,0),(1,6),(3,5),(8,0),(0,6),(-1,6),(1,5),(-1,5),(2,5),
            (-1,4),(2,4),(4,4),(0,3),(1,3),(3,3),(5,3),(6,3),(0,2),(2,2),
            (3,2),(7,2),(4,1),(1,1),(5,1),(6,1),(7,1),(1,0),(2,0),(3,0),
            (4,0),(7,0),(5,-1),(8,-1),(8,-2),(6,-2),(8,-3),(7,-3)]
        coordinates = GadgetSearch._axial_to_offset.(axial)
        positions = GadgetSearch.get_physical_positions(Triangular(), coordinates)
        graph = GadgetSearch.unit_disk_graph(positions, get_radius(Triangular()))
        @test is_gadget_replacement(cross_graph(), graph, [1,2,3,4], [1,2,3,4]) == (true, 15.0)
        rays = [(-1,0), (0,1), (0,1), (1,0)]
        @test all(check_crossing_frame(Triangular(), coordinates, coordinates[1:4], rays))
        patch, placed, _ = GadgetSearch._embed_induced_graph(graph, [1,2,3,4], Triangular())
        @test patch !== nothing && placed == 37
    end

    @testset "public bounded search" begin
        report = search_unweighted_gadgets(
            SimpleGraph(1), [1], Square(); min_vertices=3, max_vertices=5,
            max_evaluations=50, beam_width=4, mutations_per_candidate=3,
            random_candidates_per_generation=2, rng=MersenneTwister(12),
        )
        @test report.termination_reason == :solution
        @test length(report.trace) == report.evaluated
        gadget = only(report.gadgets)
        @test is_gadget_replacement(
            report.target_graph, gadget.replacement_graph,
            report.target_boundary, gadget.boundary_vertices,
        ) == (true, gadget.constant_offset)
    end
end
