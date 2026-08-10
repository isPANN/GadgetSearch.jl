using GadgetSearch
using Graphs
using Test

function cross_graph()
    graph = SimpleGraph(4)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 2, 4)
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
        @test_throws ArgumentError check_crossing_frame(
            Square(), [(0, 0)], pins, rays,
        )
        @test_throws ArgumentError check_crossing_frame(
            Square(), [square; square[1]], pins, rays,
        )
        @test_throws ArgumentError check_crossing_frame(
            Square(), square, [pins[1], pins[1], pins[3], pins[4]], rays,
        )
        triangular = [(0, 0), (-1, 0), (0, 1), (1, 0), (-1, -1)]
        triangular_pins = triangular[2:5]
        @test all(check_crossing_frame(
            Triangular(), triangular, triangular_pins, rays,
        ))
    end

    @testset "fixed SAT positive control" begin
        coordinates = [
            (0,3),(6,6),(1,3),(2,1),(2,2),(2,3),(3,1),(3,3),
            (3,4),(3,5),(4,1),(4,2),(4,4),(4,5),(5,1),(5,2),
            (5,3),(5,4),(5,5),(6,2),(6,3),(6,5),(7,4),
        ]
        frame = (
            pins=[(5,1),(0,3),(3,5),(6,6)],
            rays=[6,3,3,2],
            allowed=coordinates,
        )
        target_reduced = vec(calculate_reduced_alpha_tensor(
            cross_graph(), [1,2,3,4],
        ))
        analysis = GadgetSearch._solve_fixed_crossing_sat(
            target_reduced, Triangular(), frame, 23, 7,
        )
        @test analysis !== nothing
        @test analysis.offset == 7
        @test nv(analysis.graph) == 23
        @test is_gadget_replacement(
            cross_graph(), analysis.graph, [1,2,3,4], analysis.boundary,
        ) == (true, 7.0)
        @test all(check_crossing_frame(
            Triangular(), analysis.patch.coordinates, analysis.patch.pins,
            GadgetSearch._patch_ray_directions(Triangular(), analysis.patch),
        ))
        shifted_target = map(value -> isfinite(value) ? value + 8 : value, target_reduced)
        negative_offset = GadgetSearch._solve_fixed_crossing_sat(
            shifted_target, Triangular(), frame, 23, -1,
        )
        @test negative_offset !== nothing
        @test negative_offset.offset == -1

        ksg_coordinates = [(0, 0), (-1, 0), (0, 1), (1, 0), (0, -1)]
        ksg_frame = (
            pins=ksg_coordinates[2:5],
            rays=[4, 7, 5, 2],
            allowed=ksg_coordinates,
        )
        ksg_patch = GadgetSearch._LatticePatch(
            ksg_coordinates, ksg_frame.pins, ksg_frame.rays,
        )
        ksg_graph, ksg_boundary, _ = GadgetSearch._materialize_lattice_patch(
            Square(), ksg_patch,
        )
        ksg_target = vec(calculate_reduced_alpha_tensor(ksg_graph, ksg_boundary))
        ksg_analysis = GadgetSearch._solve_fixed_crossing_sat(
            ksg_target, Square(), ksg_frame, 5, 0,
        )
        @test ksg_analysis !== nothing
        @test nv(ksg_analysis.graph) == 5
        @test all(check_crossing_frame(
            Square(), ksg_analysis.patch.coordinates, ksg_analysis.patch.pins,
            GadgetSearch._patch_ray_directions(Square(), ksg_analysis.patch),
        ))
        @test isnothing(GadgetSearch._solve_fixed_crossing_sat(
            ksg_target, Square(), ksg_frame, 6, 0,
        ))
    end

    @testset "public bounded search" begin
        @test_throws ArgumentError search_unweighted_gadgets(
            path_graph(3), [1,2,3], Triangular(),
        )
        @test_throws ArgumentError search_unweighted_gadgets(
            path_graph(4), [1,2,3,4]; window_side=1,
        )
        ksg = search_unweighted_gadgets(
            cross_graph(), [1,2,3,4], Square();
            min_vertices=4, max_vertices=4, max_evaluations=1,
        )
        @test ksg.lattice == :KSG
        @test ksg.evaluated <= 1
        @test ksg.termination_reason in (:budget, :frame_budget)

        exhausted = search_unweighted_gadgets(
            cross_graph(), [1,2,3,4], Triangular();
            min_vertices=4, max_vertices=4, max_evaluations=1,
        )
        @test exhausted.evaluated == 1
        @test isempty(exhausted.gadgets)
        @test exhausted.termination_reason in (:budget, :frame_budget)
    end
end
