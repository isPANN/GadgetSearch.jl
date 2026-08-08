using GadgetSearch
using Graphs
using JSON3
using Random
using Test

function _cross_graph()
    graph = SimpleGraph(4)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 2, 4)
    return graph
end

function _reconstruct_record(lattice, record)
    positions = GadgetSearch.get_physical_positions(lattice, record.lattice_coordinates)
    graph = GadgetSearch.unit_disk_graph(positions, get_radius(lattice))
    indices = Dict(coordinate => index for (index, coordinate) in enumerate(record.lattice_coordinates))
    boundary = [indices[pin] for pin in record.pin_coordinates]
    return graph, boundary
end

@testset "Unweighted Search" begin
    @testset "every state is a concrete induced lattice patch" begin
        for lattice in (Square(), Triangular())
            target = SimpleGraph(1)
            report = search_unweighted_gadgets(
                target,
                [1],
                lattice;
                min_vertices=3,
                max_vertices=5,
                max_evaluations=50,
                beam_width=4,
                mutations_per_candidate=3,
                random_candidates_per_generation=2,
                rng=MersenneTwister(12),
            )

            @test report isa UnweightedSearchResult
            @test report.lattice == (lattice isa Square ? :KSG : :triangular)
            @test report.target_graph == target
            @test report.target_boundary == [1]
            @test !isempty(report.gadgets)
            @test report.evaluated <= 50
            @test length(report.trace) == report.evaluated

            for record in report.trace
                graph, boundary = _reconstruct_record(lattice, record)
                @test nv(graph) == record.vertices
                @test ne(graph) == record.edges
                @test boundary == record.boundary_vertices
                @test is_connected(graph)
            end

            gadget = only(report.gadgets)
            reconstructed = GadgetSearch.unit_disk_graph(gadget.pos, get_radius(lattice))
            @test reconstructed == gadget.replacement_graph
            @test length(gadget.pin_rays) == 1
            @test any(record ->
                record.lattice_coordinates == gadget.lattice_coordinates &&
                record.boundary_vertices == gadget.boundary_vertices,
                report.trace,
            )
            @test is_gadget_replacement(
                target,
                gadget.replacement_graph,
                [1],
                gadget.boundary_vertices,
            ) == (true, gadget.constant_offset)
        end
    end

    @testset "four-pin results require the complete crossing frame" begin
        target = _cross_graph()
        for (lattice, seed, budget) in ((Square(), 2, 400), (Triangular(), 2027, 1_000))
            report = search_unweighted_gadgets(
                target,
                [1, 2, 3, 4],
                lattice;
                min_vertices=5,
                max_vertices=17,
                max_evaluations=budget,
                beam_width=32,
                mutations_per_candidate=8,
                random_candidates_per_generation=8,
                rng=MersenneTwister(seed),
            )

            @test report.evaluated == budget
            lattice isa Square && @test !isempty(report.gadgets)
            @test all(record -> 0 <= record.frame_violations <= 4, report.trace)
            @test all(record -> record.frame_violations == 0, report.trace)
            for gadget in report.gadgets
                checks = check_crossing_frame(
                    lattice,
                    gadget.lattice_coordinates,
                    gadget.lattice_coordinates[gadget.boundary_vertices],
                    gadget.pin_rays,
                )
                @test all(checks)
                @test is_gadget_replacement(
                    target,
                    gadget.replacement_graph,
                    [1, 2, 3, 4],
                    gadget.boundary_vertices,
                ) == (true, gadget.constant_offset)
            end
            @test any(record -> record.parent_key !== nothing, report.trace)
            @test all(record -> record.lattice == report.lattice, report.trace)
        end
    end

    @testset "checks G1-G4 exactly" begin
        square_coordinates = [(0, 0), (-1, 0), (0, 1), (1, 0), (0, -1)]
        square_pins = [(-1, 0), (0, 1), (1, 0), (0, -1)]
        square_rays = [(-1, 0), (0, 1), (1, 0), (0, -1)]
        @test all(check_crossing_frame(Square(), square_coordinates, square_pins, square_rays))

        triangular_coordinates = [(0, 0), (-1, 0), (0, 1), (1, 0), (-1, -1)]
        triangular_pins = [(-1, 0), (0, 1), (1, 0), (-1, -1)]
        triangular_rays = [(-1, 0), (0, 1), (1, 0), (0, -1)]
        @test all(check_crossing_frame(
            Triangular(), triangular_coordinates, triangular_pins, triangular_rays,
        ))

        blocked_coordinates = [square_coordinates; (-2, 1)]
        blocked = check_crossing_frame(Square(), blocked_coordinates, square_pins, square_rays)
        @test !blocked.G4
        @test_throws ArgumentError check_crossing_frame(
            Triangular(), triangular_coordinates, triangular_pins,
            [(2, 0); triangular_rays[2:4]],
        )
    end

    @testset "records self-contained dynamic transitions" begin
        report = search_unweighted_gadgets(
            _cross_graph(),
            [1, 2, 3, 4],
            Triangular();
            min_vertices=5,
            max_vertices=9,
            max_evaluations=80,
            beam_width=6,
            mutations_per_candidate=6,
            random_candidates_per_generation=3,
            max_results=4,
            rng=MersenneTwister(8),
        )

        evaluated_keys = Set(record.key for record in report.trace)
        @test all(
            record.parent_key === nothing || record.parent_key in evaluated_keys
            for record in report.trace
        )

        path = tempname()
        try
            @test save_unweighted_trace(path, report) == path
            rows = JSON3.read.(readlines(path))
            @test length(rows) == report.evaluated
            @test rows[1].target_vertices == nv(report.target_graph)
            @test Tuple.(rows[1].target_edges) == [(1, 3), (2, 4)]
            @test rows[1].lattice == "triangular"
            @test Tuple.(rows[1].lattice_coordinates) == report.trace[1].lattice_coordinates
            @test Tuple.(rows[1].pin_coordinates) == report.trace[1].pin_coordinates
            @test Tuple.(rows[1].pin_rays) == report.trace[1].pin_rays
        finally
            isfile(path) && rm(path)
        end
    end

    @testset "validates the explicit search budget" begin
        target = SimpleGraph(1)
        @test_throws ArgumentError search_unweighted_gadgets(
            target,
            [1],
            Square();
            min_vertices=2,
            max_vertices=1,
        )
        @test_throws ArgumentError search_unweighted_gadgets(
            target,
            [1],
            Square();
            max_evaluations=0,
        )
        @test_throws ArgumentError search_unweighted_gadgets(
            target,
            [1],
            Square();
            exploration_fraction=1.0,
        )
    end

    @testset "verifier behavior remains covered" begin
        @test GadgetSearch.inf_mask([0.0, -Inf, 3.0, -Inf]) == BigInt(10)
        @test GadgetSearch.inf_mask(fill(-Inf, 4)) == BigInt(15)
        reduced = calculate_reduced_alpha_tensor(_cross_graph(), [1, 2, 3, 4])
        @test GadgetSearch.inf_mask(reduced) == BigInt(60576)
    end
end
