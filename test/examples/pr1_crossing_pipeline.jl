using Test
using GadgetSearch
using Graphs

include(joinpath(pkgdir(GadgetSearch), "examples", "pr1_crossing_pipeline.jl"))

@testset "PR1 crossing pipeline example" begin
    @testset "build_crossing_demo_candidates returns crossing variants" begin
        target_graph = build_cross_graph()
        candidates = build_crossing_demo_candidates(target_graph)

        @test length(candidates) == 3
        @test candidates[1] == target_graph
        @test nv(candidates[2]) == 5
        @test nv(candidates[3]) == 6
    end

    @testset "build_loader_from_candidates validates non-empty input" begin
        @test_throws ArgumentError build_loader_from_candidates(SimpleGraph{Int}[], [1, 2, 3, 4])
    end

    @testset "run_search_demo searches provided candidates" begin
        target_graph = build_cross_graph()
        target_boundary = [1, 2, 3, 4]
        candidates = build_crossing_demo_candidates(target_graph)
        results = run_search_demo(target_graph, target_boundary, candidates)

        @test !isempty(results)
        @test all(r -> r.pattern_graph == target_graph, results)
    end
end
