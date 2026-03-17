using GadgetSearch
using Graphs
using GraphIO
using GenericTensorNetworks: content
using Test

function _cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    return g
end

function _batoidea_graph()
    g = SimpleGraph(11)
    add_edge!(g, 1, 5);  add_edge!(g, 1, 9)
    add_edge!(g, 2, 5);  add_edge!(g, 2, 6);  add_edge!(g, 2, 7)
    add_edge!(g, 3, 8)
    add_edge!(g, 4, 9);  add_edge!(g, 4, 10); add_edge!(g, 4, 11)
    add_edge!(g, 5, 6);  add_edge!(g, 5, 9);  add_edge!(g, 5, 10)
    add_edge!(g, 6, 7);  add_edge!(g, 6, 9);  add_edge!(g, 6, 10); add_edge!(g, 6, 11)
    add_edge!(g, 7, 8);  add_edge!(g, 7, 10); add_edge!(g, 7, 11)
    add_edge!(g, 8, 11)
    add_edge!(g, 9, 10)
    add_edge!(g, 10, 11)
    return g
end

function _edge_graph()
    g = SimpleGraph(2)
    add_edge!(g, 1, 2)
    return g
end

function _connected_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 1, 4)
    add_edge!(g, 2, 4)
    add_edge!(g, 3, 4)
    return g
end

function _isolated_graph()
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    return g
end

function _to_g6(g)
    return GraphIO.Graph6._graphToG6String(g)[11:end]
end

@testset "Unweighted Search" begin
    @testset "equivalent_representations" begin
        cross = _cross_graph()
        reprs = equivalent_representations(cross, [1, 2, 3, 4])

        @test reprs isa Vector{Tuple{SimpleGraph{Int}, Vector{Int}}}
        @test !isempty(reprs)
        @test reprs[1] == (cross, [1, 2, 3, 4])

        boundaries = [r[2] for r in reprs]
        @test allunique(boundaries)

        reduced_tensors = [
            vec(Float64.(content.(calculate_reduced_alpha_tensor(r[1], r[2]))))
            for r in reprs
        ]
        @test allunique(reduced_tensors)
        @test all(r -> r[1] === cross, reprs)
    end

    @testset "search_unweighted_gadgets: basic" begin
        cross = _cross_graph()
        batoidea = _batoidea_graph()
        loader = GraphLoader(
            GraphDataset([_to_g6(cross), _to_g6(batoidea)]),
            pinset=[1, 2, 3, 4],
        )

        results = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader)

        @test results isa Vector{UnweightedGadget}
        @test any(r -> r.constant_offset == 0.0, results)
        @test any(r -> r.constant_offset == 2.0, results)
        @test all(r -> r.pattern_graph == cross, results)
        @test !hasproperty(UnweightedGadget, :target_index)
    end

    @testset "search_unweighted_gadgets: limit and max_results" begin
        cross = _cross_graph()
        batoidea = _batoidea_graph()
        loader = GraphLoader(
            GraphDataset([_to_g6(cross), _to_g6(batoidea)]),
            pinset=[1, 2, 3, 4],
        )

        limited = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader; limit=1)
        @test length(limited) == 1
        @test limited[1].constant_offset == 0.0

        capped = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader; max_results=1)
        @test length(capped) == 1
    end

    @testset "search_unweighted_gadgets: prefilter" begin
        loader = GraphLoader(GraphDataset([_to_g6(_cross_graph())]), pinset=[1, 3])
        edge = _edge_graph()

        results_on = search_unweighted_gadgets(edge, [1, 2], loader; prefilter=true)
        results_off = search_unweighted_gadgets(edge, [1, 2], loader; prefilter=false)

        @test isempty(results_on)
        @test length(results_off) == 1
    end

    @testset "UnweightedGadget has no target_index" begin
        @test !(:target_index in fieldnames(UnweightedGadget))
    end

    @testset "inf_mask (internal)" begin
        @test GadgetSearch.inf_mask([0.0, -Inf, 3.0, -Inf]) == BigInt(10)
        @test GadgetSearch.inf_mask(fill(-Inf, 4)) == BigInt(15)

        reduced = calculate_reduced_alpha_tensor(_cross_graph(), [1, 2, 3, 4])
        @test GadgetSearch.inf_mask(reduced) == BigInt(60576)
        @test GadgetSearch.inf_mask(reduced) == GadgetSearch.inf_mask(content.(reduced))
    end

    @testset "pins_prefilter (internal)" begin
        connected = _connected_graph()
        disconnected = _cross_graph()
        isolated = _isolated_graph()

        @test GadgetSearch.pins_prefilter(connected, [1])
        @test GadgetSearch.pins_prefilter(disconnected, [1, 2])
        @test !GadgetSearch.pins_prefilter(disconnected, [1])
        @test !GadgetSearch.pins_prefilter(isolated, [1])
        @test GadgetSearch.pins_prefilter(isolated, [1, 3])

        @test_throws ArgumentError GadgetSearch.pins_prefilter(connected, [1, 1])
        @test_throws ArgumentError GadgetSearch.pins_prefilter(connected, [0])
    end

    @testset "Triangular UDG Integration" begin
        path = tempname() * ".g6"
        try
            generate_full_grid_udg(Triangular(), 1, 1; path=path)
            loader = GraphLoader(path; pinset=[1, 2, 3, 4])
            target = loader[1]

            results = search_unweighted_gadgets(target, [1, 2, 3, 4], loader; limit=1, max_results=1)

            @test length(results) == 1
            @test results[1].constant_offset == 0.0
        finally
            isfile(path) && rm(path)
        end
    end
end
