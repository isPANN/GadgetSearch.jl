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

@testset "Unified Unweighted Search" begin
    @testset "Single Target Filter" begin
        cross = _cross_graph()
        batoidea = _batoidea_graph()

        filter_fn = make_unweighted_filter(cross, [1, 2, 3, 4])

        replacement = filter_fn(batoidea, nothing, nothing)
        @test replacement !== nothing
        @test replacement isa UnweightedGadget
        @test replacement.boundary_vertices == [1, 2, 3, 4]
        @test replacement.constant_offset == 2.0
        @test replacement.target_index == 1

        self_replacement = filter_fn(cross, nothing, [1, 2, 3, 4])
        @test self_replacement !== nothing
        @test self_replacement.constant_offset == 0.0
        @test self_replacement.target_index == 1

        @test filter_fn(SimpleGraph(3), nothing, nothing) === nothing
    end

    @testset "Multi-Target Filter" begin
        cross = _cross_graph()
        batoidea = _batoidea_graph()

        targets = [
            (cross, [1, 3, 2, 4]),
            (cross, [1, 2, 4, 3]),
            (cross, [1, 2, 3, 4]),
        ]
        filter_fn = make_unweighted_filter(targets)
        result = filter_fn(batoidea, nothing, [1, 2, 3, 4])

        @test result !== nothing
        @test result isa UnweightedGadget
        @test result.pattern_graph == cross
        @test result.target_index == 3
        @test result.constant_offset == 2.0
    end

    @testset "Multi-Target Prefilter" begin
        loader = GraphLoader(GraphDataset([_to_g6(_cross_graph())]), pinset=[1, 3])
        targets = [(_edge_graph(), [1, 2])]

        results_prefilter_on = search_unweighted_gadgets(targets, loader; prefilter=true)
        results_prefilter_off = search_unweighted_gadgets(targets, loader; prefilter=false)

        @test isempty(results_prefilter_on)
        @test length(results_prefilter_off) == 1
        @test results_prefilter_off[1].target_index == 1
    end

    @testset "Search Keywords" begin
        cross = _cross_graph()
        batoidea = _batoidea_graph()
        loader = GraphLoader(
            GraphDataset([_to_g6(cross), _to_g6(batoidea)]),
            pinset=[1, 2, 3, 4],
        )

        limited = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader; limit=1)
        @test length(limited) == 1
        @test limited[1].constant_offset == 0.0
        @test limited[1].target_index == 1

        capped = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader; max_results=1)
        @test length(capped) == 1
    end

    @testset "Single and Multi-Target Search" begin
        cross = _cross_graph()
        batoidea = _batoidea_graph()
        loader = GraphLoader(
            GraphDataset([_to_g6(cross), _to_g6(batoidea)]),
            pinset=[1, 2, 3, 4],
        )

        single_results = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader)
        @test single_results isa Vector{UnweightedGadget}
        @test any(r -> r.constant_offset == 0.0, single_results)
        @test any(r -> r.constant_offset == 2.0, single_results)
        @test all(r -> r.target_index == 1, single_results)

        multi_results = search_unweighted_gadgets(
            [
                (cross, [1, 3, 2, 4]),
                (cross, [1, 2, 4, 3]),
                (cross, [1, 2, 3, 4]),
            ],
            loader;
            max_results=2,
        )
        @test multi_results isa Vector{UnweightedGadget}
        @test length(multi_results) == 2
        @test any(r -> r.target_index == 3 && r.constant_offset == 0.0, multi_results)
        @test any(r -> r.target_index == 3 && r.constant_offset == 2.0, multi_results)
    end

    @testset "inf_mask" begin
        @test inf_mask([0.0, -Inf, 3.0, -Inf]) == BigInt(10)
        @test inf_mask(fill(-Inf, 4)) == BigInt(15)

        reduced = calculate_reduced_alpha_tensor(_cross_graph(), [1, 2, 3, 4])
        @test inf_mask(reduced) == BigInt(60576)
        @test inf_mask(reduced) == inf_mask(content.(reduced))
    end

    @testset "pins_prefilter" begin
        connected = _connected_graph()
        disconnected = _cross_graph()
        isolated = _isolated_graph()

        @test pins_prefilter(connected, [1])
        @test pins_prefilter(disconnected, [1, 2])
        @test !pins_prefilter(disconnected, [1])
        @test !pins_prefilter(isolated, [1])
        @test pins_prefilter(isolated, [1, 3])

        @test_throws ArgumentError pins_prefilter(connected, [1, 1])
        @test_throws ArgumentError pins_prefilter(connected, [0])
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
            @test results[1].target_index == 1
        finally
            isfile(path) && rm(path)
        end
    end
end
