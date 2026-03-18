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

function _subdivide_edge(g::SimpleGraph{Int}, u::Int, v::Int, count::Int=1)
    count >= 0 || throw(ArgumentError("count must be non-negative"))
    has_edge(g, u, v) || throw(ArgumentError("edge ($u, $v) does not exist"))
    count == 0 && return copy(g)

    expanded = copy(g)
    rem_edge!(expanded, u, v)
    previous = u
    for _ in 1:count
        add_vertex!(expanded)
        next_vertex = nv(expanded)
        add_edge!(expanded, previous, next_vertex)
        previous = next_vertex
    end
    add_edge!(expanded, previous, v)
    return expanded
end

function _subdivide_edges(g::SimpleGraph{Int}, specs::Vector{Tuple{Int, Int, Int}})
    expanded = copy(g)
    for (u, v, count) in specs
        expanded = _subdivide_edge(expanded, u, v, count)
    end
    return expanded
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
        @test all(r -> r[2] == [1, 2, 3, 4], reprs)

        reduced_tensors = [
            vec(Float64.(content.(calculate_reduced_alpha_tensor(r[1], r[2]))))
            for r in reprs
        ]
        @test allunique(reduced_tensors)
    end

    @testset "equivalent_representations: optional boundary permutations" begin
        cross = _cross_graph()
        reprs_fixed = equivalent_representations(cross, [1, 2, 3, 4])
        reprs_relaxed = equivalent_representations(
            cross,
            [1, 2, 3, 4];
            preserve_boundary_roles=false,
        )

        @test length(reprs_relaxed) > length(reprs_fixed)
        @test any(r -> r[2] != [1, 2, 3, 4], reprs_relaxed)
        @test all(r -> r[2] == [1, 2, 3, 4], reprs_fixed)
    end

    @testset "equivalent_representations: edge subdivision expansion" begin
        cross = _cross_graph()
        reprs = equivalent_representations(cross, [1, 2, 3, 4]; max_added_vertices=1)

        @test any(r -> nv(r[1]) == 5 && ne(r[1]) == 3 && r[2] == [1, 2, 3, 4], reprs)
    end

    @testset "equivalent_representations: multiple subdivisions" begin
        cross = _cross_graph()

        reprs_two = equivalent_representations(cross, [1, 2, 3, 4]; max_added_vertices=2)
        @test any(r -> nv(r[1]) == 6 && ne(r[1]) == 4, reprs_two)

        reprs_four = equivalent_representations(cross, [1, 2, 3, 4]; max_added_vertices=4)
        @test any(r -> nv(r[1]) == 8, reprs_four)
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

    @testset "search_unweighted_gadgets: explicit max_added_vertices" begin
        cross = _cross_graph()
        expanded_cross = _subdivide_edges(cross, [(1, 3, 1), (2, 4, 1)])
        loader = GraphLoader(
            GraphDataset([_to_g6(expanded_cross)]),
            pinset=[1, 2, 3, 4],
        )

        @test isempty(search_unweighted_gadgets(cross, [1, 2, 3, 4], loader))

        results = search_unweighted_gadgets(
            cross,
            [1, 2, 3, 4],
            loader;
            max_added_vertices=2,
        )
        @test length(results) == 1
    end

    @testset "search_unweighted_gadgets: prefilter rejects disconnected pin coverage under subdivisions" begin
        loader = GraphLoader(GraphDataset([_to_g6(_cross_graph())]), pinset=[1, 3])
        edge = _edge_graph()

        results_on = search_unweighted_gadgets(edge, [1, 2], loader; prefilter=true)
        results_off = search_unweighted_gadgets(edge, [1, 2], loader; prefilter=false)

        @test isempty(results_on)
        @test length(results_off) == 1
    end

    @testset "search_unweighted_gadgets: subdivision candidate" begin
        cross = _cross_graph()
        expanded_cross = _subdivide_edges(cross, [(1, 3, 1), (2, 4, 1)])
        loader = GraphLoader(
            GraphDataset([_to_g6(expanded_cross)]),
            pinset=[1, 2, 3, 4],
        )

        results = search_unweighted_gadgets(
            cross,
            [1, 2, 3, 4],
            loader;
            max_added_vertices=2,
        )

        @test length(results) == 1
        @test results[1].replacement_graph == expanded_cross
        @test results[1].boundary_vertices == [1, 2, 3, 4]
    end

    @testset "logical flip target variants" begin
        cross = _cross_graph()
        reprs = equivalent_representations(cross, [1, 2, 3, 4])
        target_data_no = GadgetSearch._build_target_data(reprs; include_logical_flips=false)
        target_data_yes = GadgetSearch._build_target_data(reprs; include_logical_flips=true)

        @test length(target_data_yes) > length(target_data_no)
        @test any(td -> td.flip_mask == [1], target_data_yes)

        base_tensor = Float64.(content.(calculate_reduced_alpha_tensor(cross, [1, 2, 3, 4])))
        flipped_tensor = vec(GadgetSearch.apply_flip_to_tensor(base_tensor, [1]))
        @test any(td -> td.flip_mask == [1] && td.reduced == flipped_tensor, target_data_yes)
        @test any(td -> isempty(td.flip_mask) && td.offset_from_pattern === 0.0, target_data_yes)
        @test all(td -> !isempty(td.flip_mask) ? td.offset_from_pattern === nothing : true, target_data_yes)
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
