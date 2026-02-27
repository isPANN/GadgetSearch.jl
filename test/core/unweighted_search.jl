using GadgetSearch
using Graphs
using GraphIO
using Test

# Build CROSS graph: 4 vertices, edges 1-3 and 2-4
function _cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3); add_edge!(g, 2, 4)
    return g
end

# Build BATOIDEA graph: 11-vertex replacement for CROSS (Figure 6 in the paper)
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

# Encode a graph to a Graph6 string (strips the >>graph6<< header)
function _to_g6(g)
    return GraphIO.Graph6._graphToG6String(g)[11:end]
end

@testset "make_unweighted_filter: CROSS → BATOIDEA" begin
    cross = _cross_graph()
    batoidea = _batoidea_graph()

    # Build a loader containing only BATOIDEA
    loader = GraphLoader(GraphDataset([_to_g6(batoidea)]))

    filter_fn = make_unweighted_filter(cross, [1, 2, 3, 4])
    result = filter_fn(batoidea, nothing, nothing)

    @test result !== nothing
    @test result isa UnweightedGadget
    @test result.boundary_vertices == [1, 2, 3, 4]
    @test result.constant_offset == 2.0
end

@testset "search_unweighted_gadgets: finds BATOIDEA as replacement for CROSS" begin
    cross = _cross_graph()
    batoidea = _batoidea_graph()

    # Loader with two graphs: CROSS itself (offset 0) and BATOIDEA (offset +2)
    loader = GraphLoader(GraphDataset([_to_g6(cross), _to_g6(batoidea)]),
                         pinset=[1, 2, 3, 4])

    results = search_unweighted_gadgets(cross, [1, 2, 3, 4], loader)

    @test length(results) >= 1
    # BATOIDEA should appear with constant_c == 2.0
    @test any(r -> r.constant_offset == 2.0, results)
end

