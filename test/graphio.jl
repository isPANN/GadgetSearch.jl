@testset "save and recover 1" begin
    using Graphs

    g1 = SimpleGraph(3)
    add_edge!(g1, 1, 2)
    add_edge!(g1, 2, 3)
    add_edge!(g1, 1, 3)

    g2 = SimpleGraph(4)
    add_edge!(g2, 1, 2)
    add_edge!(g2, 2, 3)
    
    g3 = SimpleGraph(5)
    add_edge!(g3, 1, 2)
    add_edge!(g3, 2, 3)
    add_edge!(g3, 3, 4)
    add_edge!(g3, 4, 5)
    add_edge!(g3, 5, 1)

    graphs = Dict("graph1" => g1, "graph2" => g2, "graph3" => g3)

    save_g6_graph(g1, pkgdir(GadgetSearch, "test", "graph.g6"))
    g1_2 = read_g6_graph(pkgdir(GadgetSearch, "test", "graph.g6"), 1)
    @test g1 == g1_2

    save_g6_graphs(graphs, pkgdir(GadgetSearch, "test", "graphs.g6"))
    graphs_2 = read_g6_graphs(pkgdir(GadgetSearch, "test", "graphs.g6"))
    @test graphs == graphs_2

    for (key, graph) in graphs_2
        gid = get_gids(key)
        @test graph == graphs["graph$gid"]
    end
    # rm(pkgdir(GadgetSearch, "test", "graph.g6"))
    # rm(pkgdir(GadgetSearch, "test", "graphs.g6"))
end

@testset "save and recover 2" begin
    graphs = read_g6_graphs(pkgdir(GadgetSearch, "data", "graphs", "graph5.g6"))
    @test length(graphs) == 34

    save_g6_graphs(graphs, pkgdir(GadgetSearch, "test", "graphs.g6"))
    graphs_2 = read_g6_graphs(pkgdir(GadgetSearch, "test", "graphs.g6"))
    @test graphs == graphs_2
    
    # rm(pkgdir(GadgetSearch, "test", "graphs.g6"))
end