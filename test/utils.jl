@testset "GraphIO" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph2.g6")
    graph_dict = readGraphDictFile(path)
    @test readGraphFile(path, 1) == graph_dict["graph1"]
end

@testset "graphPlot" begin
    vertex_num = [3]
    for i in vertex_num
        path = pkgdir(GadgetSearch, "data", "graphs", "graph$i.g6")
        graph_dict = readGraphDictFile(path)
        if !isdir(pkgdir(GadgetSearch, "data", "graphplot", "$(i)vertex/"))
            mkpath(pkgdir(GadgetSearch, "data", "graphplot", "$(i)vertex/"))
        end
        plotGraphs(graph_dict, pkgdir(GadgetSearch, "data", "graphplot", "$(i)vertex/"))
    end
end

@testset "loadJSON & locate" begin
    path = pkgdir(GadgetSearch, "datasets", "logic_gates", "2in_2out.json")
    data = loadJSONFile(path)
    @test length(data) == 256
    degeneracy = ["0011", "0101", "1010", "1100"]
    g_info = findByDegeneracy(data, degeneracy)
    plotColoredGraph(g_info, pkgdir(GadgetSearch, "data", "graphplot"))
    rm(pkgdir(GadgetSearch, "data", "graphplot", "graph.png"))
    _, _, ground_states = checkGraphMIS(g_info)
    @test Set(degeneracy) == Set(ground_states)

    path = pkgdir(GadgetSearch, "datasets", "any_constraint", "2bits_any_constraint.json")
    data = loadJSONFile(path)
    @test length(data) == 15
    path = pkgdir(GadgetSearch, "datasets", "any_constraint", "3bits_any_constraint.json")
    data = loadJSONFile(path)
    @test length(data) == 255
    
    degeneracy = ["001", "010", "100", "111"]
    g_info = findByDegeneracy(data, degeneracy)
    plotColoredGraph(g_info, pkgdir(GadgetSearch, "data", "graphplot"))
    rm(pkgdir(GadgetSearch, "data", "graphplot", "graph.png"))
    _, _, ground_states = checkGraphMIS(g_info)
    @test Set(degeneracy) == Set(ground_states)
end

@testset "genericGate" begin
    genericGate(110, 3, 1)
    showGateInfo(110, 3, 1)
end