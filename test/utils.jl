@testset "GraphIO" begin
    path = "../data/graphs/graph2.g6"
    graph_dict = readGraphDictFile(path)
    @test readGraphFile(path, 1) == graph_dict["graph1"]
end

@testset "graphPlot" begin
    vertex_num = [3]
    for i in vertex_num
        path = "../data/graphs/graph$i.g6"
        graph_dict = readGraphDictFile(path)
        if !isdir("../data/graphplot/$(i)vertex/")
            mkpath("../data/graphplot/$(i)vertex/")
        end
        plotGraphs(graph_dict, "../data/graphplot/$(i)vertex/")
    end
end

@testset "loadJSON & locate" begin
    path = "../datasets/any_constraint/2bits_any_constraint.json"
    data = loadJSONFile(path)
    @test length(data) == 15
    path = "../datasets/any_constraint/3bits_any_constraint.json"
    data = loadJSONFile(path)
    @test length(data) == 255

    degeneracy = ["001", "010", "100", "111"]
    g_info = findByDegeneracy(data, degeneracy)
    plotColoredGraph(g_info, "../data/graphplot/")
    rm("../data/graphplot/graph.png")
    _, _, ground_states = checkGraphMIS(g_info)
    @test Set(degeneracy) == Set(ground_states)
end