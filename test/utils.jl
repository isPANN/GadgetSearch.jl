@testset "GraphIO" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph2.g6")
    graph_dict = readgraphdict(path)
    @test readgraph(path, 1) == graph_dict["graph1"]

    graph_vector = read_g6_file(path)
    @test length(graph_vector) == 2
    @test graph_vector[1] == graph_dict["graph1"]
end

@testset "graphPlot" begin
    vertex_num = [3]
    for i in vertex_num
        path = pkgdir(GadgetSearch, "data", "graphs", "graph$i.g6")
        graph_dict = readgraphdict(path)
        if !isdir(pkgdir(GadgetSearch, "data", "graphplot", "$(i)vertex/"))
            mkpath(pkgdir(GadgetSearch, "data", "graphplot", "$(i)vertex/"))
        end
        plotgraphs(graph_dict, :png; saved_path=pkgdir(GadgetSearch, "data", "graphplot", "$(i)vertex/"))
    end
end

@testset "loadJSON & locate" begin
    path = pkgdir(GadgetSearch, "datasets", "logic_gates", "2in2out.json")
    data = loadjsonfile(path)
    @test length(data) == 256
    degeneracy = ["0011", "0101", "1010", "1100"]
    g_info = find_by_degeneracy(data, degeneracy)

    g_info1 = find_by_degeneracy(path, degeneracy)
    @test g_info == g_info1

    g_info2 = find_by_degeneracy(path, 4, [3, 5, 10, 12])
    @test g_info == g_info2

    plotcoloredgraph(g_info, :png; saved_path = pkgdir(GadgetSearch, "data", "graphplot"), name = "graph")
    rm(pkgdir(GadgetSearch, "data", "graphplot", "graph.png"))
    _, _, ground_states = checkgraphmis(g_info)
    @test Set(degeneracy) == Set(ground_states)

    path = pkgdir(GadgetSearch, "datasets", "any_constraint", "2bits_any_constraint.json")
    data = loadjsonfile(path)
    @test length(data) == 15
    path = pkgdir(GadgetSearch, "datasets", "any_constraint", "3bits_any_constraint.json")
    data = loadjsonfile(path)
    @test length(data) == 255
    
    degeneracy = ["001", "010", "100", "111"]
    g_info = find_by_degeneracy(data, degeneracy)
    plotcoloredgraph(g_info, :png; saved_path = pkgdir(GadgetSearch, "data", "graphplot"), name = "graph")
    rm(pkgdir(GadgetSearch, "data", "graphplot", "graph.png"))
    _, _, ground_states = checkgraphmis(g_info)
    @test Set(degeneracy) == Set(ground_states)

    data = loadjsonfile(path, :gate_id)
    g_info = find_by_gateid(data, 1)
    plotcoloredgraphs(data, :png; saved_path = pkgdir(GadgetSearch, "data", "graphplot", "test_plot"), name = "graph")
    rm(pkgdir(GadgetSearch, "data", "graphplot", "test_plot"), force = true, recursive = true)

    g_info2 = find_by_gateid(path, 1)
    @test g_info == g_info2
end

@testset "genericGate" begin
    genericgate(110, 3, 1)
    showgateinfo(110, 3, 1)
end
