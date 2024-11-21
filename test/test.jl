@testset "xx" begin
    vertex_num = [10]
    for i in vertex_num
        path = pkgdir(GadgetSearch, "data", "graphs", "graph$i.g6")
        graph_dict = readgraphdict(path)
        check_single_gate(graph_dict, 3, 1, [14,13,11,8,7,5,3,0])
    end
    # 2356
end