@testset "graphIO" begin
    vertex_num = [3]
    for i in vertex_num
        path = "../data/graphs/graph$i.g6"
        graph_dict = readGraphDictFile(path)
        if !isdir("../data/graphplot/$(i)vertex/")
            mkpath("../data/graphplot/$(i)vertex/")
        end
        plotGraph(graph_dict, "../data/graphplot/$(i)vertex/")
    end
end


@testset "MIS" begin
    vertex_num = [3, 4, 5]
    for i in vertex_num
        path = "../data/graphs/graph$i.g6"
        graph_dict = readGraphDictFile(path)
        checkSingleConstraint(graph_dict, 3, [0])
    end
    
    searchForAnyConstraint([2,3], 2, "../data/graphs/", "../")
    rm("../2bits_any_constraint.json")
    searchForSingleConstraint([2,3], 2, [1,2], "../data/graphs/", "../")
    rm("../2bits_[1, 2].json")
end