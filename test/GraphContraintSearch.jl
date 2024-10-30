@testset "searchForSingleConstraint" begin
    vertex_num = [3, 4, 5]
    for i in vertex_num
        path = "../data/graphs/graph$i.g6"
        graph_dict = readGraphDictFile(path)
        checkSingleConstraint(graph_dict, 3, [0])
    end
end

@testset "searchForAnyConstraint" begin    
    searchForAnyConstraint([2,3], 2, "../data/graphs/", "../")
    rm("../2bits_any_constraint.json")
    searchForSingleConstraint([2,3], 2, [1,2], "../data/graphs/", "../")
    rm("../2bits_[1, 2].json")
end