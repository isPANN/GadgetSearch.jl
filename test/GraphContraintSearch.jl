@testset "searchForSingleConstraint" begin
    vertex_num = [3, 4, 5]
    for i in vertex_num
        path = "../data/graphs/graph$i.g6"
        graph_dict = readGraphDictFile(path)
        checkSingleConstraint(graph_dict, 3, [0])
    end
end

@testset "searchForAnyConstraint" begin    
    fn = searchForAnyConstraint([2,3], 2, "../data/graphs/", "../")
    rm("$fn")
    fn = searchForSingleConstraint([2,3], 2, [1,2], "../data/graphs/", "../")
    rm("$fn")
    fn = searchForSingleConstraint([6],3,["111","101","011","000"],"../data/graphs/", "../")
    rm("$fn")
end

@testset "searchForLogicGates" begin
    fn = searchForSingleGate([6], 2, 2, 2, "../data/graphs/", "../")
    rm("$fn")
    fn = searchForGates([4], 2, 2, "../data/graphs/", "../")
    rm("$fn")
end