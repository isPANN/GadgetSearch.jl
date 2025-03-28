@testset "GraphIO" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph2.g6")
    graph_dict = read_graph_dict(path)
    @test read_graph(path, 1) == graph_dict["graph1"]
end

@testset "format_truth_table" begin
    truth_table1 = [0 0 1; 0 1 0; 1 0 0; 1 1 1]
    truth_table2 = BitMatrix([0 0 1; 0 1 0; 1 0 0; 1 1 1])
    @test format_truth_table(truth_table1) == format_truth_table(truth_table2)
end

@testset "gate_info" begin
    grstates = generic_rule(110, [3, 1]; show_info=true)
    @test Set(grstates) == Set([14, 13, 11, 8, 7, 5, 3, 0])
    @test 110 == reconstruct_rule_id(grstates, [3, 1])

    grstates = generic_rule(5, 2; show_info=true)
    @test Set(grstates) == Set([0, 2])
    @test 5 == reconstruct_rule_id(grstates)
end

@testset "split large file" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph3.g6")
    search_single_rule(path, 2; truth_table=[0 1; 1 0], max_file_size_mb=0, split_size=3)
    @test isdir(pkgdir(GadgetSearch, "data", "graphs", "tmp_graph3.g6_3"))
    rm(pkgdir(GadgetSearch, "data", "graphs", "tmp_graph3.g6_3"), recursive=true)
end

@testset "generate degeneracy cases" begin
    generate_degeneracy_cases(2)
    generate_degeneracy_cases(3; truth_table=true)
end