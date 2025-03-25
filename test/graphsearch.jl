@testset "search_single_gate" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph6.g6")
    gate_id = 0
    ground_states = generic_gate(0,2,2)
    truth_table = [0 0 0 0; 0 1 0 0; 1 0 0 0; 1 1 0 0]

    result1 = search_single_rule(path, [2,2], ground_states);
    result2 = search_single_rule(path, [2,2], truth_table);
    result3 = search_single_rule(path, [2,2], gate_id);
    @test result1.graph_id == result2.graph_id == result3.graph_id

    save_path = pkgdir(GadgetSearch, "test", "test.json")
    flag = save_results_to_json([result1, result2], save_path)
    @test flag == save_path
    rm(save_path)
end