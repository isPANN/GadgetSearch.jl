@testset "search single gate using path" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph6.g6")
    rule_id = 0
    ground_states = generic_rule(0,[2,2])
    truth_table = [0 0 0 0; 0 1 0 0; 1 0 0 0; 1 1 0 0]

    result1 = search_single_rule(path, [2,2]; ground_states=ground_states);
    result2 = search_single_rule(path, [2,2]; truth_table=truth_table);
    result3 = search_single_rule(path, [2,2]; rule_id=rule_id);
    @test result1.graph_id == result2.graph_id == result3.graph_id

    save_path = pkgdir(GadgetSearch, "test", "test.json")
    save_results_to_json([result1, result2], save_path)
    @test isfile(save_path)
    rm(save_path)
end

@testset "search single gate using graph" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph6.g6")
    rule_id = 0
    ground_states = generic_rule(0,[2,2])
    truth_table = [0 0 0 0; 0 1 0 0; 1 0 0 0; 1 1 0 0]

    graph_dict = read_graph_dict(path)
    result1 = search_single_rule(graph_dict, [2,2]; ground_states=ground_states);
    result2 = search_single_rule(graph_dict, [2,2]; truth_table=truth_table);
    result3 = search_single_rule(graph_dict, [2,2]; rule_id=rule_id);

    graph = read_graph(path, 122)
    result4 = search_single_rule(graph, [2,2]; ground_states=ground_states);
    result5 = search_single_rule(graph, [2,2]; truth_table=truth_table);
    result6 = search_single_rule(graph, [2,2]; rule_id=rule_id);

    @test result1.graph_id == result2.graph_id == result3.graph_id
    @test result4.graph_id == result5.graph_id == result6.graph_id

    # test greedy search
    result1 = search_single_rule(graph_dict, [2,2]; ground_states=ground_states, greedy=true);
    result2 = search_single_rule(graph_dict, [2,2]; ground_states=ground_states, threshold=1, max_samples=1);
end

@testset "search generic constraints using graph path" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph5.g6")
    ground_states = [1, 2]
    truth_table = [0 1; 1 0]
    rule_id = 6

    result1 = search_single_rule(path, 2; ground_states=ground_states);
    result2 = search_single_rule(path, 2; truth_table=truth_table);
    result3 = search_single_rule(path, 2; rule_id=rule_id);
    @test result1.graph_id == result2.graph_id == result3.graph_id
end

@testset "search with different optimizers" begin
    path = pkgdir(GadgetSearch, "data", "graphs", "graph5.g6")
    ground_states = [1, 2]

    # Test with default HiGHS optimizer
    result1 = search_single_rule(path, 2; ground_states=ground_states);

    # Test with GLPK optimizer
    result2 = search_single_rule(path, 2; ground_states=ground_states, optimizer=GLPK.Optimizer);

    # Both should find a solution (though possibly different ones)
    @test !isnothing(result1)
    @test !isnothing(result2)
end