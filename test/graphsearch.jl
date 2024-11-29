@testset "searchForSingleConstraint" begin
    path = pkgdir(GadgetSearch, "data", "graphs")
    result = check_singleconstraint(2, 2, ["01","10"], path)
    save_results_to_json([result], pkgdir(GadgetSearch, "test.json"))

    g_read_from_json = find_by_gateid(pkgdir(GadgetSearch, "test.json"), 0)
    _, _, ground_states = checkgraphmis(g_read_from_json)
    @test Set(["01","10"]) == Set(ground_states)

    result2 = check_singleconstraint(2, 2, [1, 2], path)
end

@testset "searchForLogicGates" begin
    path = pkgdir(GadgetSearch, "data", "graphs")
    check_singleconstraint(9, [3,1], 0, path)
end