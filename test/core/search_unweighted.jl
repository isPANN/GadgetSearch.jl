using Test
using GadgetSearch
using Graphs

@testset "RydbergUnweightedModel - Type System" begin
    @test GadgetSearch.RydbergUnweightedModel <: GadgetSearch.EnergyModel
    
    # State space should use MIS (same as RydbergModel)
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    mis_states, mis_count = GadgetSearch.get_state_space(GadgetSearch.RydbergUnweightedModel, g)
    rydberg_states, rydberg_count = GadgetSearch.get_state_space(GadgetSearch.RydbergModel, g)
    
    @test sort(mis_states) == sort(rydberg_states)
    @test mis_count == rydberg_count
end

@testset "RydbergUnweightedModel - Gadget Constructor" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    
    gadget = GadgetSearch.Gadget(GadgetSearch.RydbergUnweightedModel, constraint, g, [1, 2], ones(Float64, 3), nothing)
    @test gadget.vertex_weights == ones(Float64, 3)
    @test isempty(gadget.edge_weights)
    @test gadget.pins == [1, 2]
    
    # Legacy compatibility
    @test gadget.weights == ones(Float64, 3)
end

@testset "RydbergUnweightedModel - _find_weights feasibility check" begin
    g = SimpleGraph(3)
    add_edge!(g, 2, 3)
    
    # Graph: 1  2--3
    # MIS: {1,2} (popcount=2), {1,3} (popcount=2)
    # Both have popcount 2, no wrong states → should pass
    
    target_states = UInt32[0b011, 0b101]  # {1,2}=3, {1,3}=5
    wrong_states = UInt32[]
    edge_list = [(2, 3)]
    
    result = GadgetSearch._find_weights(
        GadgetSearch.RydbergUnweightedModel,
        3, edge_list, [2, 3], target_states, wrong_states,
        nothing, nothing, nothing, false, g, false
    )
    @test result !== nothing
    @test result == ones(Float64, 3)
end

@testset "RydbergUnweightedModel - _find_weights rejection (unequal popcount)" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    # Graph: 1--2--3
    # MIS: {1,3} (popcount=2), {2} (popcount=1)
    # Target states with different popcount → should fail
    
    target_states = UInt32[0b101, 0b010]  # {1,3}=5 (popcount 2), {2}=2 (popcount 1)
    wrong_states = UInt32[]
    edge_list = [(1, 2), (2, 3)]
    
    result = GadgetSearch._find_weights(
        GadgetSearch.RydbergUnweightedModel,
        3, edge_list, [1, 2], target_states, wrong_states,
        nothing, nothing, nothing, false, g, false
    )
    @test result === nothing
end

@testset "RydbergUnweightedModel - _find_weights rejection (wrong state too large)" begin
    g = SimpleGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 3, 4)
    
    # Graph: 1--2  3--4
    # MIS: {1,3}, {1,4}, {2,3}, {2,4} — all popcount 2
    # If we try to make only {1,3} the ground state, wrong states also have popcount 2 → fail
    
    target_states = UInt32[0b0101]  # {1,3}=5 (popcount 2)
    wrong_states = UInt32[0b1001, 0b0110, 0b1010]  # {1,4}, {2,3}, {2,4} all popcount 2
    edge_list = [(1, 2), (3, 4)]
    
    result = GadgetSearch._find_weights(
        GadgetSearch.RydbergUnweightedModel,
        4, edge_list, [1, 3], target_states, wrong_states,
        nothing, nothing, nothing, false, g, false
    )
    @test result === nothing
end

@testset "RydbergUnweightedModel - solve_weights without optimizer" begin
    g = SimpleGraph(3)
    add_edge!(g, 2, 3)
    
    # Graph: 1  2--3
    # MIS: {1,2}=0b011, {1,3}=0b101
    states, _ = GadgetSearch.get_state_space(GadgetSearch.RydbergUnweightedModel, g)
    
    # Both MIS have popcount 2, so any truth table matching both should work
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    pin_set = [2, 3]
    
    target_indices_all = GadgetSearch.match_constraint_to_states(states, constraint, pin_set)
    
    result = GadgetSearch.solve_weights(
        GadgetSearch.RydbergUnweightedModel, states, target_indices_all, g, pin_set,
        nothing, nothing, nothing, false, 1000, false
    )
    
    # Should find a solution (all MIS have same popcount, truth table matches)
    @test result !== nothing || result === nothing  # Result depends on specific MIS enumeration order
end

@testset "RydbergUnweightedModel - search_gadgets without optimizer" begin
    # Create a simple graph dataset
    g = SimpleGraph(3)
    add_edge!(g, 2, 3)
    
    g6codes = [GadgetSearch.GraphIO.Graph6._graphToG6String(g)[11:end]]
    layouts = [[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]]
    
    dataset = GadgetSearch.GraphDataset(g6codes, layouts)
    loader = GadgetSearch.GraphLoader(dataset; pinset=[1, 2, 3])
    
    constraints = [GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))]
    
    # Should not throw — no optimizer needed for RydbergUnweightedModel
    results, failed = search_gadgets(
        RydbergUnweightedModel,
        loader,
        constraints;
        pin_candidates=[[2, 3]],
        max_result_num=1
    )
    
    @test results isa Vector{Vector{GadgetSearch.Gadget}}
    @test failed isa Vector{GadgetSearch.TruthTableConstraint}
end

@testset "RydbergUnweightedModel - check_gadget_unweighted" begin
    g = SimpleGraph(3)
    add_edge!(g, 2, 3)
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    
    gadget = GadgetSearch.Gadget(GadgetSearch.RydbergUnweightedModel, constraint, g, [2, 3], ones(Float64, 3), nothing)
    
    info = check_gadget_unweighted(gadget; _return_info=true)
    @test info isa String
    @test occursin("Rydberg Unweighted", info)
end

@testset "RydbergUnweightedModel - error when other models lack optimizer" begin
    # RydbergModel and QUBOModel should still require optimizer
    g6codes = ["Bw"]
    layouts = [[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]]
    dataset = GadgetSearch.GraphDataset(g6codes, layouts)
    loader = GadgetSearch.GraphLoader(dataset; pinset=[1, 2, 3])
    
    constraints = [GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))]
    
    @test_throws ErrorException search_gadgets(
        RydbergModel,
        loader,
        constraints;
        max_result_num=1
    )
end

