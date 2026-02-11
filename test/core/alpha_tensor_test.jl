using Test
using GadgetSearch
using Graphs

# ============================================================================
# Core α-Tensor Computation Tests
# ============================================================================

@testset "α-Tensor Computation" begin
    @testset "compute_reduced_alpha_tensor - path graph (1--2--3)" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        # Pin vertex: [2]
        pins = [2]
        α = compute_reduced_alpha_tensor(g, pins)
        
        # Config 0 (pin 2 not selected): interior {1,3}, no edges between them → MIS = {1,3}, size 2
        @test α[UInt32(0)] == 2
        
        # Config 1 (pin 2 selected): interior {1,3}, both blocked by pin 2 → size 0
        @test α[UInt32(1)] == 0
    end
    
    @testset "compute_reduced_alpha_tensor - triangle" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        
        pins = [1]
        α = compute_reduced_alpha_tensor(g, pins)
        
        # Config 0 (pin 1 out): interior {2,3}, edge 2--3 → MIS = {2} or {3}, size 1
        @test α[UInt32(0)] == 1
        
        # Config 1 (pin 1 in): interior {2,3}, both blocked by pin 1 → size 0
        @test α[UInt32(1)] == 0
    end
    
    @testset "compute_reduced_alpha_tensor - disconnected edges" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 3, 4)
        
        pins = [1, 3]
        α = compute_reduced_alpha_tensor(g, pins)
        
        # Config 00: interior {2,4}, no edges → MIS = {2,4}, size 2
        @test α[UInt32(0b00)] == 2
        
        # Config 01 (pin 1 in): interior {2,4}, 2 blocked by 1 → feasible {4}, size 1
        @test α[UInt32(0b01)] == 1
        
        # Config 10 (pin 3 in): interior {2,4}, 4 blocked by 3 → feasible {2}, size 1
        @test α[UInt32(0b10)] == 1
        
        # Config 11 (both in): interior {2,4}, 2 blocked by 1, 4 blocked by 3 → size 0
        @test α[UInt32(0b11)] == 0
    end
    
    @testset "compute_reduced_alpha_tensor - infeasible boundary config" begin
        # Two adjacent pins: 1--2 with pins [1, 2]
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        pins = [1, 2]
        α = compute_reduced_alpha_tensor(g, pins)
        
        # Config 00: interior {3}, size 1
        @test α[UInt32(0b00)] == 1
        
        # Config 01 (pin 1 in): interior {3}, 3 not adjacent to 1 → size 1
        @test α[UInt32(0b01)] == 1
        
        # Config 10 (pin 2 in): interior {3}, 3 adjacent to 2 → blocked → size 0
        @test α[UInt32(0b10)] == 0
        
        # Config 11 (both in): pins 1 and 2 are adjacent → INFEASIBLE
        @test α[UInt32(0b11)] == GadgetSearch.INFEASIBLE_ALPHA
    end
    
    @testset "compute_reduced_alpha_tensor - star graph" begin
        # Star: center=1, leaves=2,3,4
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        
        pins = [1]
        α = compute_reduced_alpha_tensor(g, pins)
        
        # Config 0 (center out): interior {2,3,4}, no edges → MIS = {2,3,4}, size 3
        @test α[UInt32(0)] == 3
        
        # Config 1 (center in): all interior blocked → size 0
        @test α[UInt32(1)] == 0
    end
end

# ============================================================================
# Ground Configuration Tests
# ============================================================================

@testset "find_ground_configs" begin
    @testset "path graph ground configs" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        pins = [2]
        α = compute_reduced_alpha_tensor(g, pins)
        
        ground, max_total = GadgetSearch.find_ground_configs(α, 1)
        
        # Config 0: total = 2 + 0 = 2
        # Config 1: total = 0 + 1 = 1
        # Ground = {0}, max_total = 2
        @test max_total == 2
        @test ground == Set([UInt32(0)])
    end
    
    @testset "disconnected edges ground configs" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 3, 4)
        
        pins = [1, 3]
        α = compute_reduced_alpha_tensor(g, pins)
        
        ground, max_total = GadgetSearch.find_ground_configs(α, 2)
        
        # Config 00: total = 2 + 0 = 2
        # Config 01: total = 1 + 1 = 2
        # Config 10: total = 1 + 1 = 2
        # Config 11: total = 0 + 2 = 2
        # All configs have total 2 → all are ground
        @test max_total == 2
        @test length(ground) == 4
    end
    
    @testset "ground configs with infeasible" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        pins = [1, 2]
        α = compute_reduced_alpha_tensor(g, pins)
        
        ground, max_total = GadgetSearch.find_ground_configs(α, 2)
        
        # Config 00: total = 1 + 0 = 1
        # Config 01: total = 1 + 1 = 2
        # Config 10: total = 0 + 1 = 1
        # Config 11: INFEASIBLE → skipped
        # Ground = {01}, max_total = 2
        @test max_total == 2
        @test ground == Set([UInt32(0b01)])
    end
end

# ============================================================================
# Pin Config Extraction Tests
# ============================================================================

@testset "extract_pin_config" begin
    @testset "basic extraction" begin
        # State 0b1010 = vertices {2, 4} selected, pins = [2, 4]
        state = UInt32(0b1010)
        pins = [2, 4]
        config = extract_pin_config(state, pins)
        @test config == UInt32(0b11)  # Both pins selected
    end
    
    @testset "partial selection" begin
        # State 0b0010 = vertex {2} selected, pins = [2, 4]
        state = UInt32(0b0010)
        pins = [2, 4]
        config = extract_pin_config(state, pins)
        @test config == UInt32(0b01)  # Only pin 2 selected
    end
    
    @testset "no pins selected" begin
        state = UInt32(0b0001)  # vertex 1 only
        pins = [2, 3]
        config = extract_pin_config(state, pins)
        @test config == UInt32(0)
    end
end

# ============================================================================
# Alpha Equivalence Tests
# ============================================================================

@testset "check_alpha_equivalence" begin
    @testset "equivalent tensors - constant difference" begin
        α1 = Dict(UInt32(0) => 2, UInt32(1) => 0)
        α2 = Dict(UInt32(0) => 5, UInt32(1) => 3)
        
        is_equiv, c = check_alpha_equivalence(α1, α2)
        @test is_equiv == true
        @test c == 3
    end
    
    @testset "not equivalent - varying difference" begin
        α1 = Dict(UInt32(0) => 2, UInt32(1) => 0)
        α2 = Dict(UInt32(0) => 5, UInt32(1) => 2)
        
        is_equiv, c = check_alpha_equivalence(α1, α2)
        @test is_equiv == false
        @test c === nothing
    end
    
    @testset "different keys" begin
        α1 = Dict(UInt32(0) => 2, UInt32(1) => 0)
        α2 = Dict(UInt32(0) => 5, UInt32(2) => 3)
        
        is_equiv, c = check_alpha_equivalence(α1, α2)
        @test is_equiv == false
    end
    
    @testset "with infeasible configs" begin
        α1 = Dict(UInt32(0) => 2, UInt32(1) => 0, UInt32(3) => GadgetSearch.INFEASIBLE_ALPHA)
        α2 = Dict(UInt32(0) => 5, UInt32(1) => 3, UInt32(3) => GadgetSearch.INFEASIBLE_ALPHA)
        
        is_equiv, c = check_alpha_equivalence(α1, α2)
        @test is_equiv == true
        @test c == 3
    end
end

# ============================================================================
# Gadget Verification Tests
# ============================================================================

@testset "verify_gadget_via_alpha_tensor" begin
    @testset "valid verification - path graph" begin
        # Path: 1--2--3, pins = [1, 3]
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        pins = [1, 3]
        
        # α-tensor for this graph with pins [1, 3]:
        # Config 00: interior {2}, size 1. Total = 1+0 = 1
        # Config 01 (pin 1 in): interior {2}, blocked by 1 → size 0. Total = 0+1 = 1
        # Config 10 (pin 3 in): interior {2}, blocked by 3 → size 0. Total = 0+1 = 1
        # Config 11 (both in): interior {2}, blocked → size 0. Total = 0+2 = 2
        # Ground config = {11}, max_total = 2
        
        # Target: state 0b101 = {1, 3} → pin config = 0b11
        target_states = UInt32[0b101]
        
        result = verify_gadget_via_alpha_tensor(g, pins, target_states)
        @test result !== nothing
        if result !== nothing
            weights, max_total = result
            @test weights == ones(Float64, 3)
            @test max_total == 2
        end
    end
    
    @testset "invalid verification - wrong target" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        pins = [1, 3]
        
        # Ground config is {11} (both pins in), but target says {10} (only pin 3)
        # State 0b100 = {3} → pin config = 0b10
        target_states = UInt32[0b100]
        
        result = verify_gadget_via_alpha_tensor(g, pins, target_states)
        @test result === nothing
    end
end

# ============================================================================
# Type System and Integration Tests
# ============================================================================

@testset "RydbergUnweightedModel Type System" begin
    @test RydbergUnweightedModel <: GadgetSearch.EnergyModel
    
    # State space should use MIS (same as RydbergModel)
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    states, count = GadgetSearch.get_state_space(RydbergUnweightedModel, g)
    rydberg_states, rydberg_count = GadgetSearch.get_state_space(GadgetSearch.RydbergModel, g)
    
    @test sort(states) == sort(rydberg_states)
    @test count == rydberg_count
end

@testset "RydbergUnweightedModel _find_weights via α-tensor" begin
    @testset "feasible case" begin
        # Path: 1--2--3, pins [1, 3]
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        # MIS states: {1,3}=0b101, {2}=0b010
        # Ground config for pins [1,3] should be {11} (both pins in, state {1,3})
        # Target: {1,3}=0b101, Wrong: {2}=0b010
        target_states = UInt32[0b101]
        wrong_states = UInt32[0b010]
        edge_list = [(1,2), (2,3)]
        
        result = GadgetSearch._find_weights(
            RydbergUnweightedModel,
            3, edge_list, [1, 3], target_states, wrong_states,
            nothing, nothing, nothing, false, g, false
        )
        
        @test result !== nothing
        @test result == ones(Float64, 3)
    end
    
    @testset "infeasible case - wrong ground config" begin
        # Path: 1--2--3, pins [1, 3]
        # Ground config is {11}, but we claim {2}=0b010 is target (pin config = {00})
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        target_states = UInt32[0b010]  # {2}, pin config 00
        wrong_states = UInt32[0b101]   # {1,3}, pin config 11
        edge_list = [(1,2), (2,3)]
        
        result = GadgetSearch._find_weights(
            RydbergUnweightedModel,
            3, edge_list, [1, 3], target_states, wrong_states,
            nothing, nothing, nothing, false, g, false
        )
        
        @test result === nothing
    end
end

@testset "RydbergUnweightedModel Search Integration" begin
    @testset "search_gadgets without optimizer" begin
        # Create a connected graph: path 1--2--3
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
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
    
    @testset "other models still require optimizer" begin
        # Must use a connected graph so the filter doesn't short-circuit
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        g6codes = [GadgetSearch.GraphIO.Graph6._graphToG6String(g)[11:end]]
        layouts = [[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]]
        
        dataset = GadgetSearch.GraphDataset(g6codes, layouts)
        loader = GadgetSearch.GraphLoader(dataset; pinset=[1, 2, 3])
        
        constraints = [GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))]
        
        # RydbergModel without optimizer should error
        @test_throws ErrorException search_gadgets(
            RydbergModel,
            loader,
            constraints;
            pin_candidates=[[2, 3]],
            max_result_num=1
        )
    end
end
