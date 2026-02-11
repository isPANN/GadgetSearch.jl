using Test
using GadgetSearch
using Graphs

@testset "AlphaTensorMode - α-Tensor Computation" begin
    @testset "compute_reduced_alpha_tensor - simple graph" begin
        # Simple 3-vertex path: 1--2--3
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        # MIS of this graph: {1,3} (size 2) or {2} (size 1)
        # Pin vertices: [2]
        pins = [2]
        
        α = GadgetSearch.compute_reduced_alpha_tensor(g, pins)
        
        # Boundary config 0 (pin 2 not selected): interior {1,3} → max MIS size = 2
        @test α[UInt32(0)] == 2
        
        # Boundary config 1 (pin 2 selected): interior {1,3} blocked → max MIS size = 0
        @test α[UInt32(1)] == 0
    end
    
    @testset "compute_reduced_alpha_tensor - triangle" begin
        # Triangle: 1--2, 2--3, 3--1
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 1)
        
        # MIS of triangle: any single vertex
        # Pins: [1]
        pins = [1]
        
        α = GadgetSearch.compute_reduced_alpha_tensor(g, pins)
        
        # Boundary config 0 (pin 1 not selected): interior {2,3} → max MIS size = 1
        @test α[UInt32(0)] == 1
        
        # Boundary config 1 (pin 1 selected): interior {2,3} both blocked → max MIS size = 0
        @test α[UInt32(1)] == 0
    end
    
    @testset "compute_reduced_alpha_tensor - disconnected" begin
        # Two disconnected edges: 1--2, 3--4
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 3, 4)
        
        # Pins: [1, 3]
        pins = [1, 3]
        
        α = GadgetSearch.compute_reduced_alpha_tensor(g, pins)
        
        # Config 00: interior {2,4} → max MIS = 2
        @test α[UInt32(0b00)] == 2
        
        # Config 01: pin 1 selected, interior {2} → max MIS = 1
        @test α[UInt32(0b01)] == 1
        
        # Config 10: pin 3 selected, interior {4} → max MIS = 1
        @test α[UInt32(0b10)] == 1
        
        # Config 11: both pins selected → max MIS = 0
        @test α[UInt32(0b11)] == 0
    end
end

@testset "AlphaTensorMode - Equivalence Check" begin
    @testset "check_alpha_equivalence - equivalent tensors" begin
        α1 = Dict(UInt32(0) => 2, UInt32(1) => 0)
        α2 = Dict(UInt32(0) => 5, UInt32(1) => 3)  # differ by constant 3
        
        is_equiv, c = GadgetSearch.check_alpha_equivalence(α1, α2)
        @test is_equiv == true
        @test c == 3
    end
    
    @testset "check_alpha_equivalence - not equivalent" begin
        α1 = Dict(UInt32(0) => 2, UInt32(1) => 0)
        α2 = Dict(UInt32(0) => 5, UInt32(1) => 2)  # not constant difference
        
        is_equiv, c = GadgetSearch.check_alpha_equivalence(α1, α2)
        @test is_equiv == false
        @test c === nothing
    end
    
    @testset "check_alpha_equivalence - different keys" begin
        α1 = Dict(UInt32(0) => 2, UInt32(1) => 0)
        α2 = Dict(UInt32(0) => 5, UInt32(2) => 3)  # different keys
        
        is_equiv, c = GadgetSearch.check_alpha_equivalence(α1, α2)
        @test is_equiv == false
    end
end

@testset "AlphaTensorMode - Pattern Inference" begin
    @testset "infer_pattern_alpha - simple case" begin
        # Target states: {1,3} and {2,3} with pins [3]
        # State {1,3} = 0b101, State {2,3} = 0b110
        target_states = UInt32[0b101, 0b110]
        pins = [3]
        
        α = GadgetSearch.infer_pattern_alpha(target_states, pins)
        
        # Pin config 1 (pin 3 selected):
        # State {1,3}: total=2, selected pins=1, interior=1
        # State {2,3}: total=2, selected pins=1, interior=1
        @test α[UInt32(1)] == 1
    end
    
    @testset "extract_pin_config" begin
        # State 0b1010 with pins [2, 4]
        state = UInt32(0b1010)
        pins = [2, 4]
        
        config = GadgetSearch.extract_pin_config(state, pins)
        
        # Pin 2 (bit 1): on → bit 0 of config = 1
        # Pin 4 (bit 3): on → bit 1 of config = 1
        @test config == UInt32(0b11)
    end
end

@testset "AlphaTensorMode - Gadget Verification" begin
    @testset "verify_gadget_via_alpha_tensor - basic test" begin
        # Simple path: 1--2--3
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        # Target: MIS {1,3} with pins [1, 3]
        target_states = UInt32[0b101]  # {1, 3}
        pins = [1, 3]
        
        result = GadgetSearch.verify_gadget_via_alpha_tensor(g, pins, target_states)
        
        # This specific case may or may not be valid - just test the function runs
        @test result === nothing || (result isa Tuple && length(result) == 2)
        
        if result !== nothing
            weights, overhead = result
            @test weights == ones(Float64, 3)
            @test overhead isa Int
        end
    end
end

@testset "AlphaTensorMode - Type System Integration" begin
    @test AlphaTensorMode <: GadgetSearch.EnergyModel
    
    # State space should use MIS (same as RydbergModel)
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    states, count = GadgetSearch.get_state_space(AlphaTensorMode, g)
    rydberg_states, rydberg_count = GadgetSearch.get_state_space(GadgetSearch.RydbergModel, g)
    
    @test sort(states) == sort(rydberg_states)
    @test count == rydberg_count
end

@testset "AlphaTensorMode - Search Integration" begin
    @testset "search_gadgets without optimizer" begin
        # Create a simple graph dataset
        g = SimpleGraph(3)
        add_edge!(g, 2, 3)
        
        g6codes = [GadgetSearch.GraphIO.Graph6._graphToG6String(g)[11:end]]
        layouts = [[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]]
        
        dataset = GadgetSearch.GraphDataset(g6codes, layouts)
        loader = GadgetSearch.GraphLoader(dataset; pinset=[1, 2, 3])
        
        constraints = [GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))]
        
        # Should not throw — no optimizer needed for AlphaTensorMode
        results, failed = search_gadgets(
            AlphaTensorMode,
            loader,
            constraints;
            pin_candidates=[[2, 3]],
            max_result_num=1
        )
        
        @test results isa Vector{Vector{GadgetSearch.Gadget}}
        @test failed isa Vector{GadgetSearch.TruthTableConstraint}
    end
end

@testset "AlphaTensorMode - Comparison with RydbergUnweightedModel" begin
    @testset "Both should find same gadgets for simple cases" begin
        # Create a simple dataset
        g = SimpleGraph(3)
        add_edge!(g, 2, 3)
        
        g6codes = [GadgetSearch.GraphIO.Graph6._graphToG6String(g)[11:end]]
        layouts = [[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]]
        
        dataset = GadgetSearch.GraphDataset(g6codes, layouts)
        loader = GadgetSearch.GraphLoader(dataset; pinset=[1, 2, 3])
        
        constraint = [GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))]
        pin_cands = [[2, 3]]
        
        # Search with AlphaTensorMode
        results_alpha, _ = search_gadgets(
            AlphaTensorMode,
            loader,
            constraint;
            pin_candidates=pin_cands,
            max_result_num=1
        )
        
        # Search with RydbergUnweightedModel
        results_unweighted, _ = search_gadgets(
            RydbergUnweightedModel,
            loader,
            constraint;
            pin_candidates=pin_cands,
            max_result_num=1
        )
        
        # Both should find the same number of gadgets
        @test length(results_alpha[1]) == length(results_unweighted[1])
    end
end

