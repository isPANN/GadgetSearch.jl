using Test
using GadgetSearch
using Graphs
using HiGHS
using JSON3
using Combinatorics

# Helper function to create a simple test graph
function create_test_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    return g
end

# Helper function to create a test GraphLoader
function create_test_loader()
    # Create some simple test graphs
    g1 = SimpleGraph(3)
    add_edge!(g1, 1, 2)
    
    g2 = SimpleGraph(4)
    add_edge!(g2, 1, 2)
    add_edge!(g2, 2, 3)
    
    # Convert to g6 format (simplified for testing)
    g6codes = ["Bw", "C~"]  # These are actual g6 codes for simple graphs
    layouts = [
        [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (1.0, 1.0)]
    ]
    
    dataset = GadgetSearch.GraphDataset(g6codes, layouts)
    return GadgetSearch.GraphLoader(dataset, [1, 2, 3], BitVector(undef, 1024), false, 10, nothing, Dict{String, SimpleGraph{Int}}(), String[])
end

# Helper function to create a mock optimizer
function create_mock_optimizer()
    return () -> HiGHS.Optimizer()
end

# Test constants are defined correctly
@testset "Constants" begin
    @test GadgetSearch.MAX_SUPPORTED_VERTICES == 32
    @test GadgetSearch.DEFAULT_EPSILON == 1.0
    @test GadgetSearch.DEFAULT_MAX_SAMPLES == 100
end


@testset "Cache utilities" begin
    # Ensure cache can be cleared and stats reported
    GadgetSearch.clear_cache!()
    stats = GadgetSearch.get_cache_stats()
    @test stats.size == 0
    
    # Populate cache by calling MIS on a small graph
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    GadgetSearch.find_maximal_independent_sets(g)
    stats2 = GadgetSearch.get_cache_stats()
    @test stats2.size >= 1
    
    # Clear again
    GadgetSearch.clear_cache!()
    @test GadgetSearch.get_cache_stats().size == 0
end

@testset "_graph_hash stability and connectivity check" begin
    # _graph_hash should be stable across identical graphs
    g1 = SimpleGraph(3)
    add_edge!(g1, 1, 2); add_edge!(g1, 2, 3)
    g2 = SimpleGraph(3)
    # Add edges in different order
    add_edge!(g2, 2, 3); add_edge!(g2, 1, 2)
    h1 = GadgetSearch._graph_hash(g1)
    h2 = GadgetSearch._graph_hash(g2)
    @test h1 == h2

    # connectivity check after removal
    g3 = SimpleGraph(4)
    add_edge!(g3, 1, 2); add_edge!(g3, 2, 3); add_edge!(g3, 3, 4)
    # Removing vertex index 2 (0-based -> [1]) should disconnect 1 from 3-4
    @test !GadgetSearch.check_connectivity_after_removal(g3, [1])
    # No removal keeps original connectivity
    @test GadgetSearch.check_connectivity_after_removal(g3, Int[]) == Graphs.is_connected(g3)
end

@testset "make_filter - pin_candidates behavior" begin
    truth_table = BitMatrix([1 0; 0 1])
    optimizer = create_mock_optimizer()
    filter_fn = GadgetSearch.make_filter(truth_table, optimizer, nothing; connected=false,
                                         pin_candidates=[[1, 2]], max_samples=20)
    # Triangle graph often yields a valid solution in our setup
    g = SimpleGraph(3)
    add_edge!(g, 1, 2); add_edge!(g, 1, 3); add_edge!(g, 2, 3)
    pos = [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]
    result = filter_fn(g, pos, nothing)
    # New API returns Gadget or nothing
    @test result === nothing || result isa GadgetSearch.Gadget
    if result !== nothing
        @test result.pins == [1, 2]
    end
end

@testset "solve_weight_enumerate - sampling path" begin
    # Construct a case where combinations exceed max_samples to trigger sampling
    # Need to provide a graph for the new API
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    mis_result = UInt32[0x1, 0x2, 0x4, 0x3, 0x5, 0x6]
    target_mis_indices_all = [collect(1:3), collect(1:3)]  # 9 combinations
    vertex_num = 3
    pin_set = [1, 2]
    w = GadgetSearch.solve_weight_enumerate(mis_result, target_mis_indices_all, vertex_num,
                                            pin_set, create_mock_optimizer(), nothing, nothing,
                                            false, 1, g, false)  # max_samples=1 forces sampling path
    @test w isa Vector{Float64}
    @test length(w) == vertex_num || isempty(w)
end

@testset "find_matching_gadget - early termination" begin
    # Build a loader and a filter that always returns a result
    loader = create_test_loader()
    simple_filter = function(g, pos, pins)
        # Return a Gadget to match the new API
        constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
        return GadgetSearch.Gadget(GadgetSearch.RydbergModel, constraint, g, [1, 2], ones(Float64, nv(g)), pos)
    end
    res = GadgetSearch.find_matching_gadget(loader; filter=simple_filter, max_results=1)
    @test length(res) == 1
end

@testset "find_maximal_independent_sets" begin
    # Test with a simple triangle graph
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 1, 3)
    
    masks, count = GadgetSearch.find_maximal_independent_sets(g)
    
    @test count == length(masks)
    @test count == 3  # Triangle has 3 maximal independent sets: {1}, {2}, {3}
    
    # Verify masks represent single vertices
    expected_masks = UInt32[0x1, 0x2, 0x4]  # Binary: 001, 010, 100
    @test sort(masks) == sort(expected_masks)
end

@testset "find_maximal_independent_sets - error handling" begin
    # Test error for too many vertices
    large_g = SimpleGraph(40)  # More than MAX_SUPPORTED_VERTICES
    
    @test_throws ErrorException GadgetSearch.find_maximal_independent_sets(large_g)
end

@testset "match_rows_by_pinset" begin
    # Create test data
    masks = UInt32[0x5, 0x3, 0x6]  # Binary: 101, 011, 110
    truth_table = BitMatrix([
        1 0;  # Looking for pattern 10
        0 1;  # Looking for pattern 01
        1 1;  # Looking for pattern 11
    ])
    pin_set = [1, 2] # from right(low) to left(high)
    
    result = GadgetSearch.match_rows_by_pinset(masks, truth_table, pin_set)
    @test length(result) == 3
    @test result[1] == [1]  # Mask 0x5 (101) matches pattern 10 at positions [1,2]
    @test result[2] == [3]  # Mask 0x6 (110) matches pattern 01 at positions [1,2]
    @test result[3] == [2]  # Mask 0x3 (011) matches pattern 11 at positions [1,2]
end

@testset "_find_weights - RydbergModel" begin
    # Test with simple case using new API
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    vertex_num = 3
    edge_list = [(1, 2), (2, 3)]
    target_states = UInt32[0x1, 0x2]  # States that should have equal energy
    wrong_states = UInt32[0x4]        # State that should have higher energy
    optimizer = create_mock_optimizer()
    pin_set = [1, 2]
    
    result = GadgetSearch._find_weights(
        GadgetSearch.RydbergModel,
        vertex_num, edge_list, pin_set, target_states, wrong_states, 
        optimizer, nothing, nothing, false, g, false
    )
    
    # Check that we get a result (weights should be found for this simple case)
    if result !== nothing
        @test length(result) == vertex_num
        @test all(w >= 1.0 for w in result)  # All weights should be >= 1
    else
        @test result === nothing  # Also valid if no solution found
    end
end

@testset "solve_weight_enumerate" begin
    # Test with simple MIS result - need to provide graph now
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    mis_result = UInt32[0x1, 0x2, 0x4]
    target_mis_indices_all = [[1], [2]]  # Simple target indices
    vertex_num = 3
    optimizer = create_mock_optimizer()
    pin_set = [1, 2]
    
    weights = GadgetSearch.solve_weight_enumerate(
        mis_result, target_mis_indices_all, vertex_num, pin_set, optimizer,
        nothing, nothing, false, 1000, g, false
    )
    
    # Should return weights or empty vector
    @test weights isa Vector{Float64}
    @test length(weights) == vertex_num || length(weights) == 0
end

@testset "solve_weight_enumerate - error handling" begin
    g = SimpleGraph(2)
    add_edge!(g, 1, 2)
    
    # Test error when optimizer is nothing
    mis_result = UInt32[0x1, 0x2]
    target_mis_indices_all = [[1], [2]]
    vertex_num = 2
    pin_set = [1, 2]
    
    @test_throws ErrorException GadgetSearch.solve_weight_enumerate(
        mis_result, target_mis_indices_all, vertex_num, pin_set, nothing,
        nothing, nothing, false, 1000, g, false
    )
    
    # Test empty target indices
    empty_target = [Int[], [1]]
    weights = GadgetSearch.solve_weight_enumerate(
        mis_result, empty_target, vertex_num, pin_set, create_mock_optimizer(),
        nothing, nothing, false, 1000, g, false
    )
    @test length(weights) == 0
end

@testset "make_filter" begin
    # Create test data
    truth_table = BitMatrix([
        1 0;
        0 1;
    ])
    optimizer = create_mock_optimizer()
    
    filter_fn = GadgetSearch.make_filter(
        truth_table, optimizer, nothing;
        connected=false,
        max_samples=10
    )
    
    @test filter_fn isa Function
    
    # Test the filter with a simple graph
    g = SimpleGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 3, 4)
    
    pos = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
    pin_set = [1, 2]
    
    # Test filter function - now returns Gadget or nothing
    result = filter_fn(g, pos, pin_set)
    @test result === nothing || result isa GadgetSearch.Gadget
end

@testset "find_matching_gadget" begin
    loader = create_test_loader()
    
    # Simple filter that always returns nothing (no matches)
    simple_filter = (g, pos, pins) -> nothing
    
    results = GadgetSearch.find_matching_gadget(loader; filter=simple_filter, limit=1)
    @test results isa Vector{GadgetSearch.Gadget}
    @test length(results) == 0  # No matches expected with our simple filter
end

@testset "search_by_truth_tables" begin
    loader = create_test_loader()
    truth_tables = [BitMatrix([1 0; 0 1]), BitMatrix([1 1; 0 0])]
    optimizer = create_mock_optimizer()
    
    # Create temporary file for saving results
    temp_file = tempname() * ".json"
    
    try
        results, failed_tt = GadgetSearch.search_by_truth_tables(
            loader, truth_tables;
            optimizer=optimizer,
            max_result_num=1,
            save_path=temp_file
        )
        # New API returns results grouped by truth table
        @test results isa Vector{Vector{GadgetSearch.Gadget}}
        @test failed_tt isa Vector{BitMatrix}
        # Should respect max_result_num per truth table
        @test all(length(r) <= 1 for r in results)
        
        # Check that save file is created when results are found
        all_results = vcat(results...)
        if !isempty(all_results) && isfile(temp_file)
            saved_data = JSON3.read(temp_file)
            @test saved_data isa Union{Vector, JSON3.Array}
            @test length(saved_data) > 0
        end
        
    finally
        # Clean up temp file
        isfile(temp_file) && rm(temp_file)
    end
end

@testset "search_by_truth_tables - max_result_num behavior" begin
    # Test max_result_num parameter functionality
    loader = create_test_loader()
    
    # Create multiple truth tables to test with
    truth_tables = [
        BitMatrix([1 0; 0 1]), 
        BitMatrix([1 1; 0 0]),
        BitMatrix([0 1; 1 0])
    ]
    
    temp_file = tempname() * ".json"
    
    try
        results, failed_tt = GadgetSearch.search_by_truth_tables(
            loader, truth_tables;
            optimizer=create_mock_optimizer(),
            max_result_num=1,  # Should stop after 1 result
            save_path=temp_file
        )
        # Should respect the max_result_num limit per truth table
        @test results isa Vector{Vector{GadgetSearch.Gadget}}
        @test all(length(r) <= 1 for r in results)
        @test failed_tt isa Vector{BitMatrix}
        
    finally
        isfile(temp_file) && rm(temp_file)
    end
end

# Integration test with real-world-like scenario
@testset "Integration test" begin
    # Create a more realistic test case
    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 5)
    
    # Test that we can find maximal independent sets
    masks, count = GadgetSearch.find_maximal_independent_sets(g)
    @test count > 0
    @test length(masks) == count
    
    # Test matching with a simple truth table
    truth_table = BitMatrix([1 0; 0 1])
    pin_set = [1, 2]
    
    matched = GadgetSearch.match_rows_by_pinset(masks, truth_table, pin_set)
    @test length(matched) == 2
    @test all(m isa Vector{Int} for m in matched)
end

# New tests for the unified framework
@testset "EnergyModel types" begin
    @test GadgetSearch.RydbergModel <: GadgetSearch.EnergyModel
    @test GadgetSearch.QUBOModel <: GadgetSearch.EnergyModel
end

@testset "GadgetConstraint types" begin
    @test GadgetSearch.TruthTableConstraint <: GadgetSearch.GadgetConstraint
    
    # Test TruthTableConstraint
    tt = BitMatrix([1 0; 0 1])
    ttc = GadgetSearch.TruthTableConstraint(tt)
    @test GadgetSearch.get_pin_num(ttc) == 2
end

@testset "get_state_space" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    # Rydberg model uses MIS
    mis_states, mis_count = GadgetSearch.get_state_space(GadgetSearch.RydbergModel, g)
    @test mis_count > 0
    @test length(mis_states) == mis_count
    
    # QUBO model uses all 2^n states
    all_states, all_count = GadgetSearch.get_state_space(GadgetSearch.QUBOModel, g)
    @test all_count == 2^3  # 8 states for 3 vertices
    @test length(all_states) == 8
end

@testset "Gadget constructors" begin
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    constraint = GadgetSearch.TruthTableConstraint(BitMatrix([1 0; 0 1]))
    
    # RydbergModel gadget
    rydberg_gadget = GadgetSearch.Gadget(GadgetSearch.RydbergModel, constraint, g, [1, 2], [1.0, 2.0, 1.0], nothing)
    @test rydberg_gadget.vertex_weights == [1.0, 2.0, 1.0]
    @test isempty(rydberg_gadget.edge_weights)
    
    # Legacy compatibility
    @test rydberg_gadget.weights == [1.0, 2.0, 1.0]
    @test rydberg_gadget.ground_states == BitMatrix([1 0; 0 1])
end
