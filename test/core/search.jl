using Test
using GadgetSearch
using Graphs
using HiGHS
using JSON3
using Combinatorics

# Test constants are defined correctly
@testset "Constants" begin
    @test GadgetSearch.MAX_SUPPORTED_VERTICES == 32
    @test GadgetSearch.DEFAULT_EPSILON == 1.0
    @test GadgetSearch.DEFAULT_MAX_SAMPLES == 100
end

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

@testset "_find_weight" begin
    # Test with simple case
    vertex_num = 3
    target_masks = UInt32[0x1, 0x2]  # Masks that should have equal energy
    wrong_masks = UInt32[0x4]        # Mask that should have higher energy
    optimizer = create_mock_optimizer()
    
    weights = GadgetSearch._find_weight(
        vertex_num, target_masks, wrong_masks, optimizer, nothing, nothing, false
    )
    
    # Check that we get a result (weights should be found for this simple case)
    @test length(weights) == vertex_num || length(weights) == 0
    
    if length(weights) > 0
        @test all(w >= 1.0 for w in weights)  # All weights should be >= 1
    end
end

@testset "solve_weight_enumerate" begin
    # Test with simple MIS result
    mis_result = UInt32[0x1, 0x2, 0x4]
    target_mis_indices_all = [[1], [2]]  # Simple target indices
    vertex_num = 3
    optimizer = create_mock_optimizer()
    
    weights = GadgetSearch.solve_weight_enumerate(
        mis_result, target_mis_indices_all, vertex_num, optimizer
    )
    
    # Should return weights or empty vector
    @test weights isa Vector{Float64}
    @test length(weights) == vertex_num || length(weights) == 0
end

@testset "solve_weight_enumerate - error handling" begin
    # Test error when optimizer is nothing
    mis_result = UInt32[0x1, 0x2]
    target_mis_indices_all = [[1], [2]]
    vertex_num = 2
    
    @test_throws ErrorException GadgetSearch.solve_weight_enumerate(
        mis_result, target_mis_indices_all, vertex_num, nothing
    )
    
    # Test empty target indices
    empty_target = [Int[], [1]]
    weights = GadgetSearch.solve_weight_enumerate(
        mis_result, empty_target, vertex_num, create_mock_optimizer()
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
    
    # Test filter function (may return nothing if no solution found)
    result = filter_fn(g, pos, pin_set)
    @test result isa Tuple{Union{Nothing, Vector{Float64}}, BitMatrix, Union{Nothing, Vector{Int}}}
end

@testset "find_matching_gadget" begin
    loader = create_test_loader()
    
    # Simple filter that always returns nothing (no matches)
    simple_filter = (g, pos, pins) -> (nothing, BitMatrix([1 0; 0 1]), nothing)
    
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
        
        @test results isa Vector{GadgetSearch.Gadget}
        @test failed_tt isa Vector{BitMatrix}
        @test length(results) <= 1  # Should respect max_result_num
        
        # Check that save file is created when results are found
        if !isempty(results) && isfile(temp_file)
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
        
        # Should respect the max_result_num limit
        @test length(results) <= 1
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


