using Test
using GadgetSearch
using Graphs
using JSON3

# Create a temporary test file with some graph data
function create_test_graph_file()
    test_data = """
    graph1 DEk
    graph2 DEo
    graph3 DEa
    """
    test_file = tempname()
    write(test_file, test_data)
    return test_file
end

# Create a temporary layout file
function create_test_layout_file()
    test_layout = Dict(
        "graph1" => [1.0, 2.0, 3.0, 4.0],
        "graph2" => [2.0, 3.0, 4.0, 1.0],
        "graph3" => [3.0, 4.0, 1.0, 2.0]
    )
    test_file = tempname()
    write(test_file, JSON3.write(test_layout))
    return test_file
end

@testset "GraphLoader Tests" begin
    # Test GraphDataset creation
    @testset "GraphDataset" begin
        test_file = create_test_graph_file()
        ds = GraphDataset(test_file)
        
        @test length(ds.keys) == 3
        @test ds.keys[1] == "graph1"
        @test ds.keys[2] == "graph2"
        @test ds.keys[3] == "graph3"
        @test ds.key_to_index["graph1"] == 1
        @test ds.key_to_index["graph2"] == 2
        @test ds.key_to_index["graph3"] == 3
        
        rm(test_file)
    end

    # Test GraphLoader basic functionality
    @testset "GraphLoader Basic" begin
        test_file = create_test_graph_file()
        loader = GraphLoader(test_file)
        
        @test length(loader) == 3
        @test collect(keys(loader)) == ["graph1", "graph2", "graph3"]
        
        # Test graph loading
        g1 = loader["graph1"]
        @test g1 isa SimpleGraph{Int}
        
        rm(test_file)
    end

    # Test caching functionality
    @testset "GraphLoader Cache" begin
        test_file = create_test_graph_file()
        cache_file = tempname()
        
        # Test with cache enabled
        loader = GraphLoader(test_file; enable_cache=true, cachepath=cache_file, max_cached=2)
        
        # Load graphs to fill cache
        g1 = loader["graph1"]
        g2 = loader["graph2"]
        g3 = loader["graph3"]
        
        # Test cache size limit
        @test length(loader.parsed) <= 2
        
        # Save and reload cache
        save_cache(loader)
        new_loader = GraphLoader(test_file; enable_cache=true, cachepath=cache_file)
        @test length(new_loader.parsed) > 0
        
        rm(test_file)
        rm(cache_file)
    end

    # Test layout functionality
    @testset "GraphLoader Layout" begin
        test_file = create_test_graph_file()
        layout_file = create_test_layout_file()
        
        loader = GraphLoader(test_file; layoutfile=layout_file)
        
        # Test layout access
        @test loader.layout["graph1"] == [(1.0, 2.0), (3.0, 4.0)]
        @test loader.layout["graph2"] == [(2.0, 3.0), (4.0, 1.0)]
        @test loader.layout["graph3"] == [(3.0, 4.0), (1.0, 2.0)]
        
        # Test layout cache
        @test !isempty(loader.layoutcache)
        
        rm(test_file)
        rm(layout_file)
    end
end 