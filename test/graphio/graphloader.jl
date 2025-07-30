using Test
using GadgetSearch
using Graphs

# Create a temporary test file with some graph data
function create_test_graph_file()
    test_data = """A_ (0.0, 0.0); (1.0, 0.0)
B_ (0.0, 0.0); (1.0, 0.0); (0.5, 1.0)
C_ (0.0, 0.0); (1.0, 0.0); (0.5, 1.0); (0.5, 0.5)"""
    test_file = tempname()
    write(test_file, test_data)
    return test_file
end

# Test GraphDataset creation
@testset "GraphDataset" begin
    test_file = create_test_graph_file()
    ds = GraphDataset(test_file)
    
    @test ds.n == 3
    @test length(ds.g6codes) == 3
    @test ds.g6codes[1] == "A_"
    @test ds.g6codes[2] == "B_"
    @test ds.g6codes[3] == "C_"
    
    # Test layouts
    @test ds.layouts[1] == [(0.0, 0.0), (1.0, 0.0)]
    @test ds.layouts[2] == [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]
    @test ds.layouts[3] == [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.5, 0.5)]
    
    rm(test_file)
end

# Test GraphLoader basic functionality
@testset "GraphLoader Basic" begin
    test_file = create_test_graph_file()
    loader = GraphLoader(test_file)
    
    @test length(loader) == 3
    @test collect(keys(loader)) == 1:3
    
    # Test graph loading
    g1 = loader[1]
    @test g1 isa SimpleGraph{Int}
    
    # Test string indexing
    g1_str = loader["1"]
    @test g1_str isa SimpleGraph{Int}
    
    rm(test_file)
end

# Test caching functionality
@testset "GraphLoader Cache" begin
    test_file = create_test_graph_file()
    cache_file = tempname()
    
    # Test with cache enabled
    loader = GraphLoader(test_file; enable_cache=true, cachepath=cache_file, max_cached=2)
    
    # Load graphs to fill cache
    g1 = loader[1]
    g2 = loader[2]
    g3 = loader[3]
    
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
    loader = GraphLoader(test_file)
    
    # Test layout access
    @test loader.layout[1] == [(0.0, 0.0), (1.0, 0.0)]
    @test loader.layout[2] == [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0)]
    @test loader.layout[3] == [(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.5, 0.5)]
    
    # Test string key access
    @test loader.layout["1"] == [(0.0, 0.0), (1.0, 0.0)]
    
    rm(test_file)
end

# Test flexible constructors
@testset "GraphDataset Flexible Constructors" begin
    # Test constructor with just g6 codes
    g6_codes = ["A_", "B_"]
    ds1 = GraphDataset(g6_codes)
    
    @test ds1.n == 2
    @test length(ds1.g6codes) == 2
    @test ds1.g6codes[1] == "A_"
    @test all(layout === nothing for layout in ds1.layouts)
    
    # Test constructor with g6 codes and layouts
    layouts = [[(0.0, 0.0), (1.0, 0.0)], [(0.0, 0.0), (1.0, 1.0)]]
    ds2 = GraphDataset(g6_codes, layouts)
    
    @test ds2.n == 2
    @test ds2.layouts[1] == [(0.0, 0.0), (1.0, 0.0)]
    @test ds2.layouts[2] == [(0.0, 0.0), (1.0, 1.0)]
    
    # Test GraphLoader with GraphDataset
    loader = GraphLoader(ds2)
    @test length(loader) == 2
    @test loader.layout[1] == [(0.0, 0.0), (1.0, 0.0)]
end


# Test Graph6 parsing functions
@testset "Graph6 Parsing" begin
    # TODO
    # Test _parse_vertex_count function
    # @testset "Vertex Count Parsing" begin
    #     # Small graphs (n <= 62): single character encoding
    #     nv, pos = GadgetSearch._parse_vertex_count("A")  # 'A' = 65, so n = 65-63 = 2
    #     @test nv == 2
    #     @test pos == 2
        
    #     nv, pos = GadgetSearch._parse_vertex_count("?")  # '?' = 63, so n = 63-63 = 0
    #     @test nv == 0
    #     @test pos == 2
        
    #     nv, pos = GadgetSearch._parse_vertex_count("~")  # '~' = 126, so n = 126-63 = 63
    #     @test nv == 63
    #     @test pos == 2
        
    #     # Medium graphs (63 <= n <= 258047): 4-character encoding starting with "~"
    #     # For n=100: encode as ~AAd where 100 = 1*4096 + 1*64 + 36*1
    #     nv, pos = GadgetSearch._parse_vertex_count("~AAd")
    #     @test nv == 100
    #     @test pos == 5
        
    #     # Large graphs (n >= 258048): 7-character encoding starting with "~~"
    #     # This is rarely used but let's test the boundary
    #     nv, pos = GadgetSearch._parse_vertex_count("~~?????")  # Minimum large graph encoding
    #     @test nv == 258048
    #     @test pos == 9
    # end
    
    # @testset "G6 String Parsing" begin
    #     temp_bitvec = BitVector(undef, 100)
        
    #     # Test empty graph
    #     g = GadgetSearch._parse_g6_string("?", temp_bitvec)
    #     @test nv(g) == 0
    #     @test ne(g) == 0
        
    #     # Test single vertex graph
    #     g = GadgetSearch._parse_g6_string("@", temp_bitvec)  # '@' = 64, so n = 64-63 = 1
    #     @test nv(g) == 1
    #     @test ne(g) == 0
        
    #     # Test two vertex graphs
    #     g = GadgetSearch._parse_g6_string("A_", temp_bitvec)  # No edge between vertices 1,2
    #     @test nv(g) == 2
    #     @test ne(g) == 0
    #     @test !has_edge(g, 1, 2)
        
    #     g = GadgetSearch._parse_g6_string("A?", temp_bitvec)  # Edge between vertices 1,2
    #     @test nv(g) == 2
    #     @test ne(g) == 1
    #     @test has_edge(g, 1, 2)
        
    #     # Test three vertex graphs
    #     g = GadgetSearch._parse_g6_string("B_", temp_bitvec)  # Triangle with no edges
    #     @test nv(g) == 3
    #     @test ne(g) == 0
        
    #     g = GadgetSearch._parse_g6_string("BW", temp_bitvec)  # Complete triangle
    #     @test nv(g) == 3
    #     @test ne(g) == 3
    #     @test has_edge(g, 1, 2)
    #     @test has_edge(g, 1, 3)
    #     @test has_edge(g, 2, 3)
        
    #     # Test four vertex graphs
    #     g = GadgetSearch._parse_g6_string("C_", temp_bitvec)  # 4 vertices, no edges
    #     @test nv(g) == 4
    #     @test ne(g) == 0
        
    #     # Test graph6 header handling
    #     g = GadgetSearch._parse_g6_string(">>graph6<<A?", temp_bitvec)
    #     @test nv(g) == 2
    #     @test ne(g) == 1
    #     @test has_edge(g, 1, 2)
        
    #     # Test error handling
    #     @test_throws ArgumentError GadgetSearch._parse_g6_string("", temp_bitvec)
    # end
    
    # @testset "G6 Integration with GraphLoader" begin
    #     # Test that GraphLoader correctly uses the g6 parsing
    #     g6_codes = ["A_", "A?", "BW", "C_"]
    #     dataset = GraphDataset(g6_codes)
    #     loader = GraphLoader(dataset)
        
    #     # Test the parsed graphs
    #     g1 = loader[1]  # A_: 2 vertices, no edges
    #     @test nv(g1) == 2
    #     @test ne(g1) == 0
        
    #     g2 = loader[2]  # A?: 2 vertices, 1 edge
    #     @test nv(g2) == 2
    #     @test ne(g2) == 1
        
    #     g3 = loader[3]  # BW: 3 vertices, complete triangle
    #     @test nv(g3) == 3
    #     @test ne(g3) == 3
        
    #     g4 = loader[4]  # C_: 4 vertices, no edges
    #     @test nv(g4) == 4
    #     @test ne(g4) == 0
    # end
end
