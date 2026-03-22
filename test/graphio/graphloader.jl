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


# Test save_cache with no cache path
@testset "GraphLoader Cache - no cachepath warning" begin
    g6_codes = ["A_", "Bw"]
    ds = GraphDataset(g6_codes)
    loader = GraphLoader(ds; enable_cache=false)
    # save_cache with no cachepath should warn, not error
    @test_logs (:warn, "No cache path provided; cannot save cache.") save_cache(loader)
end

# Test cache file corruption handling
@testset "GraphLoader Cache - corrupted cache file" begin
    test_file = create_test_graph_file()
    bad_cache = tempname()
    write(bad_cache, "this is not a valid serialized file")

    loader = @test_logs (:info,) (:warn,) match_mode=:any GraphLoader(test_file; enable_cache=true, cachepath=bad_cache)
    @test length(loader.parsed) == 0  # corrupted cache ignored

    rm(test_file)
    rm(bad_cache)
end

# Test cache file with wrong type
@testset "GraphLoader Cache - wrong type in cache" begin
    test_file = create_test_graph_file()
    bad_cache = tempname()
    using Serialization
    Serialization.serialize(bad_cache, [1, 2, 3])  # wrong type

    loader = @test_logs (:info,) (:warn, "Cache file format unexpected, ignoring cache.") GraphLoader(test_file; enable_cache=true, cachepath=bad_cache)
    @test length(loader.parsed) == 0

    rm(test_file)
    rm(bad_cache)
end

# Test GraphLoader string key error
@testset "GraphLoader - invalid string key" begin
    g6_codes = ["A_"]
    ds = GraphDataset(g6_codes)
    loader = GraphLoader(ds)
    @test_throws ErrorException loader["abc"]  # not a valid integer
end

# Test GraphDataset mismatched lengths
@testset "GraphDataset - mismatched g6codes and layouts" begin
    @test_throws ArgumentError GraphDataset(["A_", "Bw"], [[(0.0, 0.0), (1.0, 0.0)]])
end

# Test GraphLoader show method
@testset "GraphLoader show" begin
    g6_codes = ["A_", "Bw"]
    ds = GraphDataset(g6_codes)
    loader = GraphLoader(ds)
    str = sprint(show, loader)
    @test occursin("2 graphs", str)
    @test occursin("no cache", str)

    loader2 = GraphLoader(ds; enable_cache=true)
    str2 = sprint(show, loader2)
    @test occursin("cache enabled", str2)
end

# Test GraphLoader cache eviction
@testset "GraphLoader Cache - eviction" begin
    g6_codes = ["A_", "Bw", "C~"]
    ds = GraphDataset(g6_codes)
    loader = GraphLoader(ds; enable_cache=true, max_cached=2)

    _ = loader[1]
    _ = loader[2]
    @test length(loader.parsed) == 2

    _ = loader[3]
    @test length(loader.parsed) == 2  # oldest evicted

    # Access 1 again (should not be in cache, re-parsed)
    g = loader[1]
    @test nv(g) == 2
end

# Test file with invalid coordinates
@testset "GraphDataset - invalid coordinates in file" begin
    test_data = "A_ not_a_coordinate"
    test_file = tempname()
    write(test_file, test_data)

    ds = GraphDataset(test_file)
    @test ds.n == 1
    @test ds.g6codes[1] == "A_"
    @test ds.layouts[1] === nothing  # invalid coords → nothing

    rm(test_file)
end

# Test file with g6-only lines (no coordinates)
@testset "GraphDataset - g6 only lines" begin
    test_data = "A_\nBw"
    test_file = tempname()
    write(test_file, test_data)

    ds = GraphDataset(test_file)
    @test ds.n == 2
    @test ds.layouts[1] === nothing
    @test ds.layouts[2] === nothing

    rm(test_file)
end

# Test file with blank/whitespace lines
@testset "GraphDataset - blank lines in file" begin
    test_data = "A_ (0.0,0.0);(1.0,0.0)\n\n   \nBw (0.0,0.0);(1.0,0.0);(0.5,1.0)"
    test_file = tempname()
    write(test_file, test_data)

    ds = GraphDataset(test_file)
    @test ds.n == 2  # blank lines skipped

    rm(test_file)
end

# Test GraphLoader with pinset
@testset "GraphLoader - pinset" begin
    g6_codes = ["A_"]
    ds = GraphDataset(g6_codes)
    loader = GraphLoader(ds; pinset=[1, 2])
    @test loader.pinset == [1, 2]
end

# Test Graph6 parsing functions
@testset "Graph6 Parsing" begin
    @testset "Vertex Count Parsing" begin
        # Small graphs: single character encoding, n = char - 63
        n, pos = GadgetSearch._parse_vertex_count("?")   # '?' = 63, n = 0
        @test n == 0 && pos == 2

        n, pos = GadgetSearch._parse_vertex_count("@")   # '@' = 64, n = 1
        @test n == 1 && pos == 2

        n, pos = GadgetSearch._parse_vertex_count("A")   # 'A' = 65, n = 2
        @test n == 2 && pos == 2

        n, pos = GadgetSearch._parse_vertex_count("B")   # 'B' = 66, n = 3
        @test n == 3 && pos == 2

        n, pos = GadgetSearch._parse_vertex_count("C")   # 'C' = 67, n = 4
        @test n == 4 && pos == 2
    end

    @testset "G6 String Parsing" begin
        temp_bitvec = BitVector(undef, 100)

        # n=0: "?" returns empty graph
        g = GadgetSearch._parse_g6_string("?", temp_bitvec)
        @test nv(g) == 0
        @test ne(g) == 0

        # n=1: "@" returns single vertex, no edges
        g = GadgetSearch._parse_g6_string("@", temp_bitvec)
        @test nv(g) == 1
        @test ne(g) == 0

        # n=2: "A?" encodes no edge (val=0)
        g = GadgetSearch._parse_g6_string("A?", temp_bitvec)
        @test nv(g) == 2
        @test ne(g) == 0
        @test !has_edge(g, 1, 2)

        # n=2: "A_" encodes edge 1-2 ('_'=95, val=32=0b100000, bit5=1)
        g = GadgetSearch._parse_g6_string("A_", temp_bitvec)
        @test nv(g) == 2
        @test ne(g) == 1
        @test has_edge(g, 1, 2)

        # n=3: "Bw" encodes complete triangle ('w'=119, val=56=0b111000)
        g = GadgetSearch._parse_g6_string("Bw", temp_bitvec)
        @test nv(g) == 3
        @test ne(g) == 3
        @test has_edge(g, 1, 2)
        @test has_edge(g, 1, 3)
        @test has_edge(g, 2, 3)

        # Header stripping: ">>graph6<<A_" should behave like "A_"
        g = GadgetSearch._parse_g6_string(">>graph6<<A_", temp_bitvec)
        @test nv(g) == 2
        @test ne(g) == 1
        @test has_edge(g, 1, 2)

        # Empty string throws ArgumentError
        @test_throws ArgumentError GadgetSearch._parse_g6_string("", temp_bitvec)
    end

    @testset "G6 Integration with GraphLoader" begin
        # "A?" = 2 vertices, no edge; "A_" = 2 vertices, edge 1-2; "Bw" = triangle
        g6_codes = ["A?", "A_", "Bw"]
        dataset = GraphDataset(g6_codes)
        loader = GraphLoader(dataset)

        g1 = loader[1]   # A?: 2 vertices, no edges
        @test nv(g1) == 2
        @test ne(g1) == 0

        g2 = loader[2]   # A_: 2 vertices, 1 edge
        @test nv(g2) == 2
        @test ne(g2) == 1

        g3 = loader[3]   # Bw: complete triangle
        @test nv(g3) == 3
        @test ne(g3) == 3
    end
end

@testset "Graph6 Encoding" begin
    @testset "roundtrip without header" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 4)

        g6 = graph_to_g6(g)
        @test !startswith(g6, ">>graph6<<")

        parsed = GadgetSearch._parse_g6_string(g6, BitVector(undef, 0))
        @test nv(parsed) == nv(g)
        @test ne(parsed) == ne(g)
        @test collect(edges(parsed)) == collect(edges(g))
    end

    @testset "header option" begin
        g = SimpleGraph(2)
        add_edge!(g, 1, 2)

        with_header = graph_to_g6(g; include_header=true)
        @test startswith(with_header, ">>graph6<<")

        parsed = GadgetSearch._parse_g6_string(with_header, BitVector(undef, 0))
        @test nv(parsed) == 2
        @test ne(parsed) == 1
        @test has_edge(parsed, 1, 2)
    end

    @testset "extended vertex-count encoding" begin
        g = SimpleGraph(70) # requires extended graph6 vertex-count prefix
        add_edge!(g, 1, 70)
        add_edge!(g, 2, 69)

        g6 = graph_to_g6(g)
        parsed = GadgetSearch._parse_g6_string(g6, BitVector(undef, 0))
        @test nv(parsed) == 70
        @test ne(parsed) == 2
        @test has_edge(parsed, 1, 70)
        @test has_edge(parsed, 2, 69)
    end
end
