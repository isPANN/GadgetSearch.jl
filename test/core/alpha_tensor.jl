using GadgetSearch
using Graphs
using GenericTensorNetworks: content, Tropical

# Helper to build the CROSS graph: 4 vertices with two crossing edges {1-3, 2-4}
function make_cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    return g
end

@testset "calculate_alpha_tensor" begin
    g = make_cross_graph()
    boundary = [1, 2, 3, 4]
    alpha = calculate_alpha_tensor(g, boundary)
    vals = content.(alpha)

    # Verify against Table 2 from the paper (CROSS alpha tensor)
    # Indexing: (s1+1, s2+1, s3+1, s4+1) for boundary config s1s2s3s4
    @test vals[1,1,1,1] == 0.0   # 0000
    @test vals[2,1,1,1] == 1.0   # 1000
    @test vals[1,2,1,1] == 1.0   # 0100
    @test vals[1,1,2,1] == 1.0   # 0010
    @test vals[1,1,1,2] == 1.0   # 0001
    @test vals[2,2,1,1] == 2.0   # 1100
    @test vals[1,2,2,1] == 2.0   # 0110
    @test vals[1,1,2,2] == 2.0   # 0011
    @test vals[2,1,1,2] == 2.0   # 1001
    @test vals[2,1,2,1] == -Inf  # 1010 (vertices 1,3 adjacent)
    @test vals[1,2,1,2] == -Inf  # 0101 (vertices 2,4 adjacent)
    @test vals[2,2,2,1] == -Inf  # 1110
    @test vals[2,1,2,2] == -Inf  # 1011
    @test vals[2,2,1,2] == -Inf  # 1101
    @test vals[1,2,2,2] == -Inf  # 0111
    @test vals[2,2,2,2] == -Inf  # 1111
end
