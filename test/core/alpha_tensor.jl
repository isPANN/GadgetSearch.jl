using GadgetSearch
using Graphs
using GenericTensorNetworks: content, Tropical
using Test

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

function make_test_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 1, 4)
    add_edge!(g, 3, 4)
    add_edge!(g, 2, 4)
    return g
end

@testset "make_test_graph" begin
    g = make_test_graph()
    boundary = [1, 2, 3, 4]
    alpha = calculate_alpha_tensor(g, boundary)
    vals = content.(alpha)
    @test vals[1,1,1,1] == 0.0
    @test vals[1,1,1,2] == 1.0
    @test vals[1,1,2,1] == 1.0
    @test vals[1,2,1,1] == 1.0
    @test vals[2,1,1,1] == 1.0
    @test vals[1,1,2,2] == -Inf
    @test vals[1,2,1,2] == -Inf
    @test vals[1,2,2,1] == 2.0
    @test vals[2,1,1,2] == -Inf
    @test vals[2,1,2,1] == -Inf
    @test vals[2,2,1,1] == 2.0
    @test vals[1,2,2,2] == -Inf
    @test vals[2,1,2,2] == -Inf
    @test vals[2,2,1,2] == -Inf
    @test vals[2,2,2,1] == -Inf
    @test vals[2,2,2,2] == -Inf
end

# Print alpha tensor values for manual verification against Table 1
let
    g = make_test_graph()
    alpha = calculate_alpha_tensor(g, [1, 2, 3, 4])
    vals = content.(alpha)
    println("\nmake_test_graph alpha tensor (boundary=[1,2,3,4]):")
    for s1 in 0:1, s2 in 0:1, s3 in 0:1, s4 in 0:1
        println("  $s1$s2$s3$s4 => ", vals[s1+1, s2+1, s3+1, s4+1])
    end
end

# Figure 1 left graph G: right graph (R') with 4 additional internal vertices {5,6,7,8}
# representing the original gadget R before graph rewriting.
# Pins are unchanged: [1, 2, 3, 4]
function make_figure1_left_graph()
    g = SimpleGraph(8)
    # Edges from the right graph (shared boundary structure)
    add_edge!(g, 1, 3)
    add_edge!(g, 1, 4)
    add_edge!(g, 3, 4)
    add_edge!(g, 2, 4)
    # Internal red vertices 5, 6, 7, 8
    add_edge!(g, 5, 1)
    add_edge!(g, 5, 2)
    add_edge!(g, 5, 6)
    add_edge!(g, 6, 7)
    add_edge!(g, 7, 2)
    add_edge!(g, 7, 3)
    add_edge!(g, 7, 8)
    add_edge!(g, 5, 8)
    add_edge!(g, 8, 4)
    return g
end

@testset "Figure 1: graph rewriting preserves alpha tensor" begin
    g_left  = make_figure1_left_graph()
    g_right = make_test_graph()
    boundary = [1, 2, 3, 4]

    alpha_left  = calculate_alpha_tensor(g_left,  boundary)
    alpha_right = calculate_alpha_tensor(g_right, boundary)

    @test content.(alpha_left) == content.(alpha_right)
end