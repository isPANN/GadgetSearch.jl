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

function make_figure1_right_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 1, 4)
    add_edge!(g, 3, 4)
    add_edge!(g, 2, 4)
    return g
end

@testset "Figure 1 right graph: alpha tensor against Table 1" begin
    g = make_figure1_right_graph()
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

@testset "Figure 1 right graph: reduced alpha tensor against Table 1" begin
    g = make_figure1_right_graph()
    boundary = [1, 2, 3, 4]
    reduced = calculate_reduced_alpha_tensor(g, boundary)
    vals = content.(reduced)
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


# Figure 1 left graph G: right graph (R') with 4 additional internal vertices {5,6,7,8}
# representing the original gadget R before graph rewriting.
# Pins are unchanged: [1, 2, 3, 4]
function make_figure1_left_graph()
    g = SimpleGraph(8)
    # Only boundary-to-boundary edge: 1-3
    add_edge!(g, 1, 3)
    # Internal red vertices 5, 6, 7, 8 connected to boundary and each other
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

@testset "Figure 1 left graph: alpha tensor against Table 1" begin
    g = make_figure1_left_graph()
    boundary = [1, 2, 3, 4]
    alpha = calculate_alpha_tensor(g, boundary)
    vals = content.(alpha)

    # Verify against Table 1 from the paper (Figure 1 left graph R)
    # Indexing: (s1+1, s2+1, s3+1, s4+1) for boundary config s1s2s3s4
    @test vals[1,1,1,1] == 2.0   # 0000
    @test vals[1,1,1,2] == 3.0   # 0001
    @test vals[1,1,2,1] == 3.0   # 0010
    @test vals[1,2,1,1] == 3.0   # 0100
    @test vals[2,1,1,1] == 3.0   # 1000
    @test vals[1,1,2,2] == 3.0   # 0011
    @test vals[1,2,1,2] == 3.0   # 0101
    @test vals[1,2,2,1] == 4.0   # 0110
    @test vals[2,1,1,2] == 3.0   # 1001
    @test vals[2,1,2,1] == -Inf  # 1010 (vertices 1,3 adjacent)
    @test vals[2,2,1,1] == 4.0   # 1100
    @test vals[1,2,2,2] == 4.0   # 0111
    @test vals[2,1,2,2] == -Inf  # 1011
    @test vals[2,2,1,2] == 4.0   # 1101
    @test vals[2,2,2,1] == -Inf  # 1110
    @test vals[2,2,2,2] == -Inf  # 1111
end

@testset "Figure 1 left graph: reduced alpha tensor against Table 1" begin
    g = make_figure1_left_graph()
    boundary = [1, 2, 3, 4]
    reduced = calculate_reduced_alpha_tensor(g, boundary)
    vals = content.(reduced)

    # Verify against Table 1 from the paper (α̃(R) column)
    # Indexing: (s1+1, s2+1, s3+1, s4+1) for boundary config s1s2s3s4
    @test vals[1,1,1,1] == 2.0   # 0000
    @test vals[1,1,1,2] == 3.0   # 0001
    @test vals[1,1,2,1] == 3.0   # 0010
    @test vals[1,2,1,1] == 3.0   # 0100
    @test vals[2,1,1,1] == 3.0   # 1000
    @test vals[1,1,2,2] == -Inf  # 0011
    @test vals[1,2,1,2] == -Inf  # 0101
    @test vals[1,2,2,1] == 4.0   # 0110
    @test vals[2,1,1,2] == -Inf  # 1001
    @test vals[2,1,2,1] == -Inf  # 1010
    @test vals[2,2,1,1] == 4.0   # 1100
    @test vals[1,2,2,2] == -Inf  # 0111
    @test vals[2,1,2,2] == -Inf  # 1011
    @test vals[2,2,1,2] == -Inf  # 1101
    @test vals[2,2,2,1] == -Inf  # 1110
    @test vals[2,2,2,2] == -Inf  # 1111
end

@testset "is_diff_by_constant" begin
    # Tensors that differ by a constant (+2)
    t1 = [3.0, 4.0, -Inf, 5.0]
    t2 = [1.0, 2.0, -Inf, 3.0]
    valid, c = is_diff_by_constant(t1, t2)
    @test valid == true
    @test c == 2.0

    # Tensors with mismatched -Inf positions → not valid
    t3 = [3.0, -Inf, -Inf, 5.0]
    t4 = [1.0,  2.0, -Inf, 3.0]
    valid2, _ = is_diff_by_constant(t3, t4)
    @test valid2 == false

    # Tensors with inconsistent differences → not valid
    t5 = [3.0, 4.0, -Inf, 6.0]
    t6 = [1.0, 2.0, -Inf, 3.0]
    valid3, _ = is_diff_by_constant(t5, t6)
    @test valid3 == false

    # All -Inf → valid (trivially, constant is NaN)
    t7 = [-Inf, -Inf]
    t8 = [-Inf, -Inf]
    valid4, _ = is_diff_by_constant(t7, t8)
    @test valid4 == true
end

# BATOIDEA: 11-vertex unit disk replacement for CROSS (Figure 6 in the paper)
# Pins: {1, 2, 3, 4}, internal: {5, 6, 7, 8, 9, 10, 11}
function make_batoidea_graph()
    g = SimpleGraph(11)
    add_edge!(g, 1, 5);  add_edge!(g, 1, 9)
    add_edge!(g, 2, 5);  add_edge!(g, 2, 6);  add_edge!(g, 2, 7)  # 2-7, not 2-8
    add_edge!(g, 3, 8)
    add_edge!(g, 4, 9);  add_edge!(g, 4, 10); add_edge!(g, 4, 11)
    add_edge!(g, 5, 6);  add_edge!(g, 5, 9);  add_edge!(g, 5, 10)
    add_edge!(g, 6, 7);  add_edge!(g, 6, 9);  add_edge!(g, 6, 10); add_edge!(g, 6, 11)
    add_edge!(g, 7, 8);  add_edge!(g, 7, 10); add_edge!(g, 7, 11)
    add_edge!(g, 8, 11)
    add_edge!(g, 9, 10)
    add_edge!(g, 10, 11)
    return g
end

@testset "CROSS and BATOIDEA: reduced alpha tensors differ by constant (Theorem 3.7)" begin
    pins = [1, 2, 3, 4]
    valid, c = is_gadget_replacement(make_cross_graph(), make_batoidea_graph(), pins, pins)
    @test valid == true
    @test c == 2.0  # BATOIDEA has 7 internal vertices contributing offset +2
end

# Equation (3.2): CROSS + EDGE → PIRAMID, constant difference = -1
# Left graph (CROSS + EDGE): pins {1,2,3,4}, internal {5,6}
function make_cross_edge_graph()
    g = SimpleGraph(6)
    add_edge!(g, 1, 2)
    add_edge!(g, 1, 5)
    add_edge!(g, 2, 6)
    add_edge!(g, 3, 5)
    add_edge!(g, 6, 4)
    return g
end

# Right graph (PIRAMID): pins {1,2,3,4}, internal {5}
function make_piramid_graph()
    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 1)
    add_edge!(g, 1, 5)
    add_edge!(g, 2, 5)
    add_edge!(g, 3, 5)
    add_edge!(g, 4, 5)
    return g
end

@testset "CROSS+EDGE and PIRAMID: reduced alpha tensors differ by constant (Eq. 3.2)" begin
    pins = [1, 2, 3, 4]
    valid, c = is_gadget_replacement(make_cross_edge_graph(), make_piramid_graph(), pins, pins)
    @test valid == true
    @test c == -1.0  # PIRAMID has fewer internal vertices, offset = -1
end