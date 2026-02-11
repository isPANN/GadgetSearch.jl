# Example: Using AlphaTensorMode for Gadget Search
#
# This example demonstrates how to use the α-tensor verification mode,
# which is based on the UnitDiskMapping.jl framework.
#
# Reference: Liu et al., "Computer-assisted gadget design and problem reduction
# of unweighted maximum independent set", PRX Quantum 4, 010316 (2023)

using GadgetSearch

# ============================================================================
# Setup: Generate Unit Disk Graph (King's Lattice)
# ============================================================================

println("Generating unit disk graph dataset...")
generate_full_grid_udg(Triangular(), 3, 3; path="dataset_alpha.g6")
loader = GraphLoader("dataset_alpha.g6")

# ============================================================================
# Define Constraints (Truth Tables)
# ============================================================================

# Example: Logic gates as truth tables
constraints = [
    # OR gate: output = 1 if any input is 1
    TruthTableConstraint(BitMatrix([
        0 0 0;  # 00 → 0
        1 0 1;  # 01 → 1
        0 1 1;  # 10 → 1
        1 1 1   # 11 → 1
    ])),
    
    # AND gate: output = 1 only if both inputs are 1
    TruthTableConstraint(BitMatrix([
        0 0 0;  # 00 → 0
        1 0 0;  # 01 → 0
        0 1 0;  # 10 → 0
        1 1 1   # 11 → 1
    ])),
    
    # XOR gate: output = 1 if inputs differ
    TruthTableConstraint(BitMatrix([
        0 0 0;  # 00 → 0
        1 0 1;  # 01 → 1
        0 1 1;  # 10 → 1
        1 1 0   # 11 → 0
    ]))
]

# ============================================================================
# Search Using AlphaTensorMode
# ============================================================================

println("\n" * "="^60)
println("Searching with AlphaTensorMode (α-tensor verification)")
println("="^60)

# Key difference: No optimizer needed!
results_alpha, failed_alpha = search_gadgets(
    AlphaTensorMode,      # Use α-tensor verification
    loader,
    constraints;
    pin_candidates=nothing,  # Auto-detect pins
    max_result_num=3,
    check_connectivity=true
)

println("\n✅ AlphaTensorMode search completed!")
println("Found gadgets for $(length(constraints) - length(failed_alpha)) / $(length(constraints)) constraints")
println("Failed constraints: $(length(failed_alpha))")

# Display results
for (i, gadgets) in enumerate(results_alpha)
    if !isempty(gadgets)
        println("\nConstraint $i: Found $(length(gadgets)) gadget(s)")
        for (j, g) in enumerate(gadgets[1:min(2, length(gadgets))])
            println("  Gadget $j:")
            println("    Vertices: $(g.num_vertices)")
            println("    Pins: $(g.pins)")
            println("    Weights: all uniform ($(g.vertex_weights[1]))")
        end
    end
end

# ============================================================================
# Comparison: AlphaTensorMode vs RydbergUnweightedModel
# ============================================================================

println("\n" * "="^60)
println("Comparison: AlphaTensorMode vs RydbergUnweightedModel")
println("="^60)

# RydbergUnweightedModel uses simpler cardinality check
results_unweighted, failed_unweighted = search_gadgets(
    RydbergUnweightedModel,
    loader,
    constraints;
    pin_candidates=nothing,
    max_result_num=3,
    check_connectivity=true
)

println("\n✅ RydbergUnweightedModel search completed!")
println("Found gadgets for $(length(constraints) - length(failed_unweighted)) / $(length(constraints)) constraints")

# ============================================================================
# Direct α-Tensor Computation (Advanced Usage)
# ============================================================================

println("\n" * "="^60)
println("Advanced: Direct α-Tensor Computation")
println("="^60)

using Graphs

# Create a simple example graph
g = SimpleGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)

# Designate pins (boundary vertices)
pins = [1, 5]

println("\nGraph: Linear chain with 5 vertices")
println("Pins: $pins (endpoints)")

# Compute reduced α-tensor
α = compute_reduced_alpha_tensor(g, pins)

println("\nReduced α-tensor:")
for (config, value) in sort(collect(α), by=first)
    # Decode configuration
    pin_states = [((config >> i) & 0x1) == 1 ? 1 : 0 for i in 0:(length(pins)-1)]
    println("  pins=$pin_states → max interior MIS size = $value")
end

# ============================================================================
# α-Tensor Equivalence Check
# ============================================================================

println("\n" * "="^60)
println("Advanced: α-Tensor Equivalence Check")
println("="^60)

# Two α-tensors that differ by a constant are equivalent
α1 = Dict(UInt32(0) => 3, UInt32(1) => 1, UInt32(2) => 2, UInt32(3) => 0)
α2 = Dict(UInt32(0) => 8, UInt32(1) => 6, UInt32(2) => 7, UInt32(3) => 5)

is_equiv, constant = check_alpha_equivalence(α1, α2)

println("\nα1 = $α1")
println("α2 = $α2")
println("Equivalent? $is_equiv")
if is_equiv
    println("Constant difference: $constant")
    println("This means: α2[config] = α1[config] + $constant for all configs")
end

# ============================================================================
# Summary: When to Use AlphaTensorMode
# ============================================================================

println("\n" * "="^60)
println("Summary: When to Use AlphaTensorMode")
println("="^60)

println("""
AlphaTensorMode is best when:
  ✅ You want mathematically rigorous gadget verification
  ✅ You need to verify gadget equivalence using reduced α-tensors
  ✅ You're working with unweighted systems (uniform vertex weights)
  ✅ You don't need edge weights
  ✅ You want to avoid dependency on optimization solvers

Comparison with other modes:
  
  RydbergModel:
    - Needs optimizer (HiGHS, Gurobi, etc.)
    - Optimizes vertex weights
    - More flexible (non-uniform weights)
    - State space: MIS
  
  RydbergUnweightedModel:
    - No optimizer needed
    - Simpler cardinality check
    - Uniform weights only
    - Faster but less rigorous
    - State space: MIS
  
  AlphaTensorMode:
    - No optimizer needed
    - α-tensor equivalence (most rigorous)
    - Uniform weights only
    - Based on theoretical framework
    - State space: MIS
  
  QUBOModel:
    - Needs optimizer
    - Optimizes vertex + edge weights
    - Most general
    - State space: All 2^n configurations
""")

println("\nExample completed! Check results.json for saved gadgets.")

# Clean up
rm("dataset_alpha.g6", force=true)

