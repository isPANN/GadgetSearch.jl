# Example: Using RydbergUnweightedModel for Gadget Search
#
# This example demonstrates how to use the α-tensor based unweighted
# verification mode for gadget search.
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
# Search Using RydbergUnweightedModel (α-tensor verification)
# ============================================================================

println("\n" * "="^60)
println("Searching with RydbergUnweightedModel (α-tensor verification)")
println("="^60)

# Key difference from RydbergModel: No optimizer needed!
# Instead of finding weights via ILP, this mode:
#   1. Computes the reduced α-tensor for the candidate graph
#   2. Determines ground state boundary configs (maximize total IS size)
#   3. Checks if they match the truth table
results, failed = search_gadgets(
    RydbergUnweightedModel,   # Uses α-tensor verification
    loader,
    constraints;
    pin_candidates=nothing,   # Auto-detect pins
    max_result_num=3,
    check_connectivity=true
)

println("\n✅ RydbergUnweightedModel search completed!")
println("Found gadgets for $(length(constraints) - length(failed)) / $(length(constraints)) constraints")
println("Failed constraints: $(length(failed))")

# Display results
for (i, gadgets) in enumerate(results)
    if !isempty(gadgets)
        println("\nConstraint $i: Found $(length(gadgets)) gadget(s)")
        for (j, g) in enumerate(gadgets[1:min(2, length(gadgets))])
            println("  Gadget $j:")
            println("    Vertices: $(nv(g.graph))")
            println("    Pins: $(g.pins)")
            println("    Weights: all uniform (1.0)")
        end
    end
end

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
    pin_states = [((config >> i) & 0x1) == 1 ? 1 : 0 for i in 0:(length(pins)-1)]
    total = value + count_ones(config)
    status = value == GadgetSearch.INFEASIBLE_ALPHA ? "INFEASIBLE" : "α=$value, total=$total"
    println("  pins=$pin_states → $status")
end

# Find ground configurations
ground_configs, max_total = find_ground_configs(α, length(pins))
println("\nGround configs (maximize total IS size = $max_total):")
for config in sort(collect(ground_configs))
    pin_states = [((config >> i) & 0x1) == 1 ? 1 : 0 for i in 0:(length(pins)-1)]
    println("  pins=$pin_states")
end

# ============================================================================
# α-Tensor Equivalence Check (for gadget replacement verification)
# ============================================================================

println("\n" * "="^60)
println("Advanced: α-Tensor Equivalence Check")
println("="^60)

# Two α-tensors that differ by a constant are equivalent
# This is used in gadget replacement: α̃(R') = α̃(P) + c
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
# Summary
# ============================================================================

println("\n" * "="^60)
println("Summary: RydbergUnweightedModel (α-tensor)")
println("="^60)

println("""
RydbergUnweightedModel uses the reduced α-tensor framework:
  ✅ No optimizer needed (pure combinatorial verification)
  ✅ Mathematically rigorous (based on α-tensor equivalence)
  ✅ Uniform weights (all vertices have weight 1)
  ✅ Works with King's lattice / unit disk graphs
  ✅ Supports gadget replacement verification

How it works:
  1. For each candidate graph + pin assignment:
     a. Compute reduced α-tensor: α[s] = max interior MIS size for boundary config s
     b. Find ground configs: maximize total(s) = α[s] + count_ones(s)
     c. Check if ground configs match the truth table
  2. If match found → valid gadget!

Comparison with other modes:
  
  RydbergModel:
    - Needs optimizer (HiGHS, Gurobi, etc.)
    - Finds non-uniform vertex weights via ILP
    - More flexible but slower
  
  RydbergUnweightedModel:
    - No optimizer needed
    - Uses α-tensor for verification
    - All weights uniform (= 1)
    - Based on the reduced α-tensor framework
  
  QUBOModel:
    - Needs optimizer
    - Finds both vertex + edge weights
    - Most general model
""")

println("Example completed!")

# Clean up
rm("dataset_alpha.g6", force=true)
