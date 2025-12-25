```@meta
EditURL = "../../../examples/triangular_QUBO_example.jl"
```

# Search for QUBO Gadgets on Triangular Lattice

This example demonstrates how to search for **QUBO gadgets** using:
- **QUBOModel**: State space includes ALL 2^n binary configurations
- **Complete graphs**: All vertex pairs are connected
- **Vertex + Edge weights**: Energy E(σ) = Σᵢ hᵢσᵢ + Σᵢⱼ Jᵢⱼσᵢσⱼ

Notes:
- QUBO search is more general but computationally more expensive
- State constraints are specified directly as ground state strings

````@example triangular_QUBO_example
using GadgetSearch
using HiGHS
using Combinatorics
using FileIO, ImageShow
````

Define truth table constraints for QUBO
Each constraint specifies which pin configurations should be ground states
Format: TruthTableConstraint(BitMatrix) where each row is a ground state

````@example triangular_QUBO_example
constraints = [
    TruthTableConstraint(Bool[0 0 0; 0 1 1; 1 0 1; 1 1 1]),  # OR-like
    TruthTableConstraint(Bool[0 0 0; 0 1 0; 1 0 0; 1 1 1]),  # AND-like
    TruthTableConstraint(Bool[0 0 0; 0 1 1; 1 0 1; 1 1 0])   # XOR-like
]
````

Generate Complete Graph dataset on triangular lattice
All vertices are connected (no distance restriction like UDG)

````@example triangular_QUBO_example
generate_full_grid_graph(Triangular(), 2, 2; path=pkgdir(GadgetSearch, "examples", "qubo_dataset.g6"))

dataloader = GraphLoader(pkgdir(GadgetSearch, "examples", "qubo_dataset.g6"))
````

Search using QUBOModel explicitly
- QUBOModel uses ALL 2^n states (not just MIS)
- Both vertex weights (h) and edge weights (J) are optimized
- Objective minimizes sum of h and J

````@example triangular_QUBO_example
results, failed = search_gadgets(
    QUBOModel,              # Explicitly use QUBO model
    dataloader,
    constraints;
    optimizer=HiGHS.Optimizer,
    allow_defect=true,
    save_path=joinpath(pkgdir(GadgetSearch, "examples"), "triangular_QUBO_results.json"),
    max_result_num=5,
    max_samples=5000,
    check_connectivity=true
)

@show failed
println("================================")
println("    QUBO MODEL SEARCH RESULTS")
println("================================")
@info "Results: $(length(results)) constraints searched"
labels = ["OR-like", "AND-like", "XOR-like"]
for (i, label) in enumerate(labels)
    @info "Constraint '$label': $(length(results[i])) gadgets found"
end

@info "Cache statistics: $(GadgetSearch.get_cache_stats())"
GadgetSearch.clear_cache!()
````

## Examine found QUBO gadgets
QUBO gadgets have both vertex and edge weights

````@example triangular_QUBO_example
for (i, label) in enumerate(labels)
    if !isempty(results[i])
        gadget = results[i][1]
        @info """===== Found QUBO gadget for '$label' =====
        Constraint: $(gadget.ground_states)
        Pins: $(gadget.pins)
        Vertex weights (h): $(gadget.vertex_weights)
        Edge weights (J): $(gadget.edge_weights)
        Edge list: $(gadget.edge_list)
        """
        outpath = pkgdir(GadgetSearch, "examples", "qubo_gadget_$(label).png")
        GadgetSearch.plot_gadget(gadget, outpath;
            show_weights=true,           # Show vertex weights
            show_edge_weights=true,      # Show edge weights (QUBO specific)
            round_weights=true
        )
        @info "Saved to $outpath"
    end
end
````

## Check correctness using QUBO model
The check_gadget_qubo function validates using ALL states (not just MIS)

````@example triangular_QUBO_example
for (i, label) in enumerate(labels)
    if !isempty(results[i])
        gadget = results[i][1]
        println("\n===== Checking '$label' gadget (QUBO/Full state space) =====")
        println(check_gadget_qubo(gadget; _return_info=true))
    end
end
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

