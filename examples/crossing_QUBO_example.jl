# # Search for Crossing Gadget (QUBO) on Triangular Lattice
#
# A crossing gadget has 4 pins representing two independent wires:
#   Wire 1: pin1 → pin3 (copy)
#   Wire 2: pin2 → pin4 (copy)
# Ground states: all configs where pin1==pin3 AND pin2==pin4

using GadgetSearch
using HiGHS

# Crossing constraint: two independent copy wires
crossing_constraint = TruthTableConstraint(Bool[
    0 0 0 0;   # both wires carry 0
    0 1 0 1;   # wire1=0, wire2=1
    1 0 1 0;   # wire1=1, wire2=0
    1 1 1 1    # both wires carry 1
])

# Generate complete graphs on triangular lattice
# Increase grid size if no results found (e.g., 4×4, 5×5)
grid_size = 3
dataset_path = joinpath(pkgdir(GadgetSearch), "examples", "crossing_qubo_dataset.jsonl")
generate_full_grid_graph(Triangular(), grid_size, grid_size; path=dataset_path)

dataloader = GraphLoader(dataset_path)

# Search for crossing gadgets
results, failed = search_gadgets(
    QUBOModel,
    dataloader,
    [crossing_constraint];
    optimizer=HiGHS.Optimizer,
    allow_defect=true,
    max_result_num=10,
    max_samples=5000,
    check_connectivity=true
)

@show failed

# Display and visualize results
if !isempty(results[1])
    for (i, gadget) in enumerate(results[1])
        println("\n===== Crossing Gadget #$i =====")
        println("Pins: ", gadget.pins)
        println("Vertex weights: ", gadget.vertex_weights)
        println("Edge weights: ", gadget.edge_weights)
        println("Edge list: ", gadget.edge_list)
        outpath = joinpath(pkgdir(GadgetSearch), "examples", "crossing_qubo_gadget_$i.png")
        GadgetSearch.plot_gadget(gadget, outpath;
            show_weights=true, show_edge_weights=true, round_weights=true)
        println("Saved to $outpath")
    end
    # Validate first result
    println("\n===== Validation (QUBO) =====")
    println(check_gadget_qubo(results[1][1]; _return_info=true))
else
    println("No crossing gadget found. Try increasing grid_size.")
end

GadgetSearch.clear_cache!()
