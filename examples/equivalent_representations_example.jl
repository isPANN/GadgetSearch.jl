using GadgetSearch
using Graphs

function cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    return g
end

graph = cross_graph()
boundary = [1, 2, 3, 4]
outpath = pkgdir(GadgetSearch, "examples", "equivalent_representations_cross.svg")

reprs = equivalent_representations(
    graph,
    boundary;
    max_added_vertices=4,
    preserve_boundary_roles=true,
)
println("Found $(length(reprs)) distinct equivalent representations.")

plot_crossing_equivalence_gallery(
    graph,
    boundary,
    outpath;
    max_added_vertices=4,
    preserve_boundary_roles=true,
    columns=2,
    panel_graph_size=180,
)

println("Saved visualization to $outpath")

