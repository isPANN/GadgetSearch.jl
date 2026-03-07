# Search for Crossing Gadgets with Logical Flip

using GadgetSearch
using Graphs
using ProgressMeter

# Generate base targets and flip patterns
base_targets = generate_extended_cross()
flip_patterns = generate_flip_patterns()

println("Base targets: $(length(base_targets))")
println("Flip patterns: $(length(flip_patterns))")
println("Total combinations: $(length(base_targets) * length(flip_patterns))")

# Create flip-aware filter
filter_fn, target_descs = make_flip_aware_multi_target_filter(base_targets, flip_patterns)

println("\nTarget descriptions:")
for (i, desc) in enumerate(target_descs)
    println("  $i: $desc")
end

# Dataset configuration: (directory, filename) pairs
datasets = [
    # Triangular lattice UDGs (main search target!)
    ("triangular_udg_subsets", "tri_2x2_min3max4_direct4.g6"),
    ("triangular_udg_subsets", "tri_2x3_min3max6_direct4.g6"),
    ("triangular_udg_subsets", "tri_3x3_min3max9_direct4.g6"),
    ("triangular_udg_subsets", "tri_3x4_min3max12_direct4.g6"),
    ("triangular_udg_subsets", "tri_4x4_min4max16_direct4.g6"),
    ("triangular_udg_subsets", "tri_3x5_min3max15_direct4.g6"),
    ("triangular_udg_subsets", "tri_2x8_min3max16_direct4.g6"),
]

all_results = MultiTargetResult[]
result_sources = String[]  # track which dataset each result came from
total_candidates = 0

for (subdir, dataset_name) in datasets
    global total_candidates
    data_path = pkgdir(GadgetSearch, "data", subdir, dataset_name)
    !isfile(data_path) && (println("Skipping $subdir/$dataset_name (not found)"); continue)

    loader = GraphLoader(data_path; pinset=[1,2,3,4])
    total_candidates += length(loader)
    println("\nSearching $subdir/$dataset_name: $(length(loader)) candidates...")

    @showprogress for key in keys(loader)
        result = filter_fn(loader[key], loader.layout[key], loader.pinset)
        if result !== nothing
            push!(all_results, result)
            push!(result_sources, "$subdir/$dataset_name")
        end
    end
end

results = all_results

println("\nFound $(length(results)) crossing gadget replacements")


# Display results
for (i, r) in enumerate(results)
    target_desc = target_descs[r.target_index]
    g = r.gadget.replacement_graph
    source = result_sources[i]
    println("Result #$i: matched '$target_desc' | $(nv(g))v/$(ne(g))e | offset=$(r.gadget.constant_offset) | from $source")
end

# Generate Typst visualization
desktop = homedir() * "/Desktop"
typst_path = joinpath(desktop, "crossing_with_flip.typ")
pdf_path = joinpath(desktop, "crossing_with_flip.pdf")

open(typst_path, "w") do io
    println(io, "#set page(width: 210mm, height: 297mm, margin: 20mm)")
    println(io, "#set text(font: \"Arial\", size: 11pt)")
    println(io, "#align(center)[#text(size: 16pt, weight: \"bold\")[Crossing Gadgets with Logical Flip]]")
    println(io, "#v(1em)")
    println(io, "Searched $total_candidates candidates with $(length(target_descs)) target variants.\n")
    println(io, "Found $(length(results)) crossing gadget replacements.\n")

    for (i, r) in enumerate(results)
        target_desc = target_descs[r.target_index]
        g = r.gadget.replacement_graph

        source = result_sources[i]
        println(io, "== Result #$i")
        println(io, "- Target: `$target_desc`")
        println(io, "- Dataset: `$source`")
        println(io, "- Vertices: $(nv(g)), Edges: $(ne(g))")
        println(io, "- Constant offset: $(r.gadget.constant_offset)")
        println(io, "- Boundary: $(r.gadget.boundary_vertices)")

        svg_path = joinpath(desktop, "crossing_flip_$i.svg")
        GadgetSearch.plot_graph(g, svg_path; pos=r.gadget.pos, plot_size=300, margin=30)
        println(io, "\n#image(\"crossing_flip_$i.svg\", width: 60%)\n")
    end
end

println("\nTypst document: $typst_path")
run(`typst compile $typst_path $pdf_path`)
println("PDF generated: $pdf_path")


