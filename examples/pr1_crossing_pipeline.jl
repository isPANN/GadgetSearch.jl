# # PR1 Crossing Pipeline (line-by-line runnable)
#
# This script demonstrates two parts of the unweighted crossing workflow:
# 1) logical flip utilities
# 2) search_unweighted_gadgets
#
# It is intentionally organized into small, single-purpose functions so users can
# execute each section line by line in the REPL and inspect output immediately.

using GadgetSearch
using Graphs
using GenericTensorNetworks: content

const OUTPUT_DIR = pkgdir(GadgetSearch, "examples", "pr1_pipeline_output")

function build_cross_graph()
    g = SimpleGraph(4)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 4)
    return g
end

function print_graph_summary(name::AbstractString, g::SimpleGraph{Int})
    println("[$name] vertices=$(nv(g)), edges=$(ne(g)), edge_list=$(collect(edges(g)))")
end

function ensure_output_dir()
    mkpath(OUTPUT_DIR)
    println("Output directory: $OUTPUT_DIR")
end

function plot_canonical_crossing(g::SimpleGraph{Int})
    path = joinpath(OUTPUT_DIR, "crossing_canonical.svg")
    GadgetSearch.plot_graph(g, path)
    println("Saved canonical crossing plot -> $path")
    return path
end

function plot_found_replacement(g::SimpleGraph{Int})
    path = joinpath(OUTPUT_DIR, "crossing_search_match.svg")
    GadgetSearch.plot_graph(g, path)
    println("Saved search-match plot -> $path")
    return path
end

function run_flip_demo(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int})
    println("\n=== Module: flip ===")
    reduced = Float64.(content.(calculate_reduced_alpha_tensor(target_graph, target_boundary)))
    patterns = generate_flip_patterns(length(target_boundary))
    println("Flip patterns: $(length(patterns))")
    for (mask, desc) in patterns
        flipped = apply_flip_to_tensor(reduced, mask)
        finite_deltas = [f - b for (f, b) in zip(vec(flipped), vec(reduced)) if isfinite(f) && isfinite(b)]
        example_delta = isempty(finite_deltas) ? "n/a" : string(first(finite_deltas))
        println("  $desc, mask=$mask, finite_delta_example=$example_delta")
    end
end

function build_loader_from_candidates(candidate_graphs::Vector{SimpleGraph{Int}}, target_boundary::Vector{Int})
    isempty(candidate_graphs) && throw(ArgumentError("candidate_graphs must be non-empty"))
    dataset = GraphDataset(graph_to_g6.(candidate_graphs))
    return GraphLoader(dataset, pinset=target_boundary)
end

function run_search_demo(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int}, candidate_graphs::Vector{SimpleGraph{Int}})
    println("\n=== Module: search ===")
    loader = build_loader_from_candidates(candidate_graphs, target_boundary)
    results = search_unweighted_gadgets(
        target_graph,
        target_boundary,
        loader;
        include_logical_flips=true,
        max_results=10,
    )
    println("Search hits: $(length(results))")
    for (i, result) in enumerate(results)
        println("  hit[$i]: boundary=$(result.boundary_vertices), offset=$(result.constant_offset), vertices=$(nv(result.replacement_graph))")
    end
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("PR1 crossing pipeline demo (line-by-line friendly)")
    ensure_output_dir()

    target_graph = build_cross_graph()
    target_boundary = [1, 2, 3, 4]

    println("\n=== Base graph ===")
    print_graph_summary("canonical", target_graph)
    plot_canonical_crossing(target_graph)

    run_flip_demo(target_graph, target_boundary)
    candidates = [target_graph]
    results = run_search_demo(target_graph, target_boundary, candidates)

    if !isempty(results)
        plot_found_replacement(results[1].replacement_graph)
    end

    println("\nDone. You can now inspect outputs in: $OUTPUT_DIR")
    println("Tip: execute function calls above one section at a time in REPL.")
end
