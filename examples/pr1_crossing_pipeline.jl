# # PR1 Crossing Pipeline (line-by-line runnable)
#
# This script demonstrates the unweighted crossing search workflow.
#
# It is intentionally organized into small, single-purpose functions so users can
# execute each section line by line in the REPL and inspect output immediately.

using GadgetSearch
using Graphs
using Random

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

function run_search_demo(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int})
    println("\n=== Module: search ===")
    report = search_unweighted_gadgets(
        target_graph,
        target_boundary,
        Triangular();
        min_vertices=5,
        max_vertices=11,
        max_evaluations=2_000,
        beam_width=48,
        mutations_per_candidate=10,
        random_candidates_per_generation=16,
        max_results=4,
        rng=MersenneTwister(2026),
    )
    println("Evaluated: $(report.evaluated) candidates in $(report.generations) generations")
    println("Termination: $(report.termination_reason)")
    println("Best distance: mask mismatches=$(report.best_mask_mismatches), offset spread=$(report.best_offset_spread)")
    println("Search hits: $(length(report.gadgets))")
    for (i, result) in enumerate(report.gadgets)
        println("  hit[$i]: lattice=$(result.lattice), coordinates=$(result.lattice_coordinates)")
        println("          boundary=$(result.boundary_vertices), offset=$(result.constant_offset), vertices=$(nv(result.replacement_graph))")
    end
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("PR1 crossing pipeline demo (line-by-line friendly)")
    ensure_output_dir()

    target_graph = build_cross_graph()
    target_boundary = [1, 2, 3, 4]

    println("\n=== Base graph ===")
    print_graph_summary("canonical", target_graph)
    plot_canonical_crossing(target_graph)

    report = run_search_demo(target_graph, target_boundary)

    if !isempty(report.gadgets)
        plot_found_replacement(report.gadgets[1].replacement_graph)
    end

    println("\nDone. You can now inspect outputs in: $OUTPUT_DIR")
    println("Tip: execute function calls above one section at a time in REPL.")
end
