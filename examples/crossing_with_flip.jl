# Search for crossing gadgets on the triangular lattice with
# fixed pin roles and inner-subgraph isomorphism deduplication.

using Dates
using GadgetSearch
using Graphs
using JSON3
using ProgressMeter

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "crossing_flip")
const LOG_PATH = joinpath(OUTPUT_DIR, "search.log")
const RESULTS_PATH = joinpath(OUTPUT_DIR, "results.json")

mkpath(OUTPUT_DIR)

function log_msg(msg::AbstractString)
    line = "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg"
    println(line)
    open(LOG_PATH, "a") do io
        println(io, line)
    end
end

function save_results(results::Vector{MultiTargetResult}, target_descs::Vector{String})
    payload = map(results) do r
        g = r.gadget.replacement_graph
        Dict(
            "target_index" => r.target_index,
            "target_desc" => target_descs[r.target_index],
            "nv" => nv(g),
            "ne" => ne(g),
            "boundary_vertices" => r.gadget.boundary_vertices,
            "constant_offset" => r.gadget.constant_offset,
            "pos" => isnothing(r.gadget.pos) ? nothing : [[p[1], p[2]] for p in r.gadget.pos],
            "edges" => [[src(e), dst(e)] for e in edges(g)],
        )
    end
    open(RESULTS_PATH, "w") do io
        write(io, JSON3.write(payload; pretty=true))
    end
end

function build_triangular_graph(points::Vector{Tuple{Int,Int}})
    g = SimpleGraph(length(points))
    for a in 1:length(points), b in a+1:length(points)
        ia, ja = points[a]
        ib, jb = points[b]
        GadgetSearch.triangular_adjacency(ia, ja, ib, jb) && add_edge!(g, a, b)
    end
    return g
end

base_targets = generate_extended_cross()
flip_patterns = generate_flip_patterns()
filter_fn, target_descs = make_flip_aware_multi_target_filter(base_targets, flip_patterns)

log_msg("Base targets: $(length(base_targets))")
log_msg("Flip patterns: $(length(flip_patterns))")
log_msg("Total target tensors: $(length(target_descs))")

# Only keep nx <= ny because swapping the grid axes gives an equivalent search.
grid_configs = [
    (4, 5, 6, 14),
    (4, 6, 7, 14),
    (5, 5, 7, 14),
    (5, 6, 8, 16),
    (6, 6, 9, 18),
]

results = MultiTargetResult[]
total_checked = 0

for (nx, ny, min_k, max_k) in grid_configs
    log_msg("Starting grid $(nx)x$(ny), k=$(min_k):$(max_k)")

    Mx = nx + 2
    Ny_ = ny + 2
    inner_grid = Tuple{Int,Int}[(x, y) for x in 2:Mx-1 for y in 2:Ny_-1]
    inner_phys = GadgetSearch.get_physical_positions(Triangular(), inner_grid)
    n_inner = length(inner_grid)
    actual_max = min(max_k, n_inner)

    pin1_cands = Tuple{Int,Int}[(Mx, y) for y in 2:Ny_-1]   # bottom
    pin2_cands = Tuple{Int,Int}[(x, Ny_) for x in 2:Mx-1]   # right
    pin3_cands = Tuple{Int,Int}[(1, y) for y in 2:Ny_-1]    # top
    pin4_cands = Tuple{Int,Int}[(x, 1) for x in 2:Mx-1]     # left
    n_pin_combos = length(pin1_cands) * length(pin2_cands) * length(pin3_cands) * length(pin4_cands)

    for k in min_k:actual_max
        n_total_k = binomial(n_inner, k)
        unique_subsets = dedup_inner_subsets(inner_grid, k)
        n_unique = length(unique_subsets)
        dedup_pct = n_total_k > 0 ? round(100 * (1 - n_unique / n_total_k); digits=1) : 0.0
        log_msg("  k=$k: $n_total_k -> $n_unique unique inner graphs ($(dedup_pct)% reduced)")

        @showprogress "$(nx)x$(ny) k=$k" for p1 in pin1_cands, p2 in pin2_cands, p3 in pin3_cands, p4 in pin4_cands
            pin_grid = Tuple{Int,Int}[p1, p2, p3, p4]
            pin_phys = GadgetSearch.get_physical_positions(Triangular(), pin_grid)

            for subset in unique_subsets
                all_grid = vcat(pin_grid, inner_grid[subset])
                all_phys = vcat(pin_phys, inner_phys[subset])
                candidate = build_triangular_graph(all_grid)

                total_checked += 1
                result = filter_fn(candidate, all_phys, [1, 2, 3, 4])
                result === nothing && continue

                push!(results, result)
                save_results(results, target_descs)
                log_msg("  Found result #$(length(results)) on $(nx)x$(ny), k=$k, target=$(target_descs[result.target_index])")
            end
        end

        log_msg("  finished k=$k over $n_pin_combos pin combinations")
    end
end

save_results(results, target_descs)

println("\nFound $(length(results)) crossing gadget replacements after checking $total_checked candidates.")
for (i, r) in enumerate(results)
    target_desc = target_descs[r.target_index]
    g = r.gadget.replacement_graph
    println("\nResult #$i:")
    println("  Target: $target_desc")
    println("  Graph: $(nv(g))v / $(ne(g))e")
    println("  Pins: $(r.gadget.boundary_vertices)")
    println("  Offset: $(r.gadget.constant_offset)")
end


