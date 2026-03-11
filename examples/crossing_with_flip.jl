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
const PROGRESS_PATH = joinpath(OUTPUT_DIR, "progress.json")
const RESUME_SEARCH = get(ENV, "CROSSING_RESUME", "1") != "0"

mkpath(OUTPUT_DIR)

function log_msg(msg::AbstractString)
    line = "[$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))] $msg"
    println(line)
    open(LOG_PATH, "a") do io
        println(io, line)
    end
end

coord_json(p::Tuple{Int,Int}) = [p[1], p[2]]
coord_json(p::Tuple{Float64,Float64}) = [p[1], p[2]]

function json_to_julia(x)
    if x isa JSON3.Object
        return Dict{String,Any}(String(k) => json_to_julia(v) for (k, v) in pairs(x))
    elseif x isa JSON3.Array
        return [json_to_julia(v) for v in x]
    else
        return x
    end
end

function stage_key(nx::Int, ny::Int, k::Int)
    return "$(nx)x$(ny):k=$(k)"
end

function load_results_state()
    if RESUME_SEARCH && isfile(RESULTS_PATH)
        parsed = JSON3.read(read(RESULTS_PATH, String))
        total_checked = hasproperty(parsed, :total_checked) ? Int(parsed.total_checked) : 0
        result_records = hasproperty(parsed, :results) ? [json_to_julia(r) for r in parsed.results] : Dict{String,Any}[]
        return result_records, total_checked
    end
    return Dict{String,Any}[], 0
end

function save_results(result_records::Vector{Dict{String,Any}}, total_checked::Int)
    open(RESULTS_PATH, "w") do io
        write(io, JSON3.write(Dict(
            "total_checked" => total_checked,
            "n_results" => length(result_records),
            "results" => result_records,
        ); pretty=true))
    end
end

function load_progress_state()
    if RESUME_SEARCH && isfile(PROGRESS_PATH)
        parsed = JSON3.read(read(PROGRESS_PATH, String))
        completed = hasproperty(parsed, :completed_stages) ? Set(String.(parsed.completed_stages)) : Set{String}()
        total_checked = hasproperty(parsed, :total_checked) ? Int(parsed.total_checked) : 0
        return completed, total_checked
    end
    return Set{String}(), 0
end

function save_progress(completed_stages::Set{String}, total_checked::Int)
    open(PROGRESS_PATH, "w") do io
        write(io, JSON3.write(Dict(
            "resume_enabled" => RESUME_SEARCH,
            "total_checked" => total_checked,
            "completed_stages" => sort!(collect(completed_stages)),
            "updated_at" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        ); pretty=true))
    end
end

function make_result_record(result, target_descs::Vector{String},
                            metadata::Dict{String,Any})
    g = result.gadget.replacement_graph
    return merge(Dict(
        "target_index" => result.target_index,
        "target_desc" => target_descs[result.target_index],
        "nv" => nv(g),
        "ne" => ne(g),
        "boundary_vertices" => result.gadget.boundary_vertices,
        "constant_offset" => result.gadget.constant_offset,
        "pos" => isnothing(result.gadget.pos) ? nothing : [coord_json(p) for p in result.gadget.pos],
        "edges" => [[src(e), dst(e)] for e in edges(g)],
    ), metadata)
end

function maybe_resume_message(result_records::Vector{Dict{String,Any}},
                              completed_stages::Set{String},
                              total_checked::Int)
    if RESUME_SEARCH && (!isempty(result_records) || !isempty(completed_stages) || total_checked > 0)
        log_msg("Resuming previous search state: $(length(result_records)) saved results, $(length(completed_stages)) completed stages, $total_checked checked candidates")
    elseif !RESUME_SEARCH
        log_msg("Resume disabled by CROSSING_RESUME=0; starting fresh run")
    end
end

"""Non-equivalent flip patterns for 4-pin CROSS graph (considering symmetry)."""
function generate_flip_patterns()
    return [
        (Int[], "no-flip"),
        ([1], "flip-pin1"),
        ([1,2], "flip-pin1-2"),
        ([1,3], "flip-pin1-3"),
        ([1,2,3,4], "flip-all")
    ]
end

"""CROSS variants with inserted nodes on edges."""
function generate_extended_cross()
    variants = Tuple{SimpleGraph{Int}, Vector{Int}, String}[]

    # Base CROSS
    cross = SimpleGraph(4)
    add_edge!(cross, 1, 3)
    add_edge!(cross, 2, 4)
    push!(variants, (cross, [1,2,3,4], "base"))

    # Variant 1: Insert node 5 on edge 1-3
    g1 = SimpleGraph(5)
    add_edge!(g1, 1, 5); add_edge!(g1, 5, 3)
    add_edge!(g1, 2, 4)
    push!(variants, (g1, [1,2,3,4], "ext5v-13"))

    # Variant 2: Insert node 5 on edge 2-4
    g2 = SimpleGraph(5)
    add_edge!(g2, 1, 3)
    add_edge!(g2, 2, 5); add_edge!(g2, 5, 4)
    push!(variants, (g2, [1,2,3,4], "ext5v-24"))

    # Variant 3: Insert nodes on both edges
    g3 = SimpleGraph(6)
    add_edge!(g3, 1, 5); add_edge!(g3, 5, 3)
    add_edge!(g3, 2, 6); add_edge!(g3, 6, 4)
    push!(variants, (g3, [1,2,3,4], "ext6v-both"))

    return variants
end

base_targets = generate_extended_cross()
flip_patterns = generate_flip_patterns()
filter_fn, target_descs = make_flip_aware_multi_target_filter(
    base_targets,
    flip_patterns;
    permute_pins=false,
)

log_msg("Base targets: $(length(base_targets))")
log_msg("Flip patterns: $(length(flip_patterns))")
log_msg("Total target tensors: $(length(target_descs))")
log_msg("Pin permutation expansion: disabled")

# Only keep nx <= ny because swapping the grid axes gives an equivalent search.
grid_configs = [
    (4, 5, 6, 14),
    (4, 6, 7, 14),
    (5, 5, 7, 14),
    (5, 6, 8, 16),
    (6, 6, 9, 18),
]

result_records, total_checked_results = load_results_state()
completed_stages, total_checked_progress = load_progress_state()
total_checked = max(total_checked_results, total_checked_progress)
maybe_resume_message(result_records, completed_stages, total_checked)

for (nx, ny, min_k, max_k) in grid_configs
    log_msg("Starting grid $(nx)x$(ny), k=$(min_k):$(max_k)")

    Mx = nx + 2
    My = ny + 2
    inner_grid = Tuple{Int,Int}[(x, y) for x in 2:Mx-1 for y in 2:My-1]
    inner_phys = GadgetSearch.get_physical_positions(Triangular(), inner_grid)
    n_inner = length(inner_grid)
    actual_max = min(max_k, n_inner)

    pin1_cands = Tuple{Int,Int}[(Mx, y) for y in 2:My-1]   # bottom
    pin2_cands = Tuple{Int,Int}[(x, My) for x in 2:Mx-1]   # right
    pin3_cands = Tuple{Int,Int}[(1, y) for y in 2:My-1]    # top
    pin4_cands = Tuple{Int,Int}[(x, 1) for x in 2:Mx-1]     # left
    n_pin_combos = length(pin1_cands) * length(pin2_cands) * length(pin3_cands) * length(pin4_cands)

    for k in min_k:actual_max
        stage = stage_key(nx, ny, k)
        if stage in completed_stages
            log_msg("  skipping completed stage $stage")
            continue
        end

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

                metadata = Dict{String,Any}(
                    "grid" => [nx, ny],
                    "k" => k,
                    "pin_grid" => [coord_json(p) for p in pin_grid],
                    "pin_phys" => [coord_json(p) for p in pin_phys],
                    "inner_subset_indices" => collect(subset),
                    "inner_subset_grid" => [coord_json(inner_grid[i]) for i in subset],
                    "candidate_index" => total_checked,
                )
                push!(result_records, make_result_record(result, target_descs, metadata))
                save_results(result_records, total_checked)
                log_msg("  Found result #$(length(result_records)) on $(nx)x$(ny), k=$k, target=$(target_descs[result.target_index])")
            end
        end

        push!(completed_stages, stage)
        save_progress(completed_stages, total_checked)
        save_results(result_records, total_checked)
        log_msg("  finished k=$k over $n_pin_combos pin combinations")
    end
end

save_progress(completed_stages, total_checked)
save_results(result_records, total_checked)

println("\nFound $(length(result_records)) crossing gadget replacements after checking $total_checked candidates.")
for (i, r) in enumerate(result_records)
    println("\nResult #$i:")
    println("  Target: $(r["target_desc"])")
    println("  Graph: $(r["nv"])v / $(r["ne"])e")
    println("  Pins: $(r["boundary_vertices"])")
    println("  Offset: $(r["constant_offset"])")
end


