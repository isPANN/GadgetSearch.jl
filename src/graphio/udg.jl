"""
    unit_disk_graph(locs::AbstractVector, unit::Real) -> SimpleGraph

Create a unit disk graph from given locations. Two vertices are connected 
if their Euclidean distance is less than the unit distance.

# Arguments
- `locs`: Vector of vertex locations
- `unit`: Unit distance threshold for edge creation

# Returns
- `SimpleGraph`: The constructed unit disk graph
"""
function unit_disk_graph(locs::AbstractVector, unit::Real)
    n = length(locs)
    g = SimpleGraph(n)
    for i = 1:n, j = i+1:n
        if sum(abs2, locs[i] .- locs[j]) < unit ^ 2
            add_edge!(g, i, j)
        end
    end
    return g
end

abstract type LatticeType end
struct Square <: LatticeType end
struct Triangular <: LatticeType end

get_radius(::Square) = 1.5
get_radius(::Triangular) = 1.1

get_physical_positions(::Square, pos::Vector{Tuple{Int, Int}}) = Vector{Tuple{Float64, Float64}}(pos)
function get_physical_positions(::Triangular, pos::Vector{Tuple{Int, Int}})
    h = sqrt(3) / 2
    return Vector{Tuple{Float64, Float64}}([(x + (isodd(y) ? 0.5 : 0.0), y * h) for (x, y) in pos])
end

function get_pin_positions(lattice::LatticeType, nx::Int, ny::Int)
    top_candidates    = get_physical_positions(lattice, Tuple{Int, Int}[(1, y)     for y in 2:ny+1])
    bottom_candidates = get_physical_positions(lattice, Tuple{Int, Int}[(nx+2, y)  for y in 2:ny+1])
    left_candidates   = get_physical_positions(lattice, Tuple{Int, Int}[(x, 1)     for x in 2:nx+1])
    right_candidates  = get_physical_positions(lattice, Tuple{Int, Int}[(x, ny+2)  for x in 2:nx+1])
    return top_candidates, bottom_candidates, left_candidates, right_candidates
end

function get_inner_points(lattice::LatticeType, nx::Int, ny::Int)
    original_points = Tuple{Int, Int}[(x, y) for x in 2:nx+1, y in 2:ny+1]
    return get_physical_positions(lattice, vec(original_points))
end

"""
    generate_full_grid_udg(lattice::LatticeType, nx::Int, ny::Int; path::String="udg.g6") -> String

Generate unit disk graphs on a grid lattice with boundary pins and save to file.

# Arguments  
- `lattice`: Type of lattice (Square or Triangular)
- `nx`: Number of inner points in x direction
- `ny`: Number of inner points in y direction
- `path`: Output file path for saving graphs (default: "udg.g6")

# Returns
- `String`: Path to the saved graph file

# Details
Generates all possible UDGs by placing pins on boundary positions and 
connecting vertices within unit distance on the specified lattice type.
"""
function generate_full_grid_udg(lattice::LatticeType, nx::Int, ny::Int; path::String="udg.g6")
    # get possible pin positions on the boundary
    top_candidates, bottom_candidates, left_candidates, right_candidates = get_pin_positions(lattice, nx, ny)

    # inner points are fixed
    inner_points = get_inner_points(lattice, nx, ny)

    radius = get_radius(lattice)

    results = Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}[]

    for top in top_candidates, bottom in bottom_candidates,
        left in left_candidates, right in right_candidates

        selected = vcat([top, right, bottom, left], inner_points)

        g = unit_disk_graph(selected, radius)

        push!(results, (g, selected))
    end
    @info "pinset in generated graphs: [1,2,3,4]"
    return _process_and_save_graphs(results, path)
end

function _call_shortg(temp_path::String, mapping_file::String)
    if Sys.which("shortg") === nothing
        # Make shortg optional: log and signal caller to fallback
        @warn "Optional tool `shortg` not found; skipping canonicalization/dedup."
        return false
    else
        @info "shortg found in PATH; running canonicalization"
    end
    run(pipeline(`shortg -v -u $(temp_path)`, stderr=mapping_file))
    return true
end

function _process_and_save_graphs(results::Vector{Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}}, path::String)
    # If shortg is unavailable, save directly without dedup/canonicalization.
    if Sys.which("shortg") === nothing
        @warn "`shortg` not found in PATH; saving graphs without deduplication."
        save_graph(results, path)
        return path
    end

    mapping_file = tempname()
    temp_path = tempname()

    original_coords = getindex.(results, 2)
    save_graph(results, temp_path; g6_only=true)

    ok = _call_shortg(temp_path, mapping_file)
    if !ok
        # Fallback path if shortg disappeared between checks or failed early
        @warn "Skipping canonicalization due to `shortg` issue; saving raw graphs."
        try
            isfile(mapping_file) && rm(mapping_file)
            isfile(temp_path) && rm(temp_path)
        catch
        end
        save_graph(results, path)
        return path
    end

    canonical_to_original, _ = _parse_shortg_mapping(mapping_file)
    
    _write_original_representatives(temp_path, canonical_to_original, original_coords, path)

    try
        isfile(mapping_file) && rm(mapping_file)
        isfile(temp_path) && rm(temp_path)
    catch
    end
    return path
end

function _parse_shortg_mapping(filepath::String)
    canonical_to_originals = Dict{Int, Vector{Int}}()
    original_to_canonical = Dict{Int, Int}()

    for line in eachline(filepath)
        if isempty(strip(line)) || startswith(line, '>') || startswith(line, 'Z')
            continue
        end

        if occursin(":", line)
            parts = split(line, ":")
            canonical = parse(Int, strip(parts[1]))
            originals = split(strip(parts[2]))
            for orig_str in originals
                original = parse(Int, orig_str)
                push!(get!(canonical_to_originals, canonical, Int[]), original)
                original_to_canonical[original] = canonical
            end
        end
    end

    return canonical_to_originals, original_to_canonical
end

function _write_original_representatives(
    original_file::String,
    canon2orig::Dict{Int, Vector{Int}},
    orig_coords::Vector{Vector{Tuple{Float64, Float64}}},
    output_file::String
)
    original_lines = readlines(original_file)
    to_be_written = Vector{Tuple{AbstractString, Vector{Tuple{Float64, Float64}}}}()
    open(output_file, "w") do io
        for canon_line in sort(collect(keys(canon2orig)))
            orig_line = canon2orig[canon_line][1]
            push!(to_be_written, (original_lines[orig_line], orig_coords[orig_line]))
        end
    end
    save_graph(to_be_written, output_file)
end
