"""
    _has_shortg() -> Bool

Check whether the `shortg` binary (from Nauty/Traces) is available on PATH.
"""
_has_shortg() = Sys.which("shortg") !== nothing

function _call_shortg(temp_path::String, mapping_file::String)
    run(pipeline(`shortg -v -u $(temp_path)`, stderr=mapping_file))
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

function _process_and_save_graphs(
    results::Vector{Tuple{SimpleGraph{T}, Vector{Tuple{Float64, Float64}}}},
    path::String,
) where T
    if !_has_shortg()
        @warn "`shortg` not found in PATH; saving graphs without deduplication."
        save_graph(results, path)
        return path
    end

    mapping_file = tempname()
    temp_path    = tempname()

    original_coords = getindex.(results, 2)
    save_graph(results, temp_path; g6_only=true)

    try
        _call_shortg(temp_path, mapping_file)
        canon2orig, _ = _parse_shortg_mapping(mapping_file)
        _write_original_representatives(temp_path, canon2orig, original_coords, path)
    finally
        isfile(mapping_file) && rm(mapping_file; force=true)
        isfile(temp_path)    && rm(temp_path; force=true)
    end
    return path
end

function _write_original_representatives(
    original_file::String,
    canon2orig::Dict{Int, Vector{Int}},
    orig_coords::Vector{Vector{Tuple{Float64, Float64}}},
    output_file::String,
)
    original_lines = readlines(original_file)
    to_be_written = Vector{Tuple{AbstractString, Vector{Tuple{Float64, Float64}}}}()
    for canon_line in sort(collect(keys(canon2orig)))
        orig_line = canon2orig[canon_line][1]
        push!(to_be_written, (original_lines[orig_line], orig_coords[orig_line]))
    end
    save_graph(to_be_written, output_file)
end
