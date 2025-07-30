struct GraphDataset
    g6codes::Vector{SubString{String}}
    layouts::Vector{Union{Nothing, Vector{Tuple{Float64, Float64}}}}
    n::Int
end

mutable struct GraphLoader
    # Graph dataset and parsing
    dataset::GraphDataset
    pinset::Union{Nothing, Vector{Int}}
    temp_bitvec::BitVector

    # Caching parsed graphs (optional)
    enable_cache::Bool
    max_cached::Int
    cachepath::Union{Nothing, String}
    parsed::Dict{String, SimpleGraph{Int}}
    parse_order::Vector{String}
end

struct LayoutAccessor
    cds::GraphLoader
end

function GraphDataset(path::String)
    raw = read(path, String)
    codes = SubString{String}[]
    layouts = Union{Nothing, Vector{Tuple{Float64, Float64}}}[]
    i = 1
    line_start = 1
    len = lastindex(raw)

    while line_start <= len
        line_end = findnext(==('\n'), raw, line_start)
        line_end === nothing && (line_end = len)

        # Trim leading/trailing space
        start = line_start
        stop  = line_end
        while start <= stop && isspace(raw[start])
            start += 1
        end
        while stop >= start && isspace(raw[stop])
            stop -= 1
        end

        if start <= stop
            line = @view raw[start:stop]
            sep = findfirst(isspace, line)
            g6 = sep === nothing ? line : @view line[1:sep-1]
            push!(codes, g6)

            if sep === nothing
                push!(layouts, nothing)
            else
                coords_str = String(line[sep+1:end])
                coord_parts = split(coords_str, ';')
                parsed_coords = Tuple{Float64, Float64}[]
                valid = true
                for part in coord_parts
                    m = match(r"\(?\s*([-+eE.\d]+)\s*,\s*([-+eE.\d]+)\s*\)?", part)
                    if m === nothing
                        valid = false
                        break
                    end
                    x = tryparse(Float64, m.captures[1])
                    y = tryparse(Float64, m.captures[2])
                    if x === nothing || y === nothing
                        valid = false
                        break
                    end
                    push!(parsed_coords, (x, y))
                end
                push!(layouts, valid ? parsed_coords : nothing)
            end

            i += 1
        end

        line_start = line_end + 1
    end

    return GraphDataset(codes, layouts, length(codes))
end


function GraphLoader(path::String; cachepath::Union{Nothing, String}=nothing, max_cached::Int=10_000, enable_cache::Bool=false, pinset::Union{Nothing, Vector{Int}}=nothing)
    ds = GraphDataset(path)
    parsed_cache = Dict{String, SimpleGraph{Int}}()
    bitvec = BitVector(undef, 0)

    if cachepath !== nothing && isfile(cachepath)
        @info "Loading parsed graph cache from $cachepath"
        try
            loaded = deserialize(cachepath)
            if isa(loaded, Dict{String, SimpleGraph{Int}})
                parsed_cache = loaded
            else
                @warn "Cache file format unexpected, ignoring cache."
            end
        catch e
            @warn "Failed to load cache file: $e"
        end
    end

    return GraphLoader(ds, pinset, bitvec, enable_cache, max_cached, cachepath, parsed_cache, String[])
end


function Base.getindex(cds::GraphLoader, key::String)
    idx = tryparse(Int, key)
    if idx === nothing
        error("GraphLoader: String key \"$key\" is not a valid integer index.")
    end
    return cds[idx]
end

function Base.getindex(cds::GraphLoader, idx::Int)
    key = string(idx)
    if cds.enable_cache && haskey(cds.parsed, key)
        return cds.parsed[key]
    else
        g6 = cds.dataset.g6codes[idx]
        g = _parse_g6_string(g6, cds.temp_bitvec)

        if cds.enable_cache
            cds.parsed[key] = g
            push!(cds.parse_order, key)
            if length(cds.parsed) > cds.max_cached
                oldest = popfirst!(cds.parse_order)
                delete!(cds.parsed, oldest)
            end
        end
        return g
    end
end

Base.getindex(l::LayoutAccessor, key::String) = getlayout(l.cds, key)

function Base.getindex(l::LayoutAccessor, key::Int)
    return getlayout(l.cds, string(key))
end

function getlayout(cds::GraphLoader, key::String)
    idx = parse(Int, key)
    return cds.dataset.layouts[idx]
end

Base.keys(cds::GraphLoader) = 1:cds.dataset.n
Base.length(cds::GraphLoader) = cds.dataset.n

# Allow property-like access for layout: cds.layout[key]
Base.getproperty(cds::GraphLoader, name::Symbol) = name === :layout ? LayoutAccessor(cds) : getfield(cds, name)

function save_cache(cds::GraphLoader)
    if cds.cachepath === nothing
        @warn "No cache path provided; cannot save cache."
        return
    end
    @info "Saving parsed graph cache to $(cds.cachepath)"
    serialize(cds.cachepath, cds.parsed)
end

function Base.show(io::IO, loader::GraphLoader)
    nkeys = length(loader.dataset.g6codes)
    print(io, "GraphLoader with $nkeys graphs",
          loader.enable_cache ? ", cache enabled" : ", no cache",
          "")
end