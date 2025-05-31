struct GraphDataset
    keys::Vector{SubString{String}}
    g6codes::Vector{SubString{String}}
    raw::String
    key_to_index::Dict{String, Int}
end

mutable struct GraphLoader
    # Graph dataset and parsing
    dataset::GraphDataset
    temp_bitvec::BitVector

    # Layout info (optional)
    layoutfile::Union{Nothing, String}
    layoutcache::Dict{String, Vector{Int}}

    # Caching parsed graphs
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
    keys = SubString{String}[]
    codes = SubString{String}[]
    key_to_index = Dict{String, Int}()

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
            sep = findfirst(==(' '), line)
            if sep !== nothing
                key = @view line[1:sep-1]
                g6  = @view line[sep+1:end]
            else
                key = "graph$i"
                g6  = line
            end
            push!(keys, key)
            push!(codes, g6)
            key_to_index[key] = i
            i += 1
        end

        line_start = line_end + 1
    end

    return GraphDataset(keys, codes, raw, key_to_index)
end


function GraphLoader(path::String; cachepath::Union{Nothing, String}=nothing, max_cached::Int=10_000, enable_cache::Bool=false, layoutfile::Union{Nothing, String}=nothing)
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

    return GraphLoader(ds, bitvec, layoutfile, Dict{String, Vector{Int}}(), enable_cache, max_cached, cachepath, parsed_cache, String[])
end

GraphLoader(path::String) = GraphLoader(path; enable_cache=false)


function Base.getindex(cds::GraphLoader, key::String)
    if cds.enable_cache && haskey(cds.parsed, key)
        return cds.parsed[key]
    else
        idx = cds.dataset.key_to_index[key]
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

function getlayout(cds::GraphLoader, key::String)::Union{Nothing, Vector{Int}}
    if cds.layoutfile === nothing
        return nothing
    end
    if isempty(cds.layoutcache)
        cds.layoutcache = JSON3.read(cds.layoutfile, Dict{String, Vector{Int}})
    end
    # Use numeric ID as the JSON key
    if startswith(key, "graph")
        numkey = match(r"graph(\d+)", key)
        if numkey !== nothing
            return get(cds.layoutcache, numkey.captures[1], nothing)
        end
    end
    return nothing
end

Base.keys(cds::GraphLoader) = cds.dataset.keys
Base.length(cds::GraphLoader) = length(cds.dataset.keys)

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
    nkeys = length(loader.dataset.keys)
    print(io, "GraphLoader with $nkeys graphs",
          loader.enable_cache ? ", cache enabled" : ", no cache",
          loader.layoutfile !== nothing ? ", layout: ✔" : "")
end