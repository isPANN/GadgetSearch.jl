"""
    get_gids(s::String)

Extract `graph_id` from a key of a `graph_dict`, e.g. "graph1" -> 1.

# Arguments
- `s::String`: The string containing a graph identifier

# Returns
- The extracted integer ID
"""
function get_gids(s::String)
    return parse.(Int, collect(eachmatch(r"\d+", s)) .|> x -> x.match)[1]
end

"""
    read_g6_graphs(path::String)::Dict{String, SimpleGraph{Int}}

Read graphs in Graph6 format from a file.

# Arguments
- `path::String`: Path to the file containing Graph6 encoded graphs

# Returns
- A dictionary mapping graph keys, like `graph1`, to SimpleGraph objects
"""
function read_g6_graphs(path::String)::Dict{String, SimpleGraph{Int}}
    graph_dict = Dict{String, SimpleGraph{Int}}()

    open(path, "r") do io
        line_num = 1
        temp_bitvec = BitVector(undef, 0)
        for line in eachline(io)
            line = strip(line)
            isempty(line) && continue
        
            parts = split(line, ' ', limit=2)
            key, g6_code = length(parts) == 2 ? (parts[1], strip(parts[2])) : ("graph$(line_num)", parts[1])
            
            try
                graph = GraphIO.Graph6._g6StringToGraph(g6_code)
                graph = _parse_g6_string(g6_code, temp_bitvec)
                graph_dict[key] = graph
            catch e
                @warn "Failed to parse graph at line $line_num: $e"
            end
        
            length(parts) == 1 && (line_num += 1)
        end
    end

    return graph_dict
end

function read_g6_graph(path::String, id::Int)
    graph_dict = read_g6_graphs(path)
    return graph_dict["graph$(id)"]
end

"""
    each_g6_graph(path::String; buffer_size::Int=1000)

Returns an iterator over `(key, graph)` pairs from a Graph6 file.

This avoids loading all graphs into memory at once.

# Arguments
- `path::String`: Path to the Graph6 file
- `buffer_size::Int`: Size of the channel buffer (default: 1000)

# Returns
- An iterator that yields `(key, graph)` pairs
"""
function each_g6_graph(path::String; buffer_size::Int=1000)
    # 预分配一个足够大的 BitVector 以避免重复分配
    temp_bitvec = BitVector(undef, 1024)
    
    # 预分配字符串缓冲区
    key_buffer = String[]
    g6_buffer = String[]
    
    return Channel{Tuple{String, SimpleGraph{Int}}}(buffer_size) do ch
        open(path, "r") do io
            line_num = 1
            for line in eachline(io)
                line = strip(line)
                isempty(line) && continue

                # 使用预分配的缓冲区
                parts = split(line, ' ', limit=2)
                key = length(parts) == 2 ? parts[1] : "graph$(line_num)"
                g6_code = length(parts) == 2 ? strip(parts[2]) : parts[1]

                try
                    graph = _parse_g6_string(g6_code, temp_bitvec)
                    put!(ch, (key, graph))
                catch e
                    @warn "Failed to parse graph at line $line_num: $e"
                end

                line_num += 1
            end
        end
    end
end

"""
    save_g6_graphs(graphs::Dict{String, SimpleGraph{Int}}, path::String)

Save multiple graphs in Graph6 format to a file.

# Arguments
- `graphs::Dict{String, SimpleGraph{Int}}`: Dictionary of graphs to save
- `path::String`: Path to save the graphs to
"""
function save_g6_graphs(graphs::Dict{String, SimpleGraph{Int}}, path::String)
    open(path, "w") do io
        # Sort keys by numeric order
        sorted_keys = sort(collect(keys(graphs)), by=k -> parse(Int, match(r"\d+", k).match))

        # Pre-allocate buffer
        buffer = IOBuffer()

        for key in sorted_keys
            # Clear buffer
            truncate(buffer, 0)
            seek(buffer, 0)

            # Write graph's g6 encoding to buffer
            savegraph(buffer, graphs[key], "", GraphIO.Graph6.Graph6Format())

            # Get g6 code and write to file
            g6_code = String(take!(buffer))
            write(io, key * " " * strip(g6_code) * "\n")
        end
    end
end

function save_g6_graph(g::SimpleGraph{Int}, path::String)
    g_dict = Dict("graph1" => g)
    save_g6_graphs(g_dict, path)
end



