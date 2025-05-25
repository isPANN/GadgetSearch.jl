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

    # Pre-allocate BitVector to reduce memory allocations
    temp_bitvec = BitVector(undef, 1024)  # Initial size, will be adjusted as needed

    open(path, "r") do io
        line_num = 1
        for line in eachline(io)
            # Skip empty lines
            isempty(strip(line)) && continue

            # Split key name and g6 code
            parts = split(line, ' ', limit=2)

            if length(parts) == 2
                key = parts[1]
                g6_code = strip(parts[2])

                # Use optimized parsing function
                graph = _parse_g6_string(g6_code, temp_bitvec)
                graph_dict[key] = graph
            else
                # If there's no space in the line, assume the entire line is a g6 code
                graph = _parse_g6_string(parts[1], temp_bitvec)
                graph_dict["graph$(line_num)"] = graph
                line_num += 1
            end
        end
    end

    return graph_dict
end

function read_g6_graph(path::String, id::Int)
    graph_dict = read_g6_graphs(path)
    return graph_dict["graph$(id)"]
end

function _parse_g6_string(s::AbstractString, temp_bitvec::BitVector)
    if startswith(s, ">>graph6<<")
        s = s[11:end]
    end

    # Process characters directly instead of converting to Vector{UInt8}
    chars = collect(s)

    # Parse number of vertices
    if chars[1] <= '~'
        nv = Int(chars[1]) - 63
        pos = 2
    elseif length(chars) > 1 && chars[2] <= '~'
        # Handle larger vertex counts
        nv = 0
        for i in 2:4
            i > length(chars) && break
            nv = (nv << 6) + (Int(chars[i]) - 63)
        end
        pos = 5
    else
        # Handle even larger vertex counts
        nv = 0
        for i in 3:8
            i > length(chars) && break
            nv = (nv << 6) + (Int(chars[i]) - 63)
        end
        pos = 9
    end

    # Create graph
    g = Graphs.SimpleGraph(nv)

    # Calculate required number of bits
    nbits = div(nv * (nv - 1), 2)

    # Ensure temp_bitvec is large enough
    if length(temp_bitvec) < nbits
        resize!(temp_bitvec, nbits)
    end

    # Parse edge information
    bit_idx = 1
    while pos <= length(chars) && bit_idx <= nbits
        c = Int(chars[pos]) - 63
        for i in 5:-1:0
            if bit_idx <= nbits
                temp_bitvec[bit_idx] = ((c >> i) & 1) == 1
                bit_idx += 1
            end
        end
        pos += 1
    end

    # Add edges based on bit vector
    bit_idx = 1
    for col in 2:nv, row in 1:(col-1)
        if bit_idx <= nbits && temp_bitvec[bit_idx]
            add_edge!(g, row, col)
        end
        bit_idx += 1
    end

    return g
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

# """
#     read_g6_graph(path::String, id::Int)

# Read a single graph in Graph6 format from a file.

# # Arguments
# - `path::String`: Path to the file containing Graph6 encoded graphs
# - `id::Int`: ID of the graph to read

# # Returns
# - A SimpleGraph object representing the graph
# """
# function read_g6_graph(path::String, id::Int)
#     graph = Graphs.loadgraph(path, "graph$(id)", GraphIO.Graph6.Graph6Format())
#     return graph
# end

# """
#     save_g6_graph(graph::SimpleGraph{Int}, path::String)

# Save a single graph in Graph6 format to a file.

# # Arguments
# - `graph::SimpleGraph{Int}`: The graph to save
# - `path::String`: Path to save the graph to
# """
# function save_g6_graph(graph::SimpleGraph{Int}, path::String)
#     open(joinpath(path), "w") do io
#         savegraph(io, graph, "", GraphIO.Graph6.Graph6Format())
#     end
# end



# """
#     save_graph_batch(graph_dict::Dict{String, SimpleGraph{Int}},
#                      features::Dict{String, Dict{String, Any}},
#                      base_path::String,
#                      batch_id::Int)

# Save a batch of graphs and their features to files.

# # Arguments
# - `graph_dict::Dict{String, SimpleGraph{Int}}`: Dictionary of graphs to save
# - `features::Dict{String, Dict{String, Any}}`: Dictionary of features for each graph
# - `base_path::String`: Base directory to save files to
# - `batch_id::Int`: Batch identifier

# # Returns
# - Tuple of paths to the saved graph file and features file
# """
# function save_graph_batch(graph_dict::Dict{String, SimpleGraph{Int}},
#                          features::Dict{String, Dict{String, Any}},
#                          base_path::String,
#                          batch_id::Int)
#     # Create directory if it doesn't exist
#     isdir(base_path) || mkpath(base_path)

#     # Build file paths
#     g6_path = joinpath(base_path, "batch_$(batch_id).g6")
#     features_path = joinpath(base_path, "batch_$(batch_id)_features.jld2")

#     # Save graphs to g6 file
#     save_g6_graphs(graph_dict, g6_path)

#     jldsave(features_path; features=features)

#     return g6_path, features_path
# end

# """
#     read_graph_batch(base_path::String, batch_id::Int)

# 读取一批图及其特征。

# # 参数
# - `base_path::String`: 基础保存路径
# - `batch_id::Int`: 批次ID

# # 返回
# - 图字典和特征字典的元组
# """
# function read_graph_batch(base_path::String, batch_id::Int)
#     # 构建文件路径
#     g6_path = joinpath(base_path, "batch_$(batch_id).g6")
#     features_path = joinpath(base_path, "batch_$(batch_id)_features.jld2")

#     # 检查文件是否存在
#     if !isfile(g6_path) || !isfile(features_path)
#         @warn "Batch $(batch_id) files not found at $(base_path)"
#         return Dict{String, SimpleGraph{Int}}(), Dict{String, Dict{String, Any}}()
#     end

#     # Read graphs
#     graph_dict = read_g6_graphs(g6_path)

#     # Read features
#     local features_data
#     jldopen(features_path, "r") do file
#         features_data = file["features"]
#     end

#     return graph_dict, features_data
# end

# """
#     GraphBatchManager(base_path::String; max_graphs_per_batch::Int=10000)

# 管理多个图批次的生成、存储和检索。

# # 字段
# - `base_path::String`: 批次文件的基础路径
# - `max_graphs_per_batch::Int`: 每批的最大图数量
# - `current_batch_id::Int`: 当前批次ID
# - `current_graph_dict::Dict{String, SimpleGraph{Int}}`: 当前批次的图字典
# - `current_features::Dict{String, Dict{String, Any}}`: 当前批次的特征字典
# - `batch_info::Vector{Dict{String, Any}}`: 所有批次的信息
# """
# mutable struct GraphBatchManager
#     base_path::String
#     max_graphs_per_batch::Int
#     current_batch_id::Int
#     current_graph_dict::Dict{String, SimpleGraph{Int}}
#     current_features::Dict{String, Dict{String, Any}}
#     batch_info::Vector{Dict{String, Any}}

#     function GraphBatchManager(base_path::String; max_graphs_per_batch::Int=10000)
#         # 创建目录（如果不存在）
#         isdir(base_path) || mkpath(base_path)

#         # 检查是否有现有批次
#         info_path = joinpath(base_path, "batch_info.json")
#         batch_info = Vector{Dict{String, Any}}()

#         if isfile(info_path)
#             try
#                 batch_info = JSON.parsefile(info_path)
#                 @info "Found $(length(batch_info)) existing batches"
#             catch e
#                 @warn "Error reading batch info: $e"
#             end
#         end

#         # 确定下一个批次ID
#         next_batch_id = isempty(batch_info) ? 1 : maximum(b["batch_id"] for b in batch_info) + 1

#         new(
#             base_path,
#             max_graphs_per_batch,
#             next_batch_id,
#             Dict{String, SimpleGraph{Int}}(),
#             Dict{String, Dict{String, Any}}(),
#             batch_info
#         )
#     end
# end

# """
#     add_graph!(manager::GraphBatchManager, graph_key::String, graph::SimpleGraph{Int}, features::Dict{String, Any})

# 向批次管理器添加一个图及其特征。如果当前批次已满，会自动保存并开始新批次。

# # 参数
# - `manager::GraphBatchManager`: 批次管理器
# - `graph_key::String`: 图的键
# - `graph::SimpleGraph{Int}`: 图对象
# - `features::Dict{String, Any}`: 图的特征

# # 返回
# - 如果添加导致批次保存，返回保存的批次ID；否则返回0
# """
# function add_graph!(manager::GraphBatchManager, graph_key::String, graph::SimpleGraph{Int}, features::Dict{String, Any})
#     # 添加图和特征到当前批次
#     manager.current_graph_dict[graph_key] = graph
#     manager.current_features[graph_key] = features

#     # 检查是否需要保存当前批次
#     if length(manager.current_graph_dict) >= manager.max_graphs_per_batch
#         saved_batch_id = save_current_batch!(manager)
#         return saved_batch_id
#     end

#     return 0  # 没有保存批次
# end

# """
#     save_current_batch!(manager::GraphBatchManager)

# 保存当前批次并开始新批次。

# # 参数
# - `manager::GraphBatchManager`: 批次管理器

# # 返回
# - 保存的批次ID
# """
# function save_current_batch!(manager::GraphBatchManager)
#     # 如果当前批次为空，不保存
#     if isempty(manager.current_graph_dict)
#         return 0
#     end

#     # 保存当前批次
#     g6_path, features_path = save_graph_batch(
#         manager.current_graph_dict,
#         manager.current_features,
#         manager.base_path,
#         manager.current_batch_id
#     )

#     # 记录批次信息
#     batch_info = Dict{String, Any}(
#         "batch_id" => manager.current_batch_id,
#         "g6_path" => g6_path,
#         "features_path" => features_path,
#         "graph_count" => length(manager.current_graph_dict),
#         "timestamp" => string(Dates.now())
#     )
#     push!(manager.batch_info, batch_info)

#     # 保存批次信息
#     info_path = joinpath(manager.base_path, "batch_info.json")
#     open(info_path, "w") do io
#         JSON.print(io, manager.batch_info)
#     end

#     # 记录已保存的批次ID
#     saved_batch_id = manager.current_batch_id

#     # 开始新批次
#     manager.current_batch_id += 1
#     manager.current_graph_dict = Dict{String, SimpleGraph{Int}}()
#     manager.current_features = Dict{String, Dict{String, Any}}()

#     @info "Saved batch $(saved_batch_id) with $(batch_info["graph_count"]) graphs"
#     return saved_batch_id
# end

# """
#     get_batch_count(manager::GraphBatchManager)

# 获取批次管理器中的批次数量。

# # 参数
# - `manager::GraphBatchManager`: 批次管理器

# # 返回
# - 批次数量
# """
# function get_batch_count(manager::GraphBatchManager)
#     return length(manager.batch_info)
# end

# """
#     get_total_graph_count(manager::GraphBatchManager)

# 获取批次管理器中的总图数量。

# # 参数
# - `manager::GraphBatchManager`: 批次管理器

# # 返回
# - 总图数量
# """
# function get_total_graph_count(manager::GraphBatchManager)
#     total = sum(b["graph_count"] for b in manager.batch_info; init=0)
#     total += length(manager.current_graph_dict)  # 加上当前未保存的批次
#     return total
# end

# """
#     process_all_batches(manager::GraphBatchManager, process_function;
#                        start_batch::Int=1, end_batch::Int=0)

# Process all batches in the batch manager.

# # Arguments
# - `manager::GraphBatchManager`: The batch manager
# - `process_function`: Processing function that accepts (graph_dict, features, batch_id) parameters
# - `start_batch::Int=1`: Starting batch ID
# - `end_batch::Int=0`: Ending batch ID, 0 means process all batches

# # Returns
# - Collection of all batch processing results
# """
# function process_all_batches(manager::GraphBatchManager, process_function;
#                            start_batch::Int=1, end_batch::Int=0)
#     # Determine batches to process
#     batch_ids = [b["batch_id"] for b in manager.batch_info]
#     if end_batch > 0
#         batch_ids = filter(id -> start_batch <= id <= end_batch, batch_ids)
#     else
#         batch_ids = filter(id -> id >= start_batch, batch_ids)
#     end

#     # If current batch has data and is in processing range, save it
#     if !isempty(manager.current_graph_dict) &&
#        (end_batch == 0 || manager.current_batch_id <= end_batch)
#         saved_id = save_current_batch!(manager)
#         if saved_id > 0 && saved_id >= start_batch
#             push!(batch_ids, saved_id)
#         end
#     end

#     # Process batches
#     results = []

#     # Sequential processing
#     for id in batch_ids
#         graph_dict, features = read_graph_batch(manager.base_path, id)
#         result = process_function(graph_dict, features, id)
#         !isnothing(result) && push!(results, result)
#     end

#     return results
# end