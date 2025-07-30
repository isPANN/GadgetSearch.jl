function unit_disk_graph(locs::AbstractVector, unit::Real)
    n = length(locs)
    g = SimpleGraph(n)
    for i=1:n, j=i+1:n
        if sum(abs2, locs[i] .- locs[j]) < unit ^ 2
            add_edge!(g, i, j)
        end
    end
    return g
end

function generate_full_ksg(nx::Int, ny::Int; path::String="udg.g6")
    # 获取边界候选
    top_candidates    = [(1, y)     for y in 2:ny+1]
    bottom_candidates = [(nx+2, y)  for y in 2:ny+1]
    left_candidates   = [(x, 1)     for x in 2:nx+1]
    right_candidates  = [(x, ny+2)  for x in 2:nx+1]

    # 中间区域固定选择
    inner_points = [(x, y) for x in 2:nx+1, y in 2:ny+1]
    inner_points = vec(inner_points)

    results = Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}[]
    seen = Set{String}()

    for top in top_candidates, bottom in bottom_candidates,
        left in left_candidates, right in right_candidates

        selected = vcat([top, right, bottom, left], inner_points)
        index_map = Dict{Tuple{Int, Int}, Int}()

        for (i, pt) in enumerate(selected)
            index_map[pt] = i
        end

        g = unit_disk_graph(selected, 1.5)

        coords = selected
        push!(results, (g, coords))
    end
    @info "pinset in generated graphs: [1,2,3,4]"
    return _process_and_save_graphs(results, path)
end

function generate_grid_udgs(m::Int, n::Int, pin_pad::Int, directions::Vector{Symbol}=[:up, :right, :down, :left]; min_num_node=min(m, n), max_num_node=m * n, radius=1.6, save_path="")
    # m - the number of rows; n - the number of columns
    unique_graphs = SimpleGraph{Int}[]
    pos_list = Vector{Tuple{Int, Int}}[]

    pivot_dict = Dict(
        :up => [(1+pin_pad, j) for j in pin_pad+2:pin_pad+n-1],
        :right => [(i, n+pin_pad) for i in pin_pad+2:pin_pad+m-1],
        :down => [(m+pin_pad, j) for j in pin_pad+2:pin_pad+n-1],
        :left => [(i, 1+pin_pad) for i in pin_pad+2:pin_pad+m-1]
    )
    pin_dict = Dict(
        :up => [(1, j) for j in pin_pad+2:pin_pad+n-1],
        :right => [(i, n+2pin_pad) for i in pin_pad+2:pin_pad+m-1],
        :down => [(m+2pin_pad, j) for j in pin_pad+2:pin_pad+n-1],
        :left => [(i, 1) for i in pin_pad+2:pin_pad+m-1]
    )
    pivot_pos = [pivot_dict[d] for d in directions]
    pin_pos = [pin_dict[d] for d in directions]

    for indices in Iterators.product((1:length(p) for p in pivot_pos)...)
        pivot_combination = [pivot_pos[i][index] for (i, index) in enumerate(indices)]
        pin_combination = [pin_pos[i][index] for (i, index) in enumerate(indices)]
        candidate_points = [
            (i, j) 
            for i in (1+pin_pad):(m+pin_pad), j in (1+pin_pad):(n+pin_pad)
            if !((i, j) in pivot_combination)
        ]
        # Removed seen dictionary and canonical_form logic
        # Sequentially increase the number of selected points in the lattice.
        for num_node = min_num_node - length(pivot_combination):max_num_node - length(pivot_combination)
            for comb in Combinatorics.combinations(candidate_points, num_node)
                full_comb = vcat(comb, pivot_combination...)
                full_pos = [pin_combination..., full_comb...]
                udg = unit_disk_graph(full_pos, radius)
                if is_connected(udg)
                    push!(unique_graphs, udg)
                    push!(pos_list, full_pos)
                end
            end
        end
    end
    results = Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}[(unique_graphs[i], pos_list[i]) for i in 1:length(unique_graphs)]
    @info "After embedding unique graphs: $(length(results))"
    if length(save_path) > 0
        return _process_and_save_graphs(results, save_path)
    else
        return results
    end
end

function triangular_lattice_graph_with_pins(m::Int, n::Int, a=1.0; top_pin_col::Int, bottom_pin_col::Int, left_pin_row::Int, right_pin_row::Int)
    h = sqrt(3) * a / 2
    rows = m + 2
    cols = n + 2

    # 计算实际需要的节点数：内部点 + pin点
    inner_points = (m) * (n)  # 内部点数量
    total_nodes = inner_points + 4  # 加上4个pin点

    positions = Vector{Tuple{Float64, Float64}}(undef, total_nodes)

    # 创建坐标映射
    coord_to_idx = Dict{Tuple{Int, Int}, Int}()
    idx = 1

    # 添加pin点
    pin_coords = [
        (1, top_pin_col + 1),        # top
        (left_pin_row + 1, 1),        # left
        (rows, bottom_pin_col + 1),     # bottom
        (right_pin_row + 1, cols),     # right
    ]
    pin_ids = Int[]
    
    for (i, j) in pin_coords
        x_offset = isodd(i) ? 0.5 : 0.0
        x = (j - 1) + x_offset
        y = (i - 1) * h
        positions[idx] = (x, y)
        coord_to_idx[(i, j)] = idx
        push!(pin_ids, idx)
        idx += 1
    end

    # 首先添加内部点
    for i in 2:rows-1
        for j in 2:cols-1
            x_offset = isodd(i) ? 0.5 : 0.0
            x = (j - 1) + x_offset
            y = (i - 1) * h
            positions[idx] = (x, y)
            coord_to_idx[(i, j)] = idx
            idx += 1
        end
    end

    return positions
end

function triangular_lattice_graph(m::Int, n::Int, a=1.0)
    h = sqrt(3) * a / 2
    rows = m
    cols = n

    # 计算实际需要的节点数：内部点 + pin点
    total_nodes = (m) * (n)  # 内部点数量

    positions = Vector{Tuple{Float64, Float64}}(undef, total_nodes)

    # 创建坐标映射
    coord_to_idx = Dict{Tuple{Int, Int}, Int}()
    idx = 1

    # 首先添加内部点
    for i in 1:rows
        for j in 1:cols
            x_offset = isodd(i) ? 0.5 : 0.0
            x = (j - 1) + x_offset
            y = (i - 1) * h
            positions[idx] = (x, y)
            coord_to_idx[(i, j)] = idx
            idx += 1
        end
    end

    return positions
end

function generate_full_triangular_graph(n::Int, m::Int; a=1.0, save_path="triangular_grid_with_pins.g6")
    results = Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}[]
    positions = triangular_lattice_graph(m, n, a)
    g = unit_disk_graph(positions, a+0.1)
    push!(results, (g, positions))
    return _process_and_save_graphs(results, save_path)
end

function generate_full_triangular_graph_with_pins(n::Int, m::Int; a=1.0, save_path="triangular_grid_with_pins.g6")
    results = Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}[]
    for top_pin_col in 1:n, bottom_pin_col in 1:n, left_pin_row in 1:m, right_pin_row in 1:m
        positions = triangular_lattice_graph_with_pins(m, n, a; top_pin_col=top_pin_col, bottom_pin_col=bottom_pin_col, left_pin_row=left_pin_row, right_pin_row=right_pin_row)
        g = unit_disk_graph(positions, a+0.1)
        push!(results, (g, positions))
    end
    return _process_and_save_graphs(results, save_path)
end


# function triangular_grid_with_pins(n::Int, m::Int, a=1.0)
#     coords = Dict{Int, Tuple{Float64, Float64}}()
#     ij_map = Dict{Int, Tuple{Int, Int}}()
#     pin_vertices = Dict(:top => Int[], :bottom => Int[], :left => Int[], :right => Int[], :corner => Int[])
#     h = sqrt(3)/2 * a
#     id = 1

#     for i in 0:n+1
#         for j in 0:m+1
#             x = j * a + (isodd(i) ? a/2 : 0.0)
#             y = i * h
#             coords[id] = (x, y)
#             ij_map[id] = (i, j)

#             # 判断边界类型（优先判断 corner）
#             if (i == 0 || i == n+1) && (j == 0 || j == m+1)
#                 # 四个角：左下、右下、左上、右上
#                 push!(pin_vertices[:corner], id)
#             elseif i == 0
#                 push!(pin_vertices[:bottom], id)
#             elseif i == n+1
#                 push!(pin_vertices[:top], id)
#             elseif j == 0
#                 push!(pin_vertices[:left], id)
#             elseif j == m+1
#                 push!(pin_vertices[:right], id)
#             end

#             id += 1
#         end
#     end

#     return coords, ij_map, pin_vertices
# end

# function generate_full_triangular_graph(n::Int, m::Int; a=1.0, path::String="triangular_grid_with_pins.g6")
#     coords, ij_map, pin_vertices = triangular_grid_with_pins(n, m, a)

#     # Step 1: 所有 pin 组合
#     combos = [(t, b, l, r) for t in pin_vertices[:top],
#                             b in pin_vertices[:bottom],
#                             l in pin_vertices[:left],
#                             r in pin_vertices[:right]]

#     # Step 2: 所有 core 点（不在 pin 里的）
#     all_ids = collect(keys(coords))
#     pin_ids = Set(vcat(pin_vertices[:top], pin_vertices[:bottom],
#                        pin_vertices[:left], pin_vertices[:right], pin_vertices[:corner]))
#     core_ids = filter(id -> !(id in pin_ids), all_ids)

#     # Step 3: 构建所有图：每个由 core + 4 个 pin 构成
#     results = Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}[]

#     for (t, b, l, r) in combos
#         node_ids = vcat([t, l, b, r], core_ids)
#         node_coords = [coords[id] for id in node_ids]
#         g = unit_disk_graph(node_coords, a+0.1)
#         push!(results, (g, node_coords))
#     end
    
#     return _process_and_save_graphs(results, path)
# end


function _process_and_save_graphs(results::Vector{Tuple{SimpleGraph, Vector{Tuple{Float64, Float64}}}}, path::String)
    mapping_file = tempname()
    temp_path = tempname()

    original_coords = getindex.(results, 2)
    save_graph(results, temp_path; g6_only=true)

    run(pipeline(`shortg -v -u $(temp_path)`, stderr=mapping_file))

    canon2orig, _ = _parse_shortg_mapping(mapping_file)
    
    _write_original_representatives(temp_path, canon2orig, original_coords, path)

    rm(mapping_file)
    rm(temp_path)
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