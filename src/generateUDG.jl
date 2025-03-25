function filter_nonisomorphic(gs::AbstractVector{<:SimpleGraph}, fixed_vertices=nothing)
    # use characteristic_number to classify graphs
    groups = Dict{UInt64, Vector{Tuple{SimpleGraph, Int}}}()
    for (idx, g) in enumerate(gs)
        cn = characteristic_number(g)
        push!(get!(groups, cn, []), (g, idx))
    end

    # filter nonisomorphic graphs
    unique_graphs = SimpleGraph[]
    unique_indices = Int[]
    for bucket in values(groups)
        graphs = [g for (g, _) in bucket]
        indices = [idx for (_, idx) in bucket]
        filtered_graphs = filter_nonisomorphic_naive(graphs, fixed_vertices)
        append!(unique_graphs, filtered_graphs)
        append!(unique_indices, [indices[findfirst(==(g), graphs)] for g in filtered_graphs])
    end

    return unique_graphs, unique_indices
end

function characteristic_number(g::SimpleGraph)
    # up to second order
    hash(sort([hash(sort(degree.(Ref(g), neighbors(g, i)))) for i = 1:nv(g)]))   # degrees
end

function filter_nonisomorphic_naive(gs::AbstractVector{<:SimpleGraph}, fixed_vertices)
    ng = length(gs)
    mask = trues(ng)
    for i = 1:ng
        mask[i] || continue # skip if already filtered
        for j = i+1:ng
            mask[j] || continue
            if fixed_vertices === nothing
                if Graphs.Experimental.has_isomorph(gs[i], gs[j])
                    mask[j] = false
                end
            else
                weights = [i <= sum(fixed_vertices) ? 1 : 0 for i in 1:nv(gs[i])]
                if Graphs.Experimental.has_isomorph(gs[i], gs[j]; vertex_relation=(u, v) -> weights[u] == weights[v])
                    mask[j] = false
                end
            end
        end
    end
    return gs[mask]
end

function canonical_form(comb, m, n)
    # generate all possible different combinations
    sorted_comb = sort(comb)

    # rotate 90, 180, 270
    r180 = sort([(m - x + 1, n - y + 1) for (x, y) in sorted_comb])
    
    # reflect x, y
    ref_x = sort([(m - x + 1, y) for (x, y) in sorted_comb])
    ref_y = sort([(x, n - y + 1) for (x, y) in sorted_comb])
    
    if n == m
        r90 = sort([(y, m - x + 1) for (x, y) in sorted_comb])
        r270 = sort([(n - y + 1, x) for (x, y) in sorted_comb])
        # reflect diagonal
        ref_diag = sort([(y, x) for (x, y) in sorted_comb])
        # return the hash of the minimum of all combinations
        return hash(minimum([sorted_comb, r90, r180, r270, ref_x, ref_y, ref_diag]))
    else
        return hash(minimum([sorted_comb, r180, ref_x, ref_y]))
    end
end

# returns a vector of grid configurations, classified by the number of nodes
function generate_grid_udgs(m::Int, n::Int, pin_pad::Int, directions::Vector{Symbol}=[:up, :right, :down, :left]; min_num_node=min(m, n), max_num_node=m * n, radius=1.6, save_path="")
    # m - the number of rows; n - the number of columns
    unique_graphs = SimpleGraph{Int}[]
    pos_list = Vector{Int}[]

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

    for indices in product((1:length(p) for p in pivot_pos)...)
        pivot_combination = [pivot_pos[i][index] for (i, index) in enumerate(indices)]
        pin_combination = [pin_pos[i][index] for (i, index) in enumerate(indices)]
        candidate_points = [
            (i, j) 
            for i in (1+pin_pad):(m+pin_pad), j in (1+pin_pad):(n+pin_pad)
            if !((i, j) in pivot_combination)
        ]
        seen = Dict{UInt64, Int}()
        # Sequentially increase the number of selected points in the lattice.
        for num_node = min_num_node - length(pivot_combination):max_num_node - length(pivot_combination)
            for comb in Combinatorics.combinations(candidate_points, num_node)
                full_comb = vcat(comb, pivot_combination...)
                can_form = canonical_form(full_comb, m+2pin_pad, n+2pin_pad)
                
                if haskey(seen, can_form)
                    continue 
                end
                seen[can_form] = length(pos_list) + 1

                full_pos = [pin_combination..., full_comb...]
                udg = unit_disk_graph(full_pos, radius)
                if is_connected(udg)
                    push!(unique_graphs, udg)
                    push!(pos_list, _tuple_to_index(full_pos, (m+2pin_pad, n+2pin_pad)))
                end
            end
        end
    end
    final_unique_graphs, idx = filter_nonisomorphic(unique_graphs, 1:length(directions))
    # for i in length(final_unique_graphs)
    #     @assert nv(final_unique_graphs[i]) == length(pos_list[idx[i]])
    # end
    @info "After embedding unique graphs: $(length(final_unique_graphs))"
    if length(save_path) > 0
        filename = "m$(m)n$(n)pad$(pin_pad)_min$(min_num_node)max$(max_num_node)_direct$(length(directions))"
        filename_g = joinpath(save_path, filename * ".g6")
        filename_pos = joinpath(save_path, filename * ".json")
        save_grid_udgs(filename_g, final_unique_graphs, filename_pos, pos_list[idx])
        return filename_g, filename_pos
    else
        return final_unique_graphs, pos_list[idx]
    end
end

function save_grid_udgs(filename_g::String, graphs::Vector{SimpleGraph}, filename_pos::String, pos_list::AbstractVector)
    open(joinpath(filename_g), "w") do io
        for (i, g) in enumerate(graphs)
            # write(io, "graph$i\n")  # key
            savegraph(io, g, "graph$i", GraphIO.Graph6.Graph6Format())  # graph
        end
    end
    pos_dict = Dict(i => pos_list[i] for i in 1:length(graphs))
    open(joinpath(filename_pos), "w") do io
        write(io, JSON.json(pos_dict))
    end
end