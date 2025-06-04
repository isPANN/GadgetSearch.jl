function run_labelg(input_g6::String, output_g6::String)
    cmd = `labelg -q $input_g6 $output_g6`
    success(run(cmd)) || error("labelg failed")
end

function read_g6_lines(file::String)::Vector{String}
    return filter(!isempty, strip.(readlines(file)))
end

function filter_nonisomorphic(gs::Vector{SimpleGraph}, fixed_vertices=nothing; temp_g6_path="__temp.g6", canon_path="__canon.g6")
    open(temp_g6_path, "w") do io
        for (i, g) in enumerate(gs)
            savegraph(io, g, "graph$i", GraphIO.Graph6.Graph6Format())
        end
    end

    run_labelg(temp_g6_path, canon_path)
    canon_lines = read_g6_lines(canon_path)

    seen = Dict{String, Tuple{SimpleGraph, Int}}()
    for (i, canon) in enumerate(canon_lines)
        seen[canon] = (gs[i], i)
    end

    unique_graphs = [v[1] for v in values(seen)]
    unique_indices = [v[2] for v in values(seen)]

    return unique_graphs, unique_indices
end

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

    for indices in Iterators.product((1:length(p) for p in pivot_pos)...)
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