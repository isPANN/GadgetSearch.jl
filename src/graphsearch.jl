function search_single_constraint(vertex_nums::Vector{Int}, bit_num::Int, degeneracy::Vector{Int}, dir_path::String, save_path::String)
    all_graph_data_for_degeneracy = []
    for vertex_num in vertex_nums
        graph_path = joinpath(dir_path, "graph$(vertex_num).g6")
        graph_dict = readgraphdict(graph_path)
        gname, candidate, weight = check_single_constraint(graph_dict, bit_num, degeneracy)
        isnothing(gname) && continue

        g = graph_dict[gname]
        nodes = [Dict("id" => v, "weight" => weight[v]) for v in Graphs.vertices(g)]
        edges = [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(g)]
    
        push!(all_graph_data_for_degeneracy, Dict(
            "degeneracy" => [join(bin(elem, bit_num), "") for elem in degeneracy],
            # "vertex_num" => vertex_num,
            "nodes" => nodes,
            "edges" => edges,
            "work_nodes" => candidate
        ))
    end
    filename = joinpath(save_path, "$(bit_num)bits_$(degeneracy).json")
    open(filename, "w") do file
        write(file, JSON.json(all_graph_data_for_degeneracy))
    end
    return filename
end

function search_single_constraint(vertex_nums::Vector{Int}, bit_num::Int, degeneracy::Vector{String}, dir_path::String, save_path::String)
    degeneracy_int = [decimal(degen) for degen in degeneracy]
    return search_single_constraint(vertex_nums, bit_num, degeneracy_int, dir_path, save_path)
end

function search_single_gate(vertex_nums::Vector{Int}, input_bits::Int, output_bits::Int, gate_id::Int, dir_path::String, save_path::String)
    degeneracy = genericgate(gate_id, input_bits, output_bits)
    all_graph_data_for_gate = []
    for vertex_num in vertex_nums
        graph_path = joinpath(dir_path, "graph$(vertex_num).g6")
        graph_dict = readgraphdict(graph_path)
        gname, candidate, weight = check_single_constraint(graph_dict, input_bits + output_bits, degeneracy)
        isnothing(gname) && continue
        
        g = graph_dict[gname]
        nodes = [Dict("id" => v, "weight" => weight[v]) for v in Graphs.vertices(g)]
        edges = [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(g)]
        push!(all_graph_data_for_gate, Dict(
            "gate_id" => gate_id,
            "degeneracy" => [join(bin(elem, input_bits + output_bits), "") for elem in degeneracy],
            # "vertex_num" => vertex_num,
            "nodes" => nodes,
            "edges" => edges,
            "work_nodes" => candidate
        ))
    end
    filename = joinpath(save_path, "$(input_bits)in_$(output_bits)out_gate$(gate_id)_$(degeneracy).json")
    open(filename, "w") do file
        write(file, JSON.json(all_graph_data_for_gate))
    end
    return filename
end

function search_gates(vertex_nums::Vector{Int}, input_bits::Int, output_bits::Int, dir_path::String, save_path::String)
    gate_num = (2^output_bits)^(2^input_bits)
    all_graph_data = []
    # gate_list = 
    filename = joinpath(save_path, "$(input_bits)in$(output_bits)out.json")
    open(filename, "w") do file
        write(file, "[")
        first_entry = true
        for gate_id in gate_list
            degeneracy = genericgate(gate_id, input_bits, output_bits)
    
            for vertex_num in vertex_nums
                graph_path = joinpath(dir_path, "graph$(vertex_num).g6")
                # graph_dict = readgraphdict(graph_path)
                g_matrix, candidate, weight = check_single_gate_traversal(graph_path, input_bits, output_bits, degeneracy)
                if isnothing(g_matrix)
                    # @info "No suitable graph found for gate $(gate_id) in Graph $(vertex_num)."
                    continue
                end
                # g = graph_dict[gname]
                g = adjacency_matrix_to_simple_graph(g_matrix)
                nodes = [Dict("id" => v, "weight" => weight[v]) for v in Graphs.vertices(g)]
                edges = [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(g)]
                # push!(all_graph_data, Dict(
                #     "gate_id" => gate_id,
                #     "degeneracy" => [join(bin(elem, input_bits + output_bits), "") for elem in degeneracy],
                #     # "vertex_num" => vertex_num,
                #     "nodes" => nodes,
                #     "edges" => edges,
                #     "work_nodes" => candidate
                # ))
                if !first_entry
                    write(file, ",")
                end
                first_entry = false
                graph_entry = Dict(
                    "gate_id" => gate_id,
                    "degeneracy" => [join(bin(elem, input_bits + output_bits), "") for elem in degeneracy],
                    # "vertex_num" => vertex_num,
                    "nodes" => nodes,
                    "edges" => edges,
                    "work_nodes" => candidate
                )
                write(file, JSON.json(graph_entry))
                # found = true
                break 
            end
        end
        write(file, "]")
        # write(file, JSON.json(all_graph_data))

    end
    return filename
end

function search_any_constraint(vertex_nums::Vector{Int}, bit_num::Int, dir_path::String, save_path::String)
    all_graph_data = []
    for (idx, degeneracy) in enumerate([collect(s) for s in IterTools.subsets(0:2^bit_num-1) if !isempty(s)])
        found = false
        for vertex_num in vertex_nums
            graph_path = joinpath(dir_path, "graph$(vertex_num).g6")
            graph_dict = readgraphdict(graph_path)
            gname, candidate, weight = check_single_constraint(graph_dict, bit_num, degeneracy)
            isnothing(gname) && continue

            g = graph_dict[gname]
            nodes = [Dict("id" => v, "weight" => weight[v]) for v in Graphs.vertices(g)]
            edges = [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(g)]
            push!(all_graph_data, Dict(
                "degeneracy" => [join(bin(elem, bit_num), "") for elem in degeneracy],
                "gate_id" => idx-1,
                "nodes" => nodes,
                "edges" => edges,
                "work_nodes" => candidate
            ))
            found = true
            break 
        end
        !found && @info "No suitable graph found for degeneracy $(degeneracy) in Graph $(vertex_nums)."
    end
    filename = joinpath(save_path, "$(bit_num)bits_any_constraint.json")
    open(filename, "w") do file
        write(file, JSON.json(all_graph_data))
    end
    return filename
end

function check_single_constraint(graphs::Dict{String, Graphs.SimpleGraphs.SimpleGraph}, bit_num::Int, degeneracy::Vector{Int})
    @assert length(degeneracy) > 0 && maximum(degeneracy) < 2^bit_num
    for gname in sort(collect(keys(graphs)))
        # Check if the graph is connected.
        is_connected(graphs[gname]) || continue
        vertex_num = nv(graphs[gname])
        # Find Maximal Independent Sets for each graph using GenericTensorNetworks.jl(https://queracomputing.github.io/GenericTensorNetworks.jl/dev/generated/MaximalIS/)
        # maximalis = MaximalIS(graphs[gname])
        # mis_problem = GenericTensorNetwork(maximalis)
        # mis_result = read_config(solve(mis_problem, ConfigsAll())[])
        ones_vertex = maximal_cliques(complement(graphs[gname]))
        mis_result = generate_bitvectors(vertex_num, ones_vertex)
        mis_num = length(mis_result)
        
        all_candidates = collect(permutations(1:vertex_num, bit_num))
        
        for candidate in all_candidates 
            candidate_value = [[parse(Int, mis_result[i][j]) for j in candidate] for i in 1:mis_num]
            candidate_value_decimal = [decimal(candidate_value[i]) for i in 1:mis_num]

            is_subset = all(x -> x in candidate_value_decimal, degeneracy)
            is_subset || continue

            target_mis_indices = [findfirst(==(x), candidate_value_decimal) for x in degeneracy]
            wrong_mis_indices = setdiff(1:mis_num, target_mis_indices)
            target_mis_sets = mis_result[target_mis_indices]
            wrong_mis_sets = mis_result[wrong_mis_indices]

            model = Model(HiGHS.Optimizer)
            set_silent(model)
            # The number of variables = the number of vertices in the graph.
            @variable(model, x[1:vertex_num])
            all_mis_combinations = collect(combinations(target_mis_sets, 2))
            for mis in all_mis_combinations
                # @info "Adding constraint for $(mis)."
                @constraint(model, sum((parse(Int, mis[1][i]) - parse(Int, mis[2][i])) * x[i] for i in 1:vertex_num) == 0)
            end 

            ϵ = -1
            for wrong_mis in wrong_mis_sets
                for mis in target_mis_sets
                    # @info "Adding constraint for $(wrong_mis)_$(mis)."
                    @constraint(model, sum((parse(Int, mis[i]) - parse(Int, wrong_mis[i])) * x[i] for i in 1:vertex_num) <= ϵ)
                end
            end

            for i in 1:vertex_num
                @constraint(model, 1 <= x[i])
            end
            
            @objective(model, Min, sum(x[i] for i in 1:vertex_num))
            optimize!(model)
            if is_solved_and_feasible(model)
                @info "Optimization successful for $(vertex_num)_$(gname)_$(candidate)."
                @info "Optimal objective value: $([value(x[i]) for i in 1:length(x)])"
                return gname, candidate, [value(x[i]) for i in 1:length(x)]
            end
            for con in all_constraints(model; include_variable_in_set_constraints = false)
                delete(model, con)
            end
        end
    end
    return nothing, nothing, nothing
end

function generate_candidates(g::Graphs.SimpleGraphs.SimpleGraph, input_num::Int, mis::Vector{Vector{Int}})
    sets = []
    transposed = reduce(hcat, mis)'
    col_sums = sum(transposed, dims=1)
    valid_columns = findall(x -> x >= 2^(input_num - 1), col_sums[:])
    if length(valid_columns) < input_num
        return sets
    end
    for subset in combinations(valid_columns, input_num)
        is_independent = true
        for (u, v) in combinations(subset, 2)
            if has_edge(g, u, v)
                is_independent = false
                break
            end
        end
        if is_independent
            push!(sets, subset)
        end
    end
    return sets
end

function generate_candidates(adj_matrix::Matrix{Int}, input_num::Int, mis::Vector{Vector{Int}})
    transposed = reduce(hcat, mis)'
    col_sums = sum(transposed, dims=1)
    valid_columns = findall(x -> x >= 2^(input_num - 1), col_sums[:])
    if length(valid_columns) < input_num
        return []
    end
    sets = Vector{Vector{Int}}()
    explicit_combinations = collect(combinations(valid_columns, input_num))

    for subset in explicit_combinations
        is_independent = true
        for (u, v) in combinations(subset, 2)
            # 检查是否有边，若有边，则不是独立集
            if adj_matrix[u, v] == 1
                is_independent = false
                break
            end
        end
        # 如果是独立集，加入结果
        if is_independent
            push!(sets, collect(subset))  # 将组合添加到结果集
        end
    end
    return sets
end

function check_single_gate_traversal_form(file_name::String, input_num::Int, output_num::Int, degeneracy::Vector{Int})
    open(file_name, "r") do io
        for line in eachline(io)
            line_io = IOBuffer(line)
            # g = loadgraph(line_io, "graph1", GraphIO.Graph6.Graph6Format())
            g = loadgraphs(line_io, GraphIO.Graph6.Graph6Format())

            gname, candidate, vertex_weight = check_single_gate(g, input_num, output_num, degeneracy)
            # for gname in collect(keys(g))
            #     ones_vertex = Graphs.maximal_cliques(Graphs.complement(g[gname]))
            #     mis_result = generate_bitvectors(9, ones_vertex)
            # end 
            !isnothing(gname) && return gname, candidate, vertex_weight
        end
    end
end

function find_maximal_cliques(adj_matrix::Matrix{Int})
    n = size(adj_matrix, 1)  # Number of vertices
    maxconn = -1
    pivotnbrs = Set{Int}()
    pivotdonenbrs = Set{Int}()
    
    # 预处理邻居信息
    nnbrs = Vector{Set{Int}}(undef, n)
    for i in 1:n
        nnbrs[i] = Set{Int}()
        for j in 1:n
            if i != j && adj_matrix[i, j] != 0
                push!(nnbrs[i], j)
            end
        end
        # 找到度数最高的点作为初始 pivot
        if length(nnbrs[i]) > maxconn
            maxconn = length(nnbrs[i])
            pivotnbrs = nnbrs[i]
        end
    end

    # 初始化搜索栈
    cand = Set{Int}(1:n)
    smallcand = setdiff(cand, pivotnbrs)
    done = Set{Int}()
    stack = Vector{Tuple{Set{Int}, Set{Int}, Set{Int}}}()
    clique_so_far = Vector{Int}()
    cliques = Vector{Vector{Int}}()

    # 主循环
    while !isempty(smallcand) || !isempty(stack)
        if !isempty(smallcand)
            v = pop!(smallcand)
        else
            # 回溯
            cand, done, smallcand = pop!(stack)
            pop!(clique_so_far)
            continue
        end

        # 将节点加入当前团
        push!(clique_so_far, v)
        delete!(cand, v)
        push!(done, v)
        new_cand = intersect(cand, nnbrs[v])
        new_done = intersect(done, nnbrs[v])

        # 检查是否形成一个团
        if isempty(new_cand)
            if isempty(new_done)
                push!(cliques, copy(clique_so_far))
            end
            pop!(clique_so_far)
            continue
        end

        # 快速路径：只剩一个候选节点
        if isempty(new_done) && length(new_cand) == 1
            push!(cliques, vcat(clique_so_far, collect(new_cand)))
            pop!(clique_so_far)
            continue
        end

        # 找到 pivot 节点
        numb_cand = length(new_cand)
        maxconndone = -1
        for u in new_done
            cn = intersect(new_cand, nnbrs[u])
            conn = length(cn)
            if conn > maxconndone
                pivotdonenbrs = cn
                maxconndone = conn
                if maxconndone == numb_cand
                    break
                end
            end
        end
        if maxconndone == numb_cand
            pop!(clique_so_far)
            continue
        end

        maxconn = -1
        for u in new_cand
            cn = intersect(new_cand, nnbrs[u])
            conn = length(cn)
            if conn > maxconn
                pivotnbrs = cn
                maxconn = conn
                if maxconn == numb_cand - 1
                    break
                end
            end
        end

        if maxconndone > maxconn
            pivotnbrs = pivotdonenbrs
        end

        # 保存搜索状态
        push!(stack, (cand, done, smallcand))
        cand = new_cand
        done = new_done
        smallcand = setdiff(cand, pivotnbrs)
    end

    return cliques
end

function is_connected(graph::Matrix{Int})
    n = size(graph, 1)  # 图的节点数
    visited = falses(n)  # 访问标记数组
    
    function dfs(node::Int)
        visited[node] = true
        for neighbor in 1:n
            if graph[node, neighbor] > 0 && !visited[neighbor]
                dfs(neighbor)
            end
        end
    end
    dfs(1)
    return all(visited)
end

function check_single_gate_traversal(file_name::String, input_num::Int, output_num::Int, degeneracy::Vector{Int})
    result = nothing  # 预先定义返回值
    open(file_name, "r") do io
        for line in eachline(io)
            g = g6string_to_matrix(line)
            candidate, vertex_weight = check_single_gate(g, input_num, output_num, degeneracy)
            if !isnothing(candidate)
                result = (g, candidate, vertex_weight)
                break  # 找到结果后立即跳出循环
            end
        end
    end
    return result === nothing ? (nothing, nothing, nothing) : result
end

function check_single_gate(graph::Matrix{Int}, input_num::Int, output_num::Int, degeneracy::Vector{Int})
    bit_num = input_num + output_num
    @assert length(degeneracy) > 0 && maximum(degeneracy) < 2^bit_num
    is_connected(graph) || return nothing, nothing
    vertex_num = size(graph, 1)
    ones_vertex = find_maximal_cliques(complement(graph))
    # @show ones_vertex
    mis_result = generate_bitvectors(vertex_num, ones_vertex)
    mis_num = length(mis_result)
    if mis_num < 2^(input_num)
        return nothing, nothing
    end
    all_candidates = generate_candidates(graph, input_num, mis_result)
    # @show all_candidates
    if isempty(all_candidates)
        return nothing, nothing
    end

    for candidate in all_candidates     
        remain_elements = setdiff(1:vertex_num, candidate)
        for output_bits in permutations(remain_elements, output_num)
            candidate_new = vcat(candidate, output_bits)
            # candidate_value = [[mis_result[i][j] for j in candidate_new] for i in 1:mis_num]
            
            # candidate_value_decimal = [decimal(candidate_value[i]) for i in 1:mis_num]
            _ , candidate_value_decimal = get_values(ones_vertex, candidate_new, vertex_num)
            is_subset = all(x -> x in candidate_value_decimal, degeneracy)
            
            is_subset || continue

            target_mis_indices = [findall(==(x), candidate_value_decimal) for x in degeneracy]
    
            for combination in Iterators.product(target_mis_indices...)
                flattened_target_mis_indices = vcat(combination...)
                # println("Current combination: ", flattened_target_mis_indices)
                wrong_mis_indices = setdiff(1:mis_num, flattened_target_mis_indices)
                target_mis_sets = mis_result[flattened_target_mis_indices]
                wrong_mis_sets = mis_result[wrong_mis_indices]

                model = Model(HiGHS.Optimizer)
                set_silent(model)
                # The number of variables = the number of vertices in the graph.
                @variable(model, x[1:vertex_num])
                all_mis_combinations = collect(combinations(target_mis_sets, 2))
                for mis in all_mis_combinations
                    # @info "Adding constraint for $(mis)."
                    @constraint(model, sum((mis[1][i] - mis[2][i]) * x[i] for i in 1:vertex_num) == 0)
                end 

                ϵ = 1
                for wrong_mis in wrong_mis_sets
                    for mis in target_mis_sets
                        # @info "Adding constraint for $(wrong_mis)_$(mis)."
                        @constraint(model, sum((mis[i] - wrong_mis[i]) * x[i] for i in 1:vertex_num) >= ϵ)
                    end
                end

                for i in 1:vertex_num
                    @constraint(model, 1 <= x[i])
                end
                
                @objective(model, Min, sum(x[i] for i in 1:vertex_num))
                optimize!(model)
                if is_solved_and_feasible(model)
                    # @info "Optimization successful for $(vertex_num)_$(candidate_new)."
                    @info "Optimization successful!"
                    @info "Optimal objective value: $([value(x[i]) for i in 1:length(x)])"
                    return candidate_new, [value(x[i]) for i in 1:length(x)]
                end            
                for con in all_constraints(model; include_variable_in_set_constraints = false)
                    delete(model, con)
                end
            end
        end
    end
    return nothing, nothing
end

function check_single_gate(graphs::Dict{String, Graphs.SimpleGraphs.SimpleGraph}, input_num::Int, output_num::Int, degeneracy::Vector{Int})
    bit_num = input_num + output_num
    @assert length(degeneracy) > 0 && maximum(degeneracy) < 2^bit_num
    for gname in collect(keys(graphs))
        # Check if the graph is connected.
        Graphs.is_connected(graphs[gname]) || continue
        
        vertex_num = nv(graphs[gname])
        ones_vertex = Graphs.maximal_cliques(Graphs.complement(graphs[gname]))
        mis_result = generate_bitvectors(vertex_num, ones_vertex)
        mis_num = length(mis_result)
        if mis_num < 2^(input_num)
            continue
        end

        all_candidates = generate_candidates(graphs[gname], input_num, mis_result)
        if isempty(all_candidates)
            continue
        end
        for candidate in all_candidates
            
            remain_elements = setdiff(1:vertex_num, candidate)

            for output_bits in permutations(remain_elements, output_num)
                candidate_new = vcat(candidate, output_bits)
                candidate_value = [[mis_result[i][j] for j in candidate_new] for i in 1:mis_num]
                
                candidate_value_decimal = [decimal(candidate_value[i]) for i in 1:mis_num]
                is_subset = all(x -> x in candidate_value_decimal, degeneracy)
                is_subset || continue
                # @show candidate_value_decimal
                
                target_mis_indices = [findall(==(x), candidate_value_decimal) for x in degeneracy]
                for combination in Iterators.product(target_mis_indices...)
                    flattened_target_mis_indices = vcat(combination...)
                    # println("Current combination: ", flattened_target_mis_indices)
                    wrong_mis_indices = setdiff(1:mis_num, flattened_target_mis_indices)
                    target_mis_sets = mis_result[flattened_target_mis_indices]
                    wrong_mis_sets = mis_result[wrong_mis_indices]

                    model = Model(HiGHS.Optimizer)
                    set_silent(model)
                    # The number of variables = the number of vertices in the graph.
                    @variable(model, x[1:vertex_num])
                    all_mis_combinations = collect(combinations(target_mis_sets, 2))
                    for mis in all_mis_combinations
                        # @info "Adding constraint for $(mis)."
                        @constraint(model, sum((mis[1][i] - mis[2][i]) * x[i] for i in 1:vertex_num) == 0)
                    end 

                    ϵ = 1
                    for wrong_mis in wrong_mis_sets
                        for mis in target_mis_sets
                            # @info "Adding constraint for $(wrong_mis)_$(mis)."
                            @constraint(model, sum((mis[i] - wrong_mis[i]) * x[i] for i in 1:vertex_num) >= ϵ)
                        end
                    end

                    for i in 1:vertex_num
                        @constraint(model, 1 <= x[i])
                    end
                    
                    @objective(model, Min, sum(x[i] for i in 1:vertex_num))
                    optimize!(model)
                    if is_solved_and_feasible(model)
                        @info "Optimization successful for $(vertex_num)_$(gname)_$(candidate_new)."
                        @info "Optimal objective value: $([value(x[i]) for i in 1:length(x)])"
                        return gname, candidate_new, [value(x[i]) for i in 1:length(x)]
                    end            
                    for con in all_constraints(model; include_variable_in_set_constraints = false)
                        delete(model, con)
                    end
                end
            end
        end
    end
    return nothing, nothing, nothing
end