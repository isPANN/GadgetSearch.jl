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
    for gate_id in 0:(gate_num - 1)
        degeneracy = genericgate(gate_id, input_bits, output_bits)
        found = false
        for vertex_num in vertex_nums
            graph_path = joinpath(dir_path, "graph$(vertex_num).g6")
            graph_dict = readgraphdict(graph_path)
            gname, candidate, weight = check_single_constraint(graph_dict, input_bits + output_bits, degeneracy)
            isnothing(gname) && continue
            
            g = graph_dict[gname]
            nodes = [Dict("id" => v, "weight" => weight[v]) for v in Graphs.vertices(g)]
            edges = [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(g)]
            push!(all_graph_data, Dict(
                "gate_id" => gate_id,
                "degeneracy" => [join(bin(elem, input_bits + output_bits), "") for elem in degeneracy],
                # "vertex_num" => vertex_num,
                "nodes" => nodes,
                "edges" => edges,
                "work_nodes" => candidate
            ))
            found = true
            break 
        end
        !found && @info "No suitable graph found for degeneracy $(degeneracy) in Graph $(vertex_nums)"
    end
    filename = joinpath(save_path, "$(input_bits)in_$(output_bits)out.json")
    open(filename, "w") do file
        write(file, JSON.json(all_graph_data))
    end
    return filename
end

function search_any_constraint(vertex_nums::Vector{Int}, bit_num::Int, dir_path::String, save_path::String)
    all_graph_data = []
    for degeneracy in [collect(s) for s in IterTools.subsets(0:2^bit_num-1) if !isempty(s)]
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
                # "vertex_num" => vertex_num,
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
        maximalis = MaximalIS(graphs[gname])
        mis_problem = GenericTensorNetwork(maximalis)
        mis_result = read_config(solve(mis_problem, ConfigsAll())[])
        mis_num = length(mis_result)
        
        all_candidates = collect(permutations(1:vertex_num, bit_num))
        
        for candidate in all_candidates 
            candidate_value = [[Int(mis_result[i][j]) for j in candidate] for i in 1:mis_num]
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
                @constraint(model, sum((Int(mis[1][i]) - Int(mis[2][i])) * x[i] for i in 1:vertex_num) == 0)
            end 

            ϵ = -1
            for wrong_mis in wrong_mis_sets
                for mis in target_mis_sets
                    # @info "Adding constraint for $(wrong_mis)_$(mis)."
                    @constraint(model, sum((Int(mis[i]) - Int(wrong_mis[i])) * x[i] for i in 1:vertex_num) <= ϵ)
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