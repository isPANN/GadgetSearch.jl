struct Node{T}
    id::T
    weight::Float64
end

struct Edge{T}
    source::T
    target::T
end

struct Result{T}
    gate_id::Int
    degeneracy::Vector{String}
    nodes::Vector{GadgetSearch.Node{T}}
    edges::Vector{GadgetSearch.Edge{T}}
    work_nodes::Vector{T}
end

check_singleconstraint_for_graphs(vertex_num::Int, bit_num::Int, degeneracy::Vector{String}, dir_path::String) = check_singleconstraint_for_graphs(vertex_num, bit_num, format_degeneracy_input(degeneracy), dir_path)

check_singleconstraint_for_graphs(vertex_num::Int, bit_num::Int, degeneracy::Vector{Vector{Int}}, dir_path::String) = check_singleconstraint_for_graphs(vertex_num, bit_num, format_degeneracy_input(degeneracy), dir_path)

function check_singleconstraint_for_graphs(vertex_num::Int, bit_num::Union{Int, Vector{Int}}, degeneracy::Vector{Int}, dir_path::String; max_file_size::Int = 1_000_000)::Union{Result, Nothing}
    graph_path = joinpath(dir_path, "graph$(vertex_num).g6")
    if filesize(graph_path) > max_file_size
        process_large_file(graph_path, bit_num, degeneracy)
    else
        normally_process(graph_path, bit_num, degeneracy)
    end
end

function process_large_file(path::String, bit_num::Union{Int, Vector{Int}}, degeneracy::Vector{Int})
    open(path, "r") do io
        for line in eachline(io)
            graph_data = strip(line)
            isempty(graph_data) && continue
            g = g6string_to_graph(graph_data)
            candidate, weight = check_singleconstraint_for_singlegraph(g, bit_num, degeneracy)
            !isnothing(candidate) && return convert_to_result(g, candidate, weight, degeneracy)
        end
    end
end

function normally_process(path::String, bit_num::Union{Int, Vector{Int}}, degeneracy::Vector{Int})
    graph_list = read_g6_file(path)
    for g in graph_list
        candidate, weight = check_singleconstraint_for_singlegraph(g, bit_num, degeneracy)
        !isnothing(candidate) && return convert_to_result(g, candidate, weight, degeneracy)
    end
end

function check_singleconstraint_for_singlegraph(graph::SimpleGraph, bit_num::Int, degeneracy::Vector{Int}) 
    @assert length(degeneracy) > 0 && maximum(degeneracy) < 2^(sum(bit_num))
    # This is a version for a generic constraint.
    Graphs.is_connected(graph) || return nothing, nothing
    vertex_num = Graphs.nv(graph)
    mis_result = find_maximal_independent_sets(graph)
    all_candidates = _generate_candidates(bit_num, vertex_num)
    
    for candidate in all_candidates 
        target_mis_indices_all = _check_degeneracy_in_candidate(mis_result, candidate, degeneracy)
        isempty(target_mis_indices_all) && continue
        for combination in Iterators.product(target_mis_indices_all...)
            target_mis_set, wrong_mis_set = _split_mis_set(mis_result, vcat(combination...))
            weight = _find_weight(vertex_num, target_mis_set, wrong_mis_set)
            isempty(weight) && continue
            return candidate, weight
        end 
    end
    return nothing, nothing
end


function check_singleconstraint_for_singlegraph(graph::SimpleGraph, bit_num::Vector{Int}, degeneracy::Vector{Int}) 
    @assert length(degeneracy) > 0 && maximum(degeneracy) < 2^(sum(bit_num))
     # This is a version for a logic gate.
     @assert length(bit_num) == 2
     input_num = bit_num[1]; output_num = bit_num[2];
     Graphs.is_connected(graph) || return nothing, nothing
     vertex_num = Graphs.nv(graph)
     mis_result = find_maximal_independent_sets(graph)
     size(mis_result, 1) < 2^(input_num) && return nothing, nothing
     all_candidates = _generate_candidates(graph, mis_result, input_num)
     isempty(all_candidates) && return nothing, nothing
     
     for candidate in all_candidates 
         remain_elements = setdiff(1:vertex_num, candidate)
         for output_bits in permutations(remain_elements, output_num)
             candidate_full = vcat(candidate, output_bits)
             target_mis_indices_all = _check_degeneracy_in_candidate(mis_result, candidate_full, degeneracy)
             isempty(target_mis_indices_all) && continue
             for combination in Iterators.product(target_mis_indices_all...)
                 target_mis_set, wrong_mis_set = _split_mis_set(mis_result, vcat(combination...))
                 weight = _find_weight(vertex_num, target_mis_set, wrong_mis_set)
                 isempty(weight) && continue
                 return candidate, weight
             end
         end
     end
     return nothing, nothing 
end


function find_maximal_independent_sets(g::SimpleGraph)::AbstractMatrix{Int}
    ones_vertex = Graphs.maximal_cliques(Graphs.complement(g))
    mis_result = generate_bitvectors(Graphs.nv(g), ones_vertex)
    return mis_result
end

function convert_to_result(g::SimpleGraph, candidate::Vector{T}, weight::Vector{Float64}, degeneracy::Vector{Int}, gate_id::Int = 0) where T
    nodes = [Node(i, weight[i]) for i in 1:length(weight)]

    edges = [Edge(src(e), dst(e)) for e in Graphs.edges(g)]
    
    degeneracy_io = format_degeneracy_output(degeneracy, length(candidate))

    return Result{T}(
        gate_id,
        degeneracy_io,
        nodes,
        edges,
        candidate
    )
end

function save_results_to_json(results::Vector{Result{T}}, file_path::String) where T
    json_data = JSON3.write(results; pretty=true)
    open(file_path, "w") do io
        write(io, json_data)
    end
end

function _check_degeneracy_in_candidate(mis_result::AbstractMatrix{Int}, candidate::Vector{Int}, degeneracy::Vector{Int})::Vector{Vector{Int}}
    candidate_degeneracy = get_candidate_degeneracy(mis_result, candidate)
    is_subset = all(x -> x in candidate_degeneracy, degeneracy)
    is_subset || return Vector{Vector{Int}}()
    return [findall(==(x), candidate_degeneracy) for x in degeneracy]
end

function _split_mis_set(mis_result::AbstractMatrix{Int}, target_indices::Vector{Int})::Tuple{AbstractMatrix{Int}, AbstractMatrix{Int}}
    wrong_indices = setdiff(1:size(mis_result, 1), target_indices)
    target_set = mis_result[target_indices, :]
    wrong_set = mis_result[wrong_indices, :]
    return target_set, wrong_set
end

function _find_weight(vertex_num::Int, target_set::AbstractMatrix{Int}, wrong_set::AbstractMatrix{Int})
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, x[1:vertex_num])
    target_diffs = [target_set[i, :] .- target_set[j, :] for (i, j) in combinations(axes(target_set, 1), 2)]
    for diff in target_diffs
        @constraint(model, sum(diff[i] * x[i] for i in 1:vertex_num) == 0)
    end 

    ϵ = 1
    wrong_target_diffs = [
        target_set[target_mis_index, :] .- wrong_set[wrong_mis_index, :]
        for wrong_mis_index in axes(wrong_set, 1)
        for target_mis_index in axes(target_set, 1)
    ]
    
    for diff in wrong_target_diffs
        @constraint(model, sum(diff[i] * x[i] for i in 1:vertex_num) >= ϵ)
    end

    for i in 1:vertex_num
        @constraint(model, 1 <= x[i])
    end
    
    @objective(model, Min, sum(x[i] for i in 1:vertex_num))
    optimize!(model)
    
    if is_solved_and_feasible(model)
        @info "Optimization successful!"
        return [value(x[i]) for i in 1:vertex_num]
    end
   
    for con in all_constraints(model; include_variable_in_set_constraints = false)
        delete(model, con)
    end
    return []
end

_generate_candidates(bit_num::Int, vertex_num::Int)::Vector{Vector{Int}} = collect(permutations(1:vertex_num, bit_num))

function _generate_candidates(g::SimpleGraph, mis::AbstractMatrix{Int}, input_num::Int)::Vector{Vector{Int}}
    col_sums = sum(mis, dims=1)
    valid_columns = findall(x -> x >= 2^(input_num - 1), col_sums[:])
    length(valid_columns) < input_num && return []
    sets = Vector{Vector{Int}}()
    explicit_combinations = collect(combinations(valid_columns, input_num))

    for subset in explicit_combinations
        is_independent = true
        for (u, v) in combinations(subset, 2)
            has_edge(g, u, v) && (is_independent = false; break)
        end
        if is_independent
            push!(sets, collect(subset))
        end
    end
    return sets
end



