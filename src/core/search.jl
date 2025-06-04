# Search for multiple truth tables by reusing search_over_dataset
function search_by_truth_tables(
    loader::GraphLoader,
    truth_tables::Vector{BitMatrix};
    bit_num,
    optimizer,
    env=nothing,
    connected::Bool=false,
    objective=nothing,
    allow_defect::Bool=false,
    limit=nothing
)
    results = Dict{BitMatrix, Gadget}()

    for tt in truth_tables
        filter_fn = make_filter(tt, bit_num, optimizer, env;
                                connected=connected,
                                objective=objective,
                                allow_defect=allow_defect)
        gadget = find_matching_gadget(loader; filter=filter_fn, limit=limit)
        if gadget !== nothing
            results[tt] = gadget
        end
    end

    return results
end

# Alternative: iterate over each graph once and test multiple truth tables per graph
function search_graphs_for_truth_tables(
    loader::GraphLoader,
    truth_tables::Vector{BitMatrix};
    bit_num,
    optimizer,
    env=nothing,
    connected::Bool=false,
    objective=nothing,
    allow_defect::Bool=false,
    limit=nothing
)
    results = Dict{String, Tuple{BitMatrix, Gadget}}()

    keys_to_search = keys(loader)
    if limit !== nothing
        keys_to_search = Iterators.take(keys_to_search, limit)
    end

    @showprogress for key in keys_to_search
        g = loader[key]
        connected && !Graphs.is_connected(g) && continue
        pin_set = loader.pinset

        vertex_num = Graphs.nv(g)
        mis_result = find_maximal_independent_sets(g)

        for tt in truth_tables
            all_candidates = generate_pin_variants(pin_set, bit_num)
            for candidate in all_candidates
                target_indices = match_rows_by_pinset(mis_result, tt, candidate)
                weights = solve_weight_enumerate(mis_result, target_indices, vertex_num, optimizer, env, objective, allow_defect)
                if !isempty(weights)
                    results[key] = (tt, Gadget(tt, g, candidate, weights, loader.layout[key]))
                    break
                end
            end
            haskey(results, key) && break
        end
    end
    return results
end

function find_matching_gadget(loader::GraphLoader; filter=nothing, limit::Union{Int,Nothing}=nothing, keys_range::Union{Nothing, Vector{Int}}=nothing)
    keys_raw = keys_range === nothing ? keys(loader) : keys_range
    keys_to_search = isa(keys_raw[1], Int) ? keys_raw : parse.(Int, keys_raw)
    total = limit === nothing ? length(keys_to_search) : min(length(keys_to_search), limit)

    @showprogress for key in Iterators.take(keys_to_search, total)
        g = loader[key]
        if filter !== nothing
            weights, truth_table, pin = filter(g, loader.pinset)
            if !isnothing(weights)
                return Gadget(truth_table, g, pin, weights, loader.layout[key])
            end
        end
    end
    return nothing
end

# Create a filter closure for search_over_dataset
function make_filter(truth_table::BitMatrix, bit_num::Union{Int, Tuple{Int, Int}}, optimizer, env; connected::Bool=false, objective=nothing, allow_defect::Bool=false)
  
    return function(graph::SimpleGraph{Int}, pin_set::Union{Nothing, Vector{Int}}=nothing)
        if !connected
            Graphs.is_connected(graph) || return nothing, truth_table, nothing
        end

        vertex_num = Graphs.nv(graph)

        mis_result = find_maximal_independent_sets(graph)

        if pin_set === nothing
            pin_set = collect(1:vertex_num)
        end
        
        all_candidates = generate_pin_variants(pin_set, bit_num)

        for candidate in all_candidates
            target_mis_indices_all = match_rows_by_pinset(mis_result, truth_table, candidate)
            weights = solve_weight_enumerate(mis_result, target_mis_indices_all, vertex_num, optimizer, env, objective, allow_defect)
            if !isempty(weights)
                return weights, truth_table, candidate
            end
        end
        return nothing, truth_table, nothing
    end
end

function generate_pin_variants(pin_set::Vector{Int}, bit_num::Union{Int, Tuple{Int, Int}})
    if isa(bit_num, Int)
        return collect(Combinatorics.permutations(pin_set, bit_num))
    else
        input_len, output_len = bit_num
        total = Vector{Vector{Int}}()

        for input in Combinatorics.combinations(pin_set, input_len)
            rest = setdiff(pin_set, input)
            for output in Combinatorics.permutations(rest, output_len)
                # pre-allocate the combined vector to avoid push! 
                combined = Vector{Int}(undef, input_len + output_len)
                copyto!(combined, 1, input, 1, input_len)
                copyto!(combined, input_len + 1, output, 1, output_len)
                push!(total, combined)
            end
        end
        return total
    end
end


function find_maximal_independent_sets(g::SimpleGraph{Int})::Vector{UInt16}
    cliques = Graphs.maximal_cliques(Graphs.complement(g))
    masks = UInt16[]
    for clique in cliques
        mask::UInt16 = 0
        for v in clique
            mask |= UInt16(1) << (v - 1)
        end
        push!(masks, mask)
    end
    return masks
end


function match_rows_by_pinset(masks::Vector{UInt16}, truth_table::BitMatrix, pin_set::Vector{Int})::Vector{Vector{Int}}
    result = Vector{Vector{Int}}(undef, size(truth_table, 1))
    for i in eachindex(axes(truth_table, 1))
        query_mask::UInt16 = 0
        for (bit_pos, pin) in enumerate(pin_set)
            query_mask |= UInt16(truth_table[i, bit_pos]) << (bit_pos - 1)
        end
        matches = Int[]
        for (j, m) in enumerate(masks)
            extracted::UInt16 = 0
            for (bit_pos, pin) in enumerate(pin_set)
                extracted |= ((m >> (pin - 1)) & 0x1) << (bit_pos - 1)
            end
            if extracted == query_mask
                push!(matches, j)
            end
        end
        result[i] = matches
    end
    return result
end

function solve_weight_enumerate(
    mis_result::Vector{UInt16},
    target_mis_indices_all::Vector{Vector{Int}},
    vertex_num::Int,
    optimizer,
    env=nothing,
    objective=nothing,
    allow_defect::Bool=false
)
    if optimizer === nothing
        error("Optimizer must be provided.")
    end

    any(isempty, target_mis_indices_all) && return Float64[]
    all_mis_indices = eachindex(mis_result)

    for target_indices in Iterators.product(target_mis_indices_all...)
        target_indices_set = collect(target_indices)
        wrong_indices = setdiff(all_mis_indices, target_indices_set)

        target_set = mis_result[target_indices_set]
        wrong_set = mis_result[wrong_indices]

        weights = _find_weight(vertex_num, target_set, wrong_set, optimizer, env, objective, allow_defect)
        if !isempty(weights)
            return weights
        end
    end
    # weights = _find_weight_mask_new(mis_result, target_mis_indices_all, vertex_num; optimizer=optimizer, env=env, objective=objective)
    # if !isempty(weights)
    #     return weights
    # end
    return Float64[]
end

function _find_weight(vertex_num::Int, target_masks::Vector{UInt16}, wrong_masks::Vector{UInt16}, optimizer, env, objective, allow_defect::Bool)
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)

    if allow_defect
        @variable(model, x[1:vertex_num] >= 0, Int)
    else
        @variable(model, x[1:vertex_num] >= 1, Int)
    end

    @variable(model, C)

    for m in target_masks
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) == C)
    end

    ϵ = 1.0
    for m in wrong_masks
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) <= C - ϵ)
    end

    if objective !== nothing
        @objective(model, Min, objective(x))
    end

    optimize!(model)
    if is_solved_and_feasible(model)
        return [value(x[v]) for v in 1:vertex_num]
    end
    return Float64[]
end

# function _find_weight_mask_new(
#     mis_masks::Vector{UInt16},
#     prefix_buckets::Vector{Vector{Int}},
#     vertex_num::Int;
#     optimizer,
#     env=nothing,
#     objective=nothing
# )
#     opt = isnothing(env) ? optimizer() : optimizer(env)
#     model = direct_model(opt)
#     set_silent(model)
#     set_string_names_on_creation(model, false)

#     y_dim = length(mis_masks)

#     @variable(model, x[1:vertex_num] >= 1, Int)
#     @variable(model, y[1:y_dim], Bin)

#     selected_indices = reduce(vcat, prefix_buckets)
#     selected_set = Set(selected_indices)

#     for idxs in prefix_buckets
#         @constraint(model, sum(y[i] for i in idxs) == 1)
#     end

#     @constraint(model, [i=1:y_dim; i ∉ selected_set], y[i] == 0)

#     row_weights = [count_ones(m) for m in mis_masks]
#     M_big = 2 * maximum(row_weights)

#     reference_indices = [first(idxs) for idxs in prefix_buckets]

#     for i in selected_indices
#         for ref_idx in reference_indices
#             @constraint(model,
#                 sum(((mis_masks[i] >> (v - 1)) & 1) * x[v] for v in 1:vertex_num) -
#                 sum(((mis_masks[ref_idx] >> (v - 1)) & 1) * x[v] for v in 1:vertex_num)
#                 <= M_big * (2 - y[i] - y[ref_idx])
#             )
#             @constraint(model,
#                 sum(((mis_masks[ref_idx] >> (v - 1)) & 1) * x[v] for v in 1:vertex_num) -
#                 sum(((mis_masks[i] >> (v - 1)) & 1) * x[v] for v in 1:vertex_num)
#                 <= M_big * (2 - y[i] - y[ref_idx])
#             )
#         end
#     end

#     for i in 1:y_dim
#         for j in selected_indices
#             @constraint(model,
#                 sum(((mis_masks[i] >> (v - 1)) & 1) * x[v] for v in 1:vertex_num)
#                 <= sum(((mis_masks[j] >> (v - 1)) & 1) * x[v] for v in 1:vertex_num)
#                 - 1 + M_big * (1 - y[j] + y[i])
#             )
#         end
#     end

#     if objective !== nothing
#         @objective(model, Min, objective(x))
#     end

#     set_time_limit_sec(model, 30)
#     optimize!(model)

#     if is_solved_and_feasible(model)
#         return [value(x[v]) for v in 1:vertex_num]
#     end
#     return Float64[]
# end