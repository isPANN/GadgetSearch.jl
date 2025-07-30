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
    limit=nothing,
    max_samples::Int=100,
    save_path::String="results.json",
    pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing
)
    results = Gadget[]
    failed_tt = BitMatrix[]

    for (i, tt) in enumerate(truth_tables)
        @info "searching for truth table: $(i-1)"
        filter_fn = make_filter(tt, bit_num, optimizer, env;
                                connected=connected,
                                objective=objective,
                                allow_defect=allow_defect,
                                max_samples=max_samples,
                                pin_candidates=pin_candidates)
        gadgets = find_matching_gadget(loader; filter=filter_fn, limit=limit)
        if !isempty(gadgets)
            push!(results, gadgets...)
            save_results_to_json(results, save_path)
        else
            push!(failed_tt, tt)
        end
    end

    return results, failed_tt
end


function find_matching_gadget(loader::GraphLoader; filter=nothing, limit::Union{Int,Nothing}=nothing, keys_range::Union{Nothing, Vector{Int}}=nothing)
    keys_raw = keys_range === nothing ? collect(keys(loader)) : keys_range
    keys_to_search = isa(keys_raw[1], Int) ? keys_raw : parse.(Int, keys_raw)
    total = limit === nothing ? length(keys_to_search) : min(length(keys_to_search), limit)

    results = Gadget[]

    @showprogress for key in Iterators.take(keys_to_search, total)
        g = loader[key]
        if filter !== nothing
            weights, truth_table, pin = filter(g, loader.layout[key], loader.pinset)
            if !isnothing(weights)
                # return Gadget(truth_table, g, pin, weights, loader.layout[key])
                push!(results, Gadget(truth_table, g, pin, weights, loader.layout[key]))
            end
        end
    end
    return results
end

# Create a filter closure for search_over_dataset
function make_filter(truth_table::BitMatrix, bit_num::Union{Int, Tuple{Int, Int}}, optimizer, env; connected::Bool=false, objective=nothing, allow_defect::Bool=false, max_samples::Int=1000, pin_candidates::Union{Nothing, Vector{Vector{Int}}}=nothing)
  
    return function(graph::SimpleGraph{Int}, pos::Vector{Tuple{Float64, Float64}}, pin_set::Union{Nothing, Vector{Int}}=nothing)
        if !connected
            Graphs.is_connected(graph) || return nothing, truth_table, nothing
        end

        vertex_num = Graphs.nv(graph)

        mis_result, count = find_maximal_independent_sets(graph)

        if pin_set === nothing
            pin_set = collect(1:vertex_num)
        end
        
        if pin_candidates !== nothing
            all_candidates = Vector{Vector{Int}}()
            for candidate in pin_candidates
                if all(pin in pin_set for pin in candidate)
                    push!(all_candidates, candidate)
                else
                    @warn "Candidate $candidate contains pins not in pin_set $pin_set, skipping"
                end
            end
        else
            all_candidates = generate_pin_variants(graph, pos, pin_set, bit_num)
        end
        
        for candidate in all_candidates
            target_mis_indices_all = match_rows_by_pinset(mis_result, truth_table, candidate)
            weights = solve_weight_enumerate(mis_result, target_mis_indices_all, vertex_num, optimizer, env, objective, allow_defect, max_samples)
            if !isempty(weights)
                return weights, truth_table, candidate
            end
        end
        return nothing, truth_table, nothing
    end
end

function generate_pin_variants(
    g::SimpleGraph{Int},
    pos::Vector{Tuple{Float64, Float64}},
    pin_set::Vector{Int},
    bit_num::Union{Int, Tuple{Int, Int}}
)
    if isa(bit_num, Int)
        return collect(Combinatorics.permutations(pin_set, bit_num))
    else
        total = Vector{Vector{Int}}()

        for comb in Combinatorics.permutations(pin_set, bit_num[1]+bit_num[2])
            has_conflict = any(has_edge(g, u, v) for (u, v) in combinations(comb, 2))
            # @show comb
            if !has_conflict && length(comb) ≥ 4
                # 判断线段是否相交：1-3 和 2-4
                p1, p2 = pos[comb[1]], pos[comb[3]]
                q1, q2 = pos[comb[2]], pos[comb[4]]
                if segments_intersect(p1, p2, q1, q2)
                    push!(total, comb)
                end
            end
        end

        return total
        # input_len, output_len = bit_num

        # for input in Combinatorics.combinations(pin_set, input_len)
        #     rest = setdiff(pin_set, input)
        #     if length(rest) ≥ output_len
        #         for output in permutations(rest, output_len)
        #             combined = vcat(input, output)

        #             has_conflict = any(has_edge(g, u, v) for (u, v) in combinations(combined, 2))

        #             if !has_conflict && length(combined) ≥ 4
        #                 # # 判断线段是否相交：1-3 和 2-4
        #                 # p1, p2 = pos[combined[1]], pos[combined[3]]
        #                 # q1, q2 = pos[combined[2]], pos[combined[4]]
        #                 # if segments_intersect(p1, p2, q1, q2)
        #                 #     push!(total, combined)
        #                 # end
        #                 push!(total, combined)
        #             end
                # end
            # end
        # end
        return total
    end
end

function segments_intersect(p1, p2, q1, q2)::Bool
    # 向量叉积判断顺时针 / 逆时针
    function ccw(a, b, c)
        return (c[2] - a[2]) * (b[1] - a[1]) > (b[2] - a[2]) * (c[1] - a[1])
    end
    return (ccw(p1, q1, q2) != ccw(p2, q1, q2)) && (ccw(p1, p2, q1) != ccw(p1, p2, q2))
end

function find_maximal_independent_sets(g::SimpleGraph{Int})
    cliques = Graphs.maximal_cliques(Graphs.complement(g))
    vertex_count = Graphs.nv(g)
    
    # 检查顶点数限制（UInt32支持32个顶点）
    if vertex_count > 32
        error("Graph has $(vertex_count) vertices, but maximum supported is 32. Consider reducing graph size or using matrix-based algorithms.")
    end
    
    # 使用UInt32进行高效位运算
    masks = UInt32[]
    sizehint!(masks, length(cliques))  # 预分配内存
    
    for clique in cliques
        mask::UInt32 = 0
        for v in clique
            mask |= UInt32(1) << (v - 1)
        end
        push!(masks, mask)
    end
    
    return masks, length(masks)
end

function match_rows_by_pinset(masks::Vector{UInt32}, truth_table::BitMatrix, pin_set::Vector{Int})::Vector{Vector{Int}}
    num_rows = size(truth_table, 1)
    result = Vector{Vector{Int}}(undef, num_rows)

    for i in 1:num_rows
        # 构建查询掩码
        query_mask::UInt32 = 0
        for (bit_pos, pin) in enumerate(pin_set)
            query_mask |= UInt32(truth_table[i, bit_pos]) << (bit_pos - 1)
        end
        
        # 查找匹配的MIS
        matches = Int[]
        for (j, m) in enumerate(masks)
            # 从MIS中提取pin位置的值
            extracted::UInt32 = 0
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
    mis_result::Vector{UInt32},
    target_mis_indices_all::Vector{Vector{Int}},
    vertex_num::Int,
    optimizer,
    env=nothing,
    objective=nothing,
    allow_defect::Bool=false,
    max_samples::Int=1000
)
    if optimizer === nothing
        error("Optimizer must be provided.")
    end

    any(isempty, target_mis_indices_all) && return Float64[]
    all_mis_indices = eachindex(mis_result)
    # filtered_mis = select_low_energy_mis(mis_result, target_mis_indices_all, vertex_num)

    # log_total = sum(log(length(indices)) for indices in filtered_mis)
    # if log_total <= log(max_samples)
    #     combinations = collect(Iterators.product(filtered_mis...))
    # else
    #     samples = []
    #     for _ in 1:max_samples
    #         sample = [rand(indices) for indices in filtered_mis]
    #         push!(samples, sample)
    #     end
    #     combinations = samples
    # end
    # for target_indices in combinations
    # for target_indices in collect(Iterators.product(target_mis_indices_all...))
    
    # 枚举所有组合来寻找权重
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

    return Float64[]
end

function _find_weight(vertex_num::Int, target_masks::Vector{UInt32}, wrong_masks::Vector{UInt32}, optimizer, env, objective, allow_defect::Bool)
    opt = isnothing(env) ? optimizer() : optimizer(env)
    model = direct_model(opt)
    set_silent(model)
    set_string_names_on_creation(model, false)

    if allow_defect
        x = @variable(model, [i=1:vertex_num], lower_bound = i <= 4 ? 1 : 0)
    else
        @variable(model, x[1:vertex_num] >= 1)
    end

    @variable(model, C)

    # 目标MIS必须有相同的能量
    for m in target_masks
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) == C)
    end

    # 错误的MIS必须有更高的能量
    ϵ = 1.0
    for m in wrong_masks
        @constraint(model, sum(((m >> (v - 1)) & 0x1) * x[v] for v in 1:vertex_num) <= C - ϵ)
    end

    if objective !== nothing
        @objective(model, Min, objective(x))
    end

    optimize!(model)
    if is_solved_and_feasible(model)
        @info "found a solution"
        return [value(x[v]) for v in 1:vertex_num]
    end
    return Float64[]
end