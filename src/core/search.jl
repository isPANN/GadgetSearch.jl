# using ProgressMeter

function search_over_dataset(loader::GraphLoader; filter=nothing, limit::Union{Int,Nothing}=nothing, keys_range::Union{Nothing, Vector{Int}}=nothing)
    keys_raw = keys_range === nothing ? keys(loader) : keys_range
    keys_to_search = isa(keys_raw[1], Int) ? keys_raw : parse.(Int, keys_raw)
    total = limit === nothing ? length(keys_to_search) : min(length(keys_to_search), limit)

    @showprogress for key in Iterators.take(keys_to_search, total)
        g = loader[key]
        if filter !== nothing && filter(g)
            return
        end
    end
end


# Create a filter closure for search_over_dataset
# function make_filter(truth_table::BitMatrix, bit_num::Int, params::SearchParameters, optimizer, env; connected::Bool=false)
  
#     return function(graph::SimpleGraph{Int}, pinset::Union{Nothing, Vector{Int}}=nothing)
#         if !connected
#             Graphs.is_connected(graph) || return false
#         end

#         vertex_num = Graphs.nv(graph)

#         mis_result, _ = find_maximal_independent_sets(graph)

#         if pinset === nothing
#             pin_set = collect(1:vertex_num)
#         end
#         @assert length(pin_set) >= bit_num

#         all_candidates = _generate_constraint_bit(pin_set, bit_num)

#         for candidate in all_candidates
#             target_mis_indices_all = _check_grstates_in_candidate(mis_result, candidate, ground_states, params.greedy)
#             isempty(target_mis_indices_all) && continue
#             weight = _find_weight_new(mis_result, target_mis_indices_all, optimizer, env)
#             isempty(weight) && continue
#             return true
#         end
#         return false
#     end
# end

# function find_maximal_independent_sets(g::SimpleGraph{Int})::Tuple{BitMatrix, Int}
#     compset = Graphs.complement(g)
#     cliques = Graphs.maximal_cliques(compset)
#     n = Graphs.nv(g)
#     result = falses(length(cliques), n)
#     for (i, clique) in enumerate(cliques)
#         for v in clique
#             result[i, v] = true
#         end
#     end
#     return result, length(cliques)
# end

function find_maximal_independent_sets(g::SimpleGraph{Int})
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

# function match_rows_by_pinset(mis_result::BitMatrix, truth_table::BitMatrix, pin_set::Vector{Int})
#     # 提取感兴趣列为普通数组，避免 BitMatrix 慢访问
#     projected = Matrix(mis_result[:, pin_set])
#     row_map = Dict{NTuple{length(pin_set), Bool}, Vector{Int}}()
#     for i in eachindex(axes(projected, 1))
#         key = ntuple(k -> projected[i, k], length(pin_set))
#         push!(get!(row_map, key, Int[]), i)
#     end

#     result = Vector{Vector{Int}}(undef, size(truth_table, 1))
#     for i in eachindex(axes(truth_table, 1))
#         key = ntuple(k -> truth_table[i, k], length(pin_set))
#         result[i] = get(row_map, key, Int[])
#     end
#     return result
# end

function match_rows_by_pinset(masks::Vector{UInt16}, truth_table::BitMatrix, pin_set::Vector{Int})
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