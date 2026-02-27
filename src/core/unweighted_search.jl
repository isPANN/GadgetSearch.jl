"""
    UnweightedGadget

Result of an unweighted gadget search (Definition 3.6 / Theorem 3.7 in the paper).

# Fields
- `pattern_graph::SimpleGraph{Int}`: The target graph R
- `replacement_graph::SimpleGraph{Int}`: The candidate graph R' that replaces R
- `boundary_vertices::Vector{Int}`: Boundary vertices of `replacement_graph`
- `constant_c::Float64`: Offset `α̃(R') - α̃(R)` (the MIS size difference introduced by the replacement)
"""
struct UnweightedGadget
    pattern_graph::SimpleGraph{Int}
    replacement_graph::SimpleGraph{Int}
    boundary_vertices::Vector{Int}
    constant_c::Float64
end

"""
    make_unweighted_filter(target_graph, target_boundary)

Create a filter closure for unweighted gadget search.

Pre-computes `α̃(target_graph)` once, then returns a closure that checks each
candidate graph by trying all boundary vertex combinations of the same size.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R

# Returns
- `Function`: Closure `(candidate, pos, pin_set) -> UnweightedGadget | nothing`
"""
function make_unweighted_filter(target_graph::SimpleGraph{Int}, target_boundary::Vector{Int})
    target_reduced = content.(calculate_reduced_alpha_tensor(target_graph, target_boundary))
    k = length(target_boundary)

    return function(candidate::SimpleGraph{Int}, pos, pin_set)
        nv(candidate) < k && return nothing
        candidates = pin_set !== nothing ? collect(Combinatorics.combinations(pin_set, k)) :
                                          collect(Combinatorics.combinations(1:nv(candidate), k))
        for boundary in candidates
            candidate_reduced = content.(calculate_reduced_alpha_tensor(candidate, boundary))
            valid, c = is_diff_by_constant(candidate_reduced, target_reduced)
            if valid && !isnan(c)
                return UnweightedGadget(target_graph, candidate, boundary, Float64(c))
            end
        end
        return nothing
    end
end

"""
    search_unweighted_gadgets(target_graph, target_boundary, loader; limit, max_results)

Search for unweighted gadget replacements of `target_graph` by iterating over a `GraphLoader`.

For each candidate graph, tries all boundary vertex combinations of size
`length(target_boundary)` and checks if the reduced alpha tensors differ by a constant
(Theorem 3.7). Returns on the first valid boundary found per candidate.

# Arguments
- `target_graph::SimpleGraph{Int}`: The pattern graph R
- `target_boundary::Vector{Int}`: Boundary vertices of R
- `loader::GraphLoader`: Graph dataset to search over

# Keywords
- `limit::Union{Int,Nothing}=nothing`: Maximum number of graphs to examine
- `max_results::Union{Int,Nothing}=nothing`: Stop after finding this many results

# Returns
- `Vector{UnweightedGadget}`
"""
function search_unweighted_gadgets(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    loader::GraphLoader;
    limit::Union{Int,Nothing}=nothing,
    max_results::Union{Int,Nothing}=nothing
)
    filter_fn = make_unweighted_filter(target_graph, target_boundary)
    results = UnweightedGadget[]
    total = limit === nothing ? length(loader) : min(length(loader), limit)

    @showprogress for key in Iterators.take(keys(loader), total)
        result = filter_fn(loader[key], loader.layout[key], loader.pinset)
        if result !== nothing
            push!(results, result)
            if max_results !== nothing && length(results) >= max_results
                break
            end
        end
    end
    return results
end

