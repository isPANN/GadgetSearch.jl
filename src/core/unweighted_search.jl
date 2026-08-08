# ============================================================================
# Unweighted lattice search
# ============================================================================

const _LatticeCoordinate = Tuple{Int, Int}

"""A verifier-accepted replacement together with its concrete lattice embedding."""
struct UnweightedGadget
    pattern_graph::SimpleGraph{Int}
    replacement_graph::SimpleGraph{Int}
    boundary_vertices::Vector{Int}
    constant_offset::Float64
    lattice::Symbol
    lattice_coordinates::Vector{_LatticeCoordinate}
    pos::Vector{Tuple{Float64, Float64}}
    pin_rays::Vector{_LatticeCoordinate}
end

"""One evaluated transition of the dynamically constructed lattice patch."""
struct UnweightedSearchRecord
    generation::Int
    stage::Symbol
    key::String
    lattice::Symbol
    graph_edges::Vector{Tuple{Int, Int}}
    lattice_coordinates::Union{Nothing, Vector{_LatticeCoordinate}}
    pin_coordinates::Union{Nothing, Vector{_LatticeCoordinate}}
    pin_rays::Vector{_LatticeCoordinate}
    boundary_vertices::Vector{Int}
    parent_key::Union{Nothing, String}
    action::Symbol
    vertices::Int
    edges::Int
    mask_mismatches::Int
    offset_spread::Float64
    constraint_defects::Int
    embedding_placed::Int
    rewrite_steps::Int
    is_solution::Bool
    constant_offset::Union{Nothing, Float64}
    selected::Bool
end

"""Outcome of a bounded dynamic search on one concrete lattice."""
struct UnweightedSearchResult
    target_graph::SimpleGraph{Int}
    target_boundary::Vector{Int}
    lattice::Symbol
    gadgets::Vector{UnweightedGadget}
    evaluated::Int
    generations::Int
    best_mask_mismatches::Int
    best_offset_spread::Float64
    termination_reason::Symbol
    trace::Vector{UnweightedSearchRecord}
end

struct _LatticePatch
    coordinates::Vector{_LatticeCoordinate}
    pins::Vector{_LatticeCoordinate}
    rays::Vector{Int}
end

struct _GraphState
    graph::SimpleGraph{Int}
    boundary::Vector{Int}
    parent_key::Union{Nothing, String}
    action::Symbol
    rewrite_steps::Int
end

"""Search logical skeletons, rewrite them exactly, then embed them on `lattice`."""
function search_unweighted_gadgets(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    lattice::LatticeType;
    min_vertices::Int=length(target_boundary) + 1,
    max_vertices::Int=min_vertices + 8,
    max_evaluations::Int=2_000,
    beam_width::Int=32,
    mutations_per_candidate::Int=8,
    random_candidates_per_generation::Int=8,
    exploration_fraction::Float64=0.25,
    max_results::Int=1,
    rng::AbstractRNG=Random.default_rng(),
)
    boundary_count = length(target_boundary)
    min_vertices >= boundary_count || throw(ArgumentError("min_vertices must be at least the number of boundary vertices"))
    max_vertices >= min_vertices || throw(ArgumentError("max_vertices must be at least min_vertices"))
    max_evaluations > 0 || throw(ArgumentError("max_evaluations must be positive"))
    beam_width > 0 || throw(ArgumentError("beam_width must be positive"))
    mutations_per_candidate > 0 || throw(ArgumentError("mutations_per_candidate must be positive"))
    random_candidates_per_generation > 0 || throw(ArgumentError("random_candidates_per_generation must be positive"))
    0.0 <= exploration_fraction < 1.0 || throw(ArgumentError("exploration_fraction must be in [0, 1)"))
    max_results > 0 || throw(ArgumentError("max_results must be positive"))

    target_reduced = vec(calculate_reduced_alpha_tensor(target_graph, target_boundary))
    all(isinf, target_reduced) && error("target graph has an entirely -Inf reduced alpha tensor")
    logical_budget = max(1, max_evaluations ÷ 2)
    skeletons, evaluated, generations, best_score, trace = _search_logical_skeletons(
        target_graph, target_boundary, target_reduced, lattice;
        min_vertices, max_vertices, max_evaluations=logical_budget, beam_width,
        mutations_per_candidate, random_candidates_per_generation,
        exploration_fraction, rng,
    )
    gadgets, rewrite_evaluations, rewrite_generations = _rewrite_and_embed_skeletons(
        target_graph, target_boundary, skeletons, lattice;
        max_vertices, max_evaluations=max_evaluations - evaluated, beam_width,
        mutations_per_candidate, exploration_fraction, max_results, rng, trace,
    )
    evaluated += rewrite_evaluations
    generations += rewrite_generations
    reason = !isempty(gadgets) ? :solution : evaluated == max_evaluations ? :budget : :search_space_exhausted
    return UnweightedSearchResult(
        target_graph, copy(target_boundary), _lattice_symbol(lattice), gadgets,
        evaluated, generations, best_score[1], best_score[2], reason, trace,
    )
end

_graph_edges(graph) = [(src(edge), dst(edge)) for edge in edges(graph)]
function _graph_state_key(state::_GraphState)
    return "n=$(nv(state.graph));pins=$(join(state.boundary,','));edges=$(join(("$a-$b" for (a, b) in _graph_edges(state.graph)),','))"
end

function _random_logical_state(rng, boundary_count, min_vertices, max_vertices)
    while true
        vertex_count = rand(rng, min_vertices:min(max_vertices, min_vertices + 8))
        graph = SimpleGraph(vertex_count)
        for vertex in 2:vertex_count
            add_edge!(graph, vertex, rand(rng, 1:vertex-1))
        end
        for first in 1:vertex_count-1, second in first+1:vertex_count
            rand(rng) < 0.16 && add_edge!(graph, first, second)
        end
        boundary_count == 4 && !_has_alternating_planar_frame(graph, 1:4) && continue
        return _GraphState(graph, collect(1:boundary_count), nothing, :logical_seed, 0)
    end
end

function _mutate_logical_state(rng, state, min_vertices, max_vertices)
    boundary_count = length(state.boundary)
    for _ in 1:32
        graph = deepcopy(state.graph)
        actions = Symbol[]
        nv(graph) >= 2 && push!(actions, :toggle_edges)
        nv(graph) < max_vertices && push!(actions, :add_vertex)
        nv(graph) > max(min_vertices, boundary_count) && push!(actions, :remove_vertex)
        action = rand(rng, actions)
        if action == :add_vertex
            add_vertex!(graph)
            degree = rand(rng, 1:min(4, nv(graph) - 1))
            for neighbor in randperm(rng, nv(graph) - 1)[1:degree]
                add_edge!(graph, nv(graph), neighbor)
            end
        elseif action == :remove_vertex
            rem_vertex!(graph, rand(rng, boundary_count+1:nv(graph)))
        else
            for _ in 1:rand(rng, 1:8)
                first, second = randperm(rng, nv(graph))[1:2]
                has_edge(graph, first, second) ? rem_edge!(graph, first, second) : add_edge!(graph, first, second)
            end
        end
        boundary_count == 4 && !_has_alternating_planar_frame(graph, state.boundary) && continue
        return _GraphState(graph, copy(state.boundary), _graph_state_key(state), action, state.rewrite_steps)
    end
    return state
end

function _has_alternating_planar_frame(graph, boundary)
    augmented = deepcopy(graph)
    for index in eachindex(boundary)
        add_edge!(augmented, boundary[index], boundary[mod1(index + 1, 4)])
    end
    add_vertex!(augmented)
    center = nv(augmented)
    for pin in boundary
        add_edge!(augmented, pin, center)
    end
    return is_planar(augmented)
end

function _evaluate_graph_state(state, target_reduced)
    raw, reduced, optimum_counts, near_counts = _search_alpha_tensors(state.graph, state.boundary)
    valid, offset = is_diff_by_constant(reduced, target_reduced)
    differences = [candidate - target for (candidate, target) in zip(reduced, target_reduced)
        if isfinite(candidate) && isfinite(target)]
    score = (
        count(isinf(candidate) != isinf(target) for (candidate, target) in zip(reduced, target_reduced)),
        Float64(maximum(differences) - minimum(differences)),
        length(connected_components(state.graph)) - 1,
        _tensor_plateau_signal(optimum_counts, near_counts, reduced, target_reduced),
        maximum(degree(state.graph)),
    )
    signature = raw .- raw[1]
    return score, signature, valid, Float64(offset), raw, reduced
end

function _tensor_plateau_signal(optimum_counts, near_counts, reduced, target)
    unwanted = sum((optimum_counts[index] for index in eachindex(reduced)
        if isfinite(reduced[index]) && !isfinite(target[index])); init=0)
    missing = sum((near_counts[index] for index in eachindex(reduced)
        if !isfinite(reduced[index]) && isfinite(target[index])); init=0)
    return unwanted - missing
end

function _search_alpha_tensors(graph, boundary)
    vertex_count = nv(graph)
    vertex_count <= 20 || error("logical skeleton search supports at most 20 vertices")
    adjacency = fill(UInt64(0), vertex_count)
    for edge in edges(graph)
        adjacency[src(edge)] |= UInt64(1) << (dst(edge) - 1)
        adjacency[dst(edge)] |= UInt64(1) << (src(edge) - 1)
    end
    raw = fill(-Inf, 1 << length(boundary))
    optimum_counts = zeros(Int, length(raw))
    one_conflict_counts = zeros(Int, length(raw), vertex_count + 1)
    for occupied in UInt64(0):(UInt64(1) << vertex_count)-1
        remaining = occupied
        conflicts = 0
        while remaining != 0
            vertex = trailing_zeros(remaining) + 1
            remaining &= remaining - 1
            conflicts += count_ones(adjacency[vertex] & remaining)
            conflicts > 1 && break
        end
        state = sum(((occupied >> (vertex - 1)) & 1) << (slot - 1) for (slot, vertex) in enumerate(boundary))
        size = count_ones(occupied)
        if conflicts == 0 && size > raw[state+1]
            raw[state+1] = size
            optimum_counts[state+1] = 1
        elseif conflicts == 0 && size == raw[state+1]
            optimum_counts[state+1] += 1
        elseif conflicts == 1
            one_conflict_counts[state+1, size+1] += 1
        end
    end
    tensor = reshape(Tropical.(raw), ntuple(_ -> 2, length(boundary)))
    reduced = vec(Float64.(content.(mis_compactify!(tensor))))
    near_counts = [isfinite(raw[state]) && raw[state] < vertex_count ?
        one_conflict_counts[state, Int(raw[state]) + 2] : 0 for state in eachindex(raw)]
    return raw, reduced, optimum_counts, near_counts
end

function _tensor_repair_states(rng, state, target_reduced, raw, reduced, limit)
    wrong = Set(index for index in eachindex(reduced)
        if isinf(reduced[index]) != isinf(target_reduced[index]))
    isempty(wrong) && return _GraphState[]
    toggles = Set{Tuple{Symbol, Tuple{Int, Int}}}()
    for occupied in UInt64(0):(UInt64(1) << nv(state.graph))-1
        index = _boundary_state(occupied, state.boundary) + 1
        index in wrong || continue
        if isfinite(target_reduced[index])
            wanted_size = isfinite(raw[index]) ? Int(raw[index]) + 1 : count_ones(index - 1)
            count_ones(occupied) == wanted_size || continue
            conflicts = [(src(edge), dst(edge)) for edge in edges(state.graph)
                if occupied & (UInt64(1) << (src(edge) - 1)) != 0 &&
                    occupied & (UInt64(1) << (dst(edge) - 1)) != 0]
            length(conflicts) == 1 && push!(toggles, (:remove_edge, only(conflicts)))
        else
            count_ones(occupied) == raw[index] || continue
            _is_independent_mask(state.graph, occupied) || continue
            selected = [vertex for vertex in vertices(state.graph)
                if occupied & (UInt64(1) << (vertex - 1)) != 0]
            for first in 1:length(selected)-1, second in first+1:length(selected)
                edge = minmax(selected[first], selected[second])
                !has_edge(state.graph, edge...) && push!(toggles, (:add_edge, edge))
            end
        end
    end
    ordered_toggles = shuffle!(rng, collect(toggles))
    children = _GraphState[]
    for (action, edge) in ordered_toggles[1:min(length(ordered_toggles), cld(limit, 2))]
        graph = deepcopy(state.graph)
        action == :add_edge ? add_edge!(graph, edge...) : rem_edge!(graph, edge...)
        length(state.boundary) == 4 && !_has_alternating_planar_frame(graph, state.boundary) && continue
        push!(children, _GraphState(
            graph, copy(state.boundary), _graph_state_key(state), action, state.rewrite_steps,
        ))
    end
    for _ in 1:8limit
        (length(children) == limit || length(ordered_toggles) < 2) && break
        graph = deepcopy(state.graph)
        count = rand(rng, 2:min(6, length(ordered_toggles)))
        for (action, edge) in ordered_toggles[randperm(rng, length(ordered_toggles))[1:count]]
            action == :add_edge ? add_edge!(graph, edge...) : rem_edge!(graph, edge...)
        end
        length(state.boundary) == 4 && !_has_alternating_planar_frame(graph, state.boundary) && continue
        push!(children, _GraphState(
            graph, copy(state.boundary), _graph_state_key(state), :repair_batch, state.rewrite_steps,
        ))
    end
    return children
end

function _boundary_state(occupied, boundary)
    return sum(Int((occupied >> (vertex - 1)) & 1) << (slot - 1)
        for (slot, vertex) in enumerate(boundary))
end

function _is_independent_mask(graph, occupied)
    return all(occupied & (UInt64(1) << (src(edge) - 1)) == 0 ||
        occupied & (UInt64(1) << (dst(edge) - 1)) == 0 for edge in edges(graph))
end

function _search_logical_skeletons(
    target_graph, target_boundary, target_reduced, lattice;
    min_vertices, max_vertices, max_evaluations, beam_width,
    mutations_per_candidate, random_candidates_per_generation,
    exploration_fraction, rng,
)
    boundary_count = length(target_boundary)
    logical_max = min(max_vertices, max(min_vertices, boundary_count + 12))
    beam = [_random_logical_state(rng, boundary_count, min_vertices, logical_max) for _ in 1:beam_width]
    if nv(target_graph) <= logical_max &&
        (boundary_count != 4 || _has_alternating_planar_frame(target_graph, target_boundary))
        push!(beam, _GraphState(deepcopy(target_graph), copy(target_boundary), nothing, :target_seed, 0))
    end
    cache = Dict{String, Tuple}()
    skeletons = _GraphState[]
    skeleton_keys = Set{String}()
    trace = UnweightedSearchRecord[]
    evaluated = 0
    generation = 0
    best_score = (typemax(Int), Inf, typemax(Int), typemax(Int), typemax(Int))
    while evaluated < max_evaluations
        generation += 1
        pool = [_GraphState(state.graph, state.boundary, state.parent_key, :retained, 0) for state in beam]
        for state in beam, _ in 1:mutations_per_candidate
            push!(pool, _mutate_logical_state(rng, state, min_vertices, logical_max))
        end
        for state in beam
            analysis = get(cache, _graph_state_key(state), nothing)
            analysis === nothing && continue
            append!(pool, _tensor_repair_states(
                rng, state, target_reduced, analysis[5], analysis[6], mutations_per_candidate,
            ))
        end
        append!(pool, [_random_logical_state(rng, boundary_count, min_vertices, logical_max) for _ in 1:random_candidates_per_generation])
        ranked = Tuple{_GraphState, Tuple{Int, Float64, Int, Int, Int}, Vector{Float64}}[]
        generation_records = Tuple{_GraphState, String, Tuple{Int, Float64, Int, Int, Int}, Bool, Float64}[]
        seen = Set{String}()
        for state in pool
            key = _graph_state_key(state)
            key in seen && continue
            push!(seen, key)
            if !haskey(cache, key)
                evaluated == max_evaluations && break
                cache[key] = _evaluate_graph_state(state, target_reduced)
                evaluated += 1
                score, _, valid, offset, _, _ = cache[key]
                push!(generation_records, (state, key, score, valid, offset))
            end
            score, signature, valid, offset, _, _ = cache[key]
            best_score = min(best_score, score)
            push!(ranked, (state, score, signature))
            if valid && is_connected(state.graph) && !(key in skeleton_keys)
                push!(skeleton_keys, key)
                push!(skeletons, state)
            end
        end
        sort!(ranked; by=item -> item[2])
        beam = _select_graph_beam(rng, ranked, beam_width, exploration_fraction)
        selected = Set(_graph_state_key(state) for state in beam)
        for (state, key, score, valid, offset) in generation_records
            push!(trace, _search_record(
                generation, :logical, key, lattice, state, score;
                is_solution=false, constant_offset=valid ? offset : nothing,
                embedding_placed=0, selected=key in selected,
            ))
        end
        isempty(generation_records) && break
        isempty(ranked) && break
    end
    return skeletons, evaluated, generation, best_score, trace
end

function _select_graph_beam(rng, ranked, beam_width, exploration_fraction)
    selected_count = min(beam_width, length(ranked))
    selected_count == 0 && return _GraphState[]
    random_count = min(floor(Int, selected_count * exploration_fraction), selected_count - 1)
    elite_count = selected_count - random_count
    selected = eltype(ranked)[]
    connected_quota = min(cld(elite_count, 3), count(item -> is_connected(item[1].graph), ranked))
    for item in ranked
        is_connected(item[1].graph) || continue
        push!(selected, item)
        length(selected) == connected_quota && break
    end
    profiles = Vector{Float64}[]
    selected_keys = Set(_graph_state_key(item[1]) for item in selected)
    for item in ranked
        _graph_state_key(item[1]) in selected_keys && continue
        item[3] in profiles && continue
        push!(selected, item)
        push!(selected_keys, _graph_state_key(item[1]))
        push!(profiles, item[3])
        length(selected) == min(elite_count, cld(selected_count, 2)) && break
    end
    for item in ranked
        key = _graph_state_key(item[1])
        key in selected_keys && continue
        push!(selected, item)
        push!(selected_keys, key)
        length(selected) == elite_count && break
    end
    remaining = [item for item in ranked if !(_graph_state_key(item[1]) in selected_keys)]
    if random_count > 0
        append!(selected, remaining[randperm(rng, length(remaining))[1:random_count]])
    end
    return [item[1] for item in selected]
end

function _search_record(
    generation, stage, key, lattice, state, score;
    patch=nothing, is_solution, constant_offset, embedding_placed, selected,
)
    pin_coordinates = patch === nothing ? nothing : copy(patch.pins)
    return UnweightedSearchRecord(
        generation, stage, key, _lattice_symbol(lattice), _graph_edges(state.graph),
        patch === nothing ? nothing : copy(patch.coordinates), pin_coordinates,
        patch === nothing ? _LatticeCoordinate[] : _patch_ray_directions(lattice, patch),
        copy(state.boundary), state.parent_key, state.action, nv(state.graph), ne(state.graph),
        score[1], score[2], score[3], embedding_placed, state.rewrite_steps,
        is_solution, constant_offset, selected,
    )
end

function _rewrite_and_embed_skeletons(
    target_graph, target_boundary, skeletons, lattice;
    max_vertices, max_evaluations, beam_width, mutations_per_candidate,
    exploration_fraction, max_results, rng, trace,
)
    isempty(skeletons) && return UnweightedGadget[], 0, 0
    beam = unique(_graph_state_key, skeletons)
    seen = Set(_graph_state_key(state) for state in beam)
    gadgets = UnweightedGadget[]
    evaluated = 0
    generation = 0
    while evaluated < max_evaluations && !isempty(beam)
        generation += 1
        ranked = Tuple{_GraphState, Tuple{Int, Int, Int, Int}, Vector{Tuple{Int, Int}}}[]
        for state in beam
            evaluated == max_evaluations && break
            defects = _local_geometry_defects(state.graph, state.boundary, lattice)
            patch, placed, conflicts = defects == 0 ? _embed_induced_graph(
                state.graph, state.boundary, lattice; node_limit=5_000,
            ) : (nothing, 0, Tuple{Int, Int}[])
            evaluated += 1
            score = (defects, -placed, nv(state.graph) - placed, nv(state.graph))
            key = _graph_state_key(state)
            solved = patch !== nothing
            push!(trace, _search_record(
                generation, :rewrite, key, lattice, state, (0, 0.0, defects, nv(state.graph), ne(state.graph));
                patch, is_solution=solved, constant_offset=nothing,
                embedding_placed=placed, selected=true,
            ))
            if solved
                graph, boundary, positions = _materialize_lattice_patch(lattice, patch)
                accepted, lattice_offset = is_gadget_replacement(
                    target_graph, graph, target_boundary, boundary,
                )
                accepted || error("embedded rewrite failed the fixed verifier")
                push!(gadgets, UnweightedGadget(
                    target_graph, graph, boundary, Float64(lattice_offset), _lattice_symbol(lattice),
                    copy(patch.coordinates), positions, _patch_ray_directions(lattice, patch),
                ))
                length(gadgets) == max_results && break
            end
            push!(ranked, (state, score, conflicts))
        end
        length(gadgets) == max_results && break
        sort!(ranked; by=item -> item[2])
        parents = ranked[1:min(beam_width, length(ranked))]
        proposals = _GraphState[]
        for (state, _, conflicts) in parents
            append!(proposals, _rewrite_proposals(
                rng, state, lattice, conflicts, max_vertices, mutations_per_candidate,
            ))
        end
        next_beam = _GraphState[]
        for state in proposals
            key = _graph_state_key(state)
            key in seen && continue
            push!(seen, key)
            push!(next_beam, state)
        end
        shuffle!(rng, next_beam)
        beam = next_beam[1:min(length(next_beam), max(beam_width, floor(Int, beam_width / (1 - exploration_fraction))))]
    end
    sort!(gadgets; by=gadget -> (nv(gadget.replacement_graph), ne(gadget.replacement_graph)))
    return gadgets, evaluated, generation
end

function _rewrite_proposals(rng, state, lattice, conflicts, max_vertices, proposal_count)
    room = max_vertices - nv(state.graph)
    room < 2 && return _GraphState[]
    defects = [
        vertex for vertex in vertices(state.graph)
        if !_ring_is_realizable(state.graph, vertex, lattice)
    ]
    proposals = _GraphState[]
    for _ in 1:proposal_count
        if !isempty(defects) && rand(rng) < 0.75
            vertex = rand(rng, defects)
            if vertex in state.boundary
                graph, boundary = _extend_boundary_pin(state.graph, state.boundary, vertex)
                action = :extend_pin
            else
                neighbors = Graphs.neighbors(state.graph, vertex)
                length(neighbors) >= 2 || continue
                shuffled = shuffle(rng, neighbors)
                split = rand(rng, 1:length(shuffled)-1)
                graph = _split_vertex(state.graph, vertex, shuffled[1:split])
                boundary = copy(state.boundary)
                action = :split_vertex
            end
        else
            edges_to_try = isempty(conflicts) ? _graph_edges(state.graph) :
                unique([conflicts[1:min(8, length(conflicts))]; _graph_edges(state.graph)])
            count = rand(rng, 1:min(6, room ÷ 2, length(edges_to_try)))
            chosen = edges_to_try[randperm(rng, length(edges_to_try))[1:count]]
            graph = _even_subdivide_edges(state.graph, chosen)
            boundary = copy(state.boundary)
            action = :subdivide_edges
        end
        nv(graph) <= max_vertices || continue
        push!(proposals, _GraphState(
            graph, boundary, _graph_state_key(state), action, state.rewrite_steps + 1,
        ))
    end
    return proposals
end

function _split_vertex(graph, vertex, first_neighbors)
    result = deepcopy(graph)
    second_neighbors = setdiff(Graphs.neighbors(graph, vertex), first_neighbors)
    add_vertex!(result)
    bridge = nv(result)
    add_vertex!(result)
    second = nv(result)
    for neighbor in second_neighbors
        rem_edge!(result, vertex, neighbor)
        add_edge!(result, second, neighbor)
    end
    add_edge!(result, vertex, bridge)
    add_edge!(result, bridge, second)
    return result
end

function _extend_boundary_pin(graph, boundary, pin)
    result = deepcopy(graph)
    add_vertex!(result)
    middle = nv(result)
    add_vertex!(result)
    endpoint = nv(result)
    add_edge!(result, pin, middle)
    add_edge!(result, middle, endpoint)
    new_boundary = copy(boundary)
    new_boundary[findfirst(==(pin), boundary)] = endpoint
    return result, new_boundary
end

function _even_subdivide_edges(graph, selected_edges)
    result = deepcopy(graph)
    for (first, second) in selected_edges
        has_edge(result, first, second) || continue
        rem_edge!(result, first, second)
        add_vertex!(result)
        middle_first = nv(result)
        add_vertex!(result)
        middle_second = nv(result)
        add_edge!(result, first, middle_first)
        add_edge!(result, middle_first, middle_second)
        add_edge!(result, middle_second, second)
    end
    return result
end

function _local_geometry_defects(graph, boundary, lattice)
    defects = count(vertex -> !_ring_is_realizable(graph, vertex, lattice), vertices(graph))
    return defects + count(pin -> degree(graph, pin) >= length(_lattice_directions(lattice)), boundary)
end

function _ring_is_realizable(graph, vertex, lattice)
    neighbors = Graphs.neighbors(graph, vertex)
    directions = _lattice_directions(lattice)
    length(neighbors) <= length(directions) || return false
    length(neighbors) <= 1 && return true
    for slots in permutations(eachindex(directions), length(neighbors))
        all(
            has_edge(graph, neighbors[first], neighbors[second]) ==
                (_lattice_distance(lattice, directions[slots[first]], directions[slots[second]]) == 1)
            for first in 1:length(neighbors)-1 for second in first+1:length(neighbors)
        ) && return true
    end
    return false
end

function _embed_induced_graph(graph, boundary, lattice; node_limit=100_000)
    is_connected(graph) || return nothing, 0, Tuple{Int, Int}[]
    directions = _lattice_directions(lattice)
    conflicts = Dict{Tuple{Int, Int}, Int}()
    best_placed = 0
    roots = sort!(collect(vertices(graph)); by=vertex -> (
        degree(graph, vertex),
        count(edge -> src(edge) in Graphs.neighbors(graph, vertex) &&
            dst(edge) in Graphs.neighbors(graph, vertex), edges(graph)),
    ), rev=true)
    for root in roots, first_neighbor in Graphs.neighbors(graph, root)
        placed = Dict(root => (0, 0), first_neighbor => directions[1])
        occupied = Set(values(placed))
        nodes = Ref(0)
        solution = Ref{Union{Nothing, _LatticePatch}}(nothing)
        function visit()
            nodes[] += 1
            nodes[] > node_limit && return false
            best_placed = max(best_placed, length(placed))
            if length(placed) == nv(graph)
                canonical = [placed[vertex] for vertex in vertices(graph)]
                coordinates = _from_canonical.(Ref(lattice), canonical)
                pins = coordinates[boundary]
                normalized = _normalize_lattice_patch(
                    lattice, coordinates, pins, fill(1, length(boundary)),
                )
                ray_choices = length(boundary) == 4 ? Iterators.product(ntuple(_ -> eachindex(directions), 4)...) : (ntuple(_ -> 1, length(boundary)),)
                for rays in ray_choices
                    patch = _LatticePatch(normalized.coordinates, normalized.pins, collect(rays))
                    all(_check_crossing_frame(lattice, patch)) && (solution[] = patch; return true)
                end
                return false
            end
            unplaced = [vertex for vertex in vertices(graph) if !haskey(placed, vertex)]
            vertex = argmax(candidate -> (
                count(neighbor -> haskey(placed, neighbor), Graphs.neighbors(graph, candidate)),
                degree(graph, candidate),
            ), unplaced)
            placed_neighbors = [neighbor for neighbor in Graphs.neighbors(graph, vertex) if haskey(placed, neighbor)]
            isempty(placed_neighbors) && return false
            candidates = Set(
                (placed[placed_neighbors[1]][1] + direction[1], placed[placed_neighbors[1]][2] + direction[2])
                for direction in directions
            )
            for neighbor in placed_neighbors[2:end]
                intersect!(candidates, Set(
                    (placed[neighbor][1] + direction[1], placed[neighbor][2] + direction[2])
                    for direction in directions
                ))
            end
            filter!(candidate -> !(candidate in occupied) && all(
                (_lattice_distance(lattice, candidate, coordinate) == 1) == has_edge(graph, vertex, other)
                for (other, coordinate) in placed
            ), candidates)
            if isempty(candidates)
                for neighbor in placed_neighbors
                    edge = minmax(vertex, neighbor)
                    conflicts[edge] = get(conflicts, edge, 0) + 1
                end
            end
            for candidate in candidates
                placed[vertex] = candidate
                push!(occupied, candidate)
                visit() && return true
                delete!(placed, vertex)
                delete!(occupied, candidate)
            end
            return false
        end
        visit()
        solution[] !== nothing && return solution[], best_placed, Tuple{Int, Int}[]
    end
    ordered_conflicts = sort!(collect(keys(conflicts)); by=edge -> conflicts[edge], rev=true)
    return nothing, best_placed, ordered_conflicts
end

function _materialize_lattice_patch(lattice::LatticeType, patch::_LatticePatch)
    positions = get_physical_positions(lattice, patch.coordinates)
    graph = unit_disk_graph(positions, get_radius(lattice))
    coordinate_index = Dict(coordinate => index for (index, coordinate) in enumerate(patch.coordinates))
    boundary = [coordinate_index[pin] for pin in patch.pins]
    return graph, boundary, positions
end

function _normalize_lattice_patch(
    ::Square,
    coordinates::Vector{_LatticeCoordinate},
    pins::Vector{_LatticeCoordinate},
    rays::Vector{Int},
)
    min_x = minimum(first, coordinates)
    min_y = minimum(last, coordinates)
    translate(point) = (point[1] - min_x, point[2] - min_y)
    return _LatticePatch(sort!(translate.(coordinates)), translate.(pins), copy(rays))
end

function _normalize_lattice_patch(
    ::Triangular,
    coordinates::Vector{_LatticeCoordinate},
    pins::Vector{_LatticeCoordinate},
    rays::Vector{Int},
)
    axial = [_offset_to_axial(point) for point in coordinates]
    pin_axial = [_offset_to_axial(point) for point in pins]
    min_q = minimum(first, axial)
    min_r = minimum(last, axial)
    translate(point) = (point[1] - min_q, point[2] - min_r)
    translated = _axial_to_offset.(translate.(axial))
    translated_pins = _axial_to_offset.(translate.(pin_axial))
    return _LatticePatch(sort!(translated), translated_pins, copy(rays))
end

_offset_to_axial(point::_LatticeCoordinate) = (point[1] - fld(point[2], 2), point[2])
_axial_to_offset(point::_LatticeCoordinate) = (point[1] + fld(point[2], 2), point[2])

_lattice_directions(::Square) = _LatticeCoordinate[
    (-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1),
]
_lattice_directions(::Triangular) = _LatticeCoordinate[
    (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1),
]

function _lattice_step(::Square, point::_LatticeCoordinate, direction::_LatticeCoordinate, distance::Int)
    return (point[1] + distance * direction[1], point[2] + distance * direction[2])
end

function _lattice_step(::Triangular, point::_LatticeCoordinate, direction::_LatticeCoordinate, distance::Int)
    q, r = _offset_to_axial(point)
    return _axial_to_offset((q + distance * direction[1], r + distance * direction[2]))
end

_lattice_symbol(::Square) = :KSG
_lattice_symbol(::Triangular) = :triangular

_patch_ray_directions(lattice::LatticeType, patch::_LatticePatch) =
    _lattice_directions(lattice)[patch.rays]

"""
    check_crossing_frame(lattice, coordinates, pins, pin_rays)

Check the four geometric crossing-frame conditions. `pin_rays[i]` is the
outward lattice direction attached to `pins[i]`. The returned named tuple reports
G1 (strict hull interfaces), G2 (alternating channels), G3 (outward rays), and
G4 (clear pairwise non-adjacent exterior corridors).
"""
function check_crossing_frame(
    lattice::LatticeType,
    coordinates::Vector{_LatticeCoordinate},
    pins::Vector{_LatticeCoordinate},
    pin_rays::Vector{_LatticeCoordinate},
)
    length(pins) == 4 || throw(ArgumentError("a crossing frame requires four ordered pins"))
    length(pin_rays) == 4 || throw(ArgumentError("a crossing frame requires four pin rays"))
    directions = _lattice_directions(lattice)
    ray_indices = [_lattice_direction_index(directions, ray) for ray in pin_rays]
    patch = _normalize_lattice_patch(lattice, coordinates, pins, ray_indices)
    checks = _check_crossing_frame(lattice, patch)
    return (G1=checks[1], G2=checks[2], G3=checks[3], G4=checks[4])
end

function _lattice_direction_index(directions, ray)
    index = findfirst(==(ray), directions)
    index === nothing && throw(ArgumentError("$ray is not a lattice direction"))
    return index
end

function _check_crossing_frame(lattice::LatticeType, patch::_LatticePatch)
    length(patch.pins) == 4 || return (true, true, true, true)
    directions = _patch_ray_directions(lattice, patch)
    interfaces = [
        _lattice_step(lattice, pin, direction, 1)
        for (pin, direction) in zip(patch.pins, directions)
    ]
    occupied_geometry = _geometry_coordinate.(Ref(lattice), patch.coordinates)
    interface_geometry = _geometry_coordinate.(Ref(lattice), interfaces)
    hull = _strict_convex_hull([occupied_geometry; interface_geometry])

    g1 = length(unique(interface_geometry)) == 4 && all(in(hull), interface_geometry)
    g2 = g1 && _interfaces_alternate(hull, interface_geometry)
    g3 = _rays_point_outward(lattice, interfaces, directions)
    g4 = _corridors_are_clear(lattice, patch.coordinates, patch.pins, interfaces, directions)
    return (g1, g2, g3, g4)
end

_canonical_coordinate(::Square, point::_LatticeCoordinate) = point
_canonical_coordinate(::Triangular, point::_LatticeCoordinate) = _offset_to_axial(point)
_from_canonical(::Square, point::_LatticeCoordinate) = point
_from_canonical(::Triangular, point::_LatticeCoordinate) = _axial_to_offset(point)
_lattice_distance(::Square, first, second) = max(abs(first[1] - second[1]), abs(first[2] - second[2]))
_lattice_distance(::Triangular, first, second) = max(abs(first[1] - second[1]), abs(first[2] - second[2]), abs(sum(first) - sum(second)))
_geometry_coordinate(::Square, point::_LatticeCoordinate) = point
function _geometry_coordinate(::Triangular, point::_LatticeCoordinate)
    q, r = _offset_to_axial(point)
    return (2q + r, r)
end

_orientation(a, b, c) =
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])

function _strict_convex_hull(points::Vector{_LatticeCoordinate})
    sorted_points = sort!(unique(points))
    length(sorted_points) <= 2 && return sorted_points
    lower = _LatticeCoordinate[]
    for point in sorted_points
        while length(lower) >= 2 && _orientation(lower[end-1], lower[end], point) <= 0
            pop!(lower)
        end
        push!(lower, point)
    end
    upper = _LatticeCoordinate[]
    for point in Iterators.reverse(sorted_points)
        while length(upper) >= 2 && _orientation(upper[end-1], upper[end], point) <= 0
            pop!(upper)
        end
        push!(upper, point)
    end
    return [lower[1:end-1]; upper[1:end-1]]
end

function _interfaces_alternate(
    hull::Vector{_LatticeCoordinate},
    interfaces::Vector{_LatticeCoordinate},
)
    labels = Dict(point => label for (label, point) in enumerate(interfaces))
    order = [labels[point] for point in hull if haskey(labels, point)]
    length(order) == 4 || return false
    return all(isodd(order[index]) != isodd(order[mod1(index + 1, 4)]) for index in 1:4)
end

function _rays_point_outward(
    lattice::LatticeType,
    interfaces::Vector{_LatticeCoordinate},
    directions::Vector{_LatticeCoordinate},
)
    canonical_interfaces = _canonical_coordinate.(Ref(lattice), interfaces)
    sum_q = sum(first, canonical_interfaces)
    sum_r = sum(last, canonical_interfaces)
    for (interface, direction) in zip(canonical_interfaces, directions)
        out_q = 4interface[1] - sum_q
        out_r = 4interface[2] - sum_r
        if lattice isa Square
            out_q * direction[1] + out_r * direction[2] > 0 || return false
        else
            out_x = 2out_q + out_r
            direction_x = 2direction[1] + direction[2]
            out_x * direction_x + 3out_r * direction[2] > 0 || return false
        end
    end
    return true
end

function _corridors_are_clear(
    lattice::LatticeType,
    coordinates::Vector{_LatticeCoordinate},
    pins::Vector{_LatticeCoordinate},
    interfaces::Vector{_LatticeCoordinate},
    directions::Vector{_LatticeCoordinate},
)
    occupied = _canonical_coordinate.(Ref(lattice), coordinates)
    canonical_pins = _canonical_coordinate.(Ref(lattice), pins)
    starts = _canonical_coordinate.(Ref(lattice), interfaces)
    adjacency_offsets = [_LatticeCoordinate[(0, 0)]; _lattice_directions(lattice)]

    for index in eachindex(starts), site in occupied
        site == canonical_pins[index] && continue
        for offset in adjacency_offsets
            _point_on_ray((site[1] + offset[1], site[2] + offset[2]), starts[index], directions[index]) &&
                return false
        end
    end
    for first in 1:3, second in first+1:4, offset in adjacency_offsets
        _rays_touch(
            starts[first], directions[first], starts[second], directions[second], offset,
        ) && return false
    end
    return true
end

function _point_on_ray(point, start, direction)
    displacement = (point[1] - start[1], point[2] - start[2])
    multiple = _direction_multiple(displacement, direction)
    return multiple !== nothing && multiple >= 0
end

function _direction_multiple(displacement, direction)
    if direction[1] != 0
        rem(displacement[1], direction[1]) == 0 || return nothing
        multiple = div(displacement[1], direction[1])
    else
        direction[2] != 0 || error("zero ray direction")
        rem(displacement[2], direction[2]) == 0 || return nothing
        multiple = div(displacement[2], direction[2])
    end
    displacement == (multiple * direction[1], multiple * direction[2]) || return nothing
    return multiple
end

function _rays_touch(start1, direction1, start2, direction2, offset)
    right = (start2[1] + offset[1] - start1[1], start2[2] + offset[2] - start1[2])
    determinant = direction1[2] * direction2[1] - direction1[1] * direction2[2]
    if determinant != 0
        first_numerator = right[2] * direction2[1] - right[1] * direction2[2]
        second_numerator = direction1[1] * right[2] - direction1[2] * right[1]
        rem(first_numerator, determinant) == 0 || return false
        rem(second_numerator, determinant) == 0 || return false
        return div(first_numerator, determinant) >= 0 && div(second_numerator, determinant) >= 0
    end

    if direction1 == direction2
        return _direction_multiple(right, direction1) !== nothing
    end
    multiple = _direction_multiple(right, direction1)
    return multiple !== nothing && multiple >= 0
end

# ============================================================================
# Alpha Tensor Functions
# ============================================================================

"""
    calculate_alpha_tensor(graph, boundary_vertices) -> Array{<:Tropical}

Compute the alpha tensor α(R) via tropical tensor network contraction.
Element `α(R)_s` is the maximum independent set size with boundary fixed to `s`,
or `-Inf` if `s` violates the independent set constraint.
"""
function calculate_alpha_tensor(graph::SimpleGraph{Int}, boundary_vertices::Vector{Int})
    return solve(GenericTensorNetwork(IndependentSet(graph), openvertices=boundary_vertices), SizeMax())
end

"""
    calculate_reduced_alpha_tensor(graph, boundary_vertices) -> Array{Float64}

Compute the reduced alpha tensor α̃(R) by applying `mis_compactify!` to `α(R)`.
Dominated boundary configurations (where a subset achieves equal or better MIS size)
are set to `-Inf`.
"""
function calculate_reduced_alpha_tensor(graph::SimpleGraph{Int}, boundary_vertices::Vector{Int})
    return Float64.(content.(mis_compactify!(calculate_alpha_tensor(graph, boundary_vertices))))
end

# ============================================================================
# Tensor Utility Functions
# ============================================================================

"""
    inf_mask(tensor) -> BigInt

Bitmask encoding `-Inf` positions in `tensor` (LSB = first linear index).
"""
function inf_mask(tensor::AbstractArray)
    mask = BigInt(0)
    for (i, v) in enumerate(tensor)
        v == -Inf && (mask |= BigInt(1) << (i - 1))
    end
    return mask
end

"""
    pins_prefilter(g, pins)

Return `true` when every connected component of `g` contains at least one pin.
"""
function pins_prefilter(g::SimpleGraph{Int}, pins::AbstractVector{<:Integer})
    isempty(pins) && return false
    n = Graphs.nv(g)
    unique_pins = unique(pins)
    length(unique_pins) == length(pins) || error("pins must be unique")
    all(p -> 1 <= p <= n, unique_pins) || error("pins must be valid vertex indices for a graph with $n vertices")
    pinset = Set(unique_pins)
    return all(component -> any(in(pinset), component), Graphs.connected_components(g))
end

"""
    is_diff_by_constant(t1, t2) -> (Bool, Real)

Return `(true, c)` if `t1 - t2 == c` at all finite entries and both tensors share
the same `-Inf` pattern; `(false, 0)` otherwise.
"""
function is_diff_by_constant(t1::AbstractArray{T}, t2::AbstractArray{T}) where T <: Real
    size(t1) == size(t2) || throw(DimensionMismatch("input tensors must have the same size, got $(size(t1)) and $(size(t2))"))
    any(isinf(a) ⊻ isinf(b) for (a, b) in zip(t1, t2)) && return false, zero(T)
    c = nothing
    for (a, b) in zip(t1, t2)
        isfinite(a) || continue
        d = a - b
        c === nothing ? (c = d) : (d == c || return false, zero(T))
    end
    c === nothing && throw(ArgumentError("input tensors must contain at least one finite entry"))
    return true, c
end

"""
    is_gadget_replacement(g1, g2, open_vertices1, open_vertices2) -> (Bool, Real)

Check whether `g2` is a valid MIS replacement for `g1`, i.e., their reduced alpha
tensors differ only by a constant. Returns `(is_valid, α̃(g2) - α̃(g1))`.
"""
function is_gadget_replacement(g1::SimpleGraph{Int}, g2::SimpleGraph{Int},
                                open_vertices1::Vector{Int}, open_vertices2::Vector{Int})
    t1 = calculate_reduced_alpha_tensor(g1, open_vertices1)
    t2 = calculate_reduced_alpha_tensor(g2, open_vertices2)
    return is_diff_by_constant(t2, t1)
end
