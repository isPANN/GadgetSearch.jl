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
    port_crossing_penalty::Int
end

"""One evaluated transition of the dynamically constructed lattice patch."""
struct UnweightedSearchRecord
    generation::Int
    key::String
    lattice::Symbol
    lattice_coordinates::Vector{_LatticeCoordinate}
    pin_coordinates::Vector{_LatticeCoordinate}
    graph6::String
    boundary_vertices::Vector{Int}
    parent_key::Union{Nothing, String}
    action::Symbol
    vertices::Int
    edges::Int
    mask_mismatches::Int
    offset_spread::Float64
    port_crossing_penalty::Int
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
end

struct _UnweightedProposal
    patch::_LatticePatch
    parent_key::Union{Nothing, String}
    action::Symbol
end

struct _UnweightedEvaluation
    proposal::_UnweightedProposal
    key::String
    graph::SimpleGraph{Int}
    boundary::Vector{Int}
    score::Tuple{Int, Float64, Int, Int, Int}
    valid::Bool
    constant_offset::Float64
end

"""
    search_unweighted_gadgets(target_graph, target_boundary, lattice; kwargs...)

Dynamically grow and reshape an induced patch of `Square()` (KSG) or
`Triangular()`. There is no fixed canvas and no abstract-graph stage: every
evaluated state is an explicit lattice coordinate set, and its edges are derived
from the selected lattice sites. Two-site arm extension and crowded-site split
are proposal moves, while the existing reduced-alpha-tensor verifier remains the
only acceptance criterion. The full evaluation budget is used so later accepted
states can improve the soft four-port crossing preference.
"""
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

    initial_max = min(max_vertices, min_vertices + 2)
    beam = [_random_lattice_patch(rng, lattice, boundary_count, min_vertices, initial_max) for _ in 1:beam_width]
    cache = Dict{String, Tuple{Tuple{Int, Float64, Int, Int, Int}, Float64, Bool}}()
    gadgets = UnweightedGadget[]
    gadget_keys = Set{String}()
    evaluated = 0
    generations = 0
    best_score = (typemax(Int), Inf, typemax(Int), typemax(Int), typemax(Int))
    trace = UnweightedSearchRecord[]

    while evaluated < max_evaluations
        generations += 1
        evaluated_before_generation = evaluated
        pool = [_UnweightedProposal(patch, nothing, generations == 1 ? :seed : :retained) for patch in beam]
        for patch in beam, _ in 1:mutations_per_candidate
            mutated, action = _mutate_lattice_patch(rng, lattice, patch, min_vertices, max_vertices)
            push!(pool, _UnweightedProposal(mutated, _lattice_patch_key(lattice, patch), action))
        end
        for _ in 1:random_candidates_per_generation
            restart = _random_lattice_patch(rng, lattice, boundary_count, min_vertices, initial_max)
            push!(pool, _UnweightedProposal(restart, nothing, :restart))
        end

        ranked = Tuple{_LatticePatch, Tuple{Int, Float64, Int, Int, Int}}[]
        ranked_keys = Set{String}()
        generation_evaluations = _UnweightedEvaluation[]
        for proposal in pool
            key = _lattice_patch_key(lattice, proposal.patch)
            key in ranked_keys && continue
            push!(ranked_keys, key)
            graph, boundary, positions = _materialize_lattice_patch(lattice, proposal.patch)
            if !haskey(cache, key)
                evaluated == max_evaluations && break
                candidate_reduced = vec(calculate_reduced_alpha_tensor(graph, boundary))
                valid, constant_offset = is_diff_by_constant(candidate_reduced, target_reduced)
                port_crossing_penalty = _port_crossing_penalty(positions, boundary)
                score = _unweighted_tensor_distance(
                    candidate_reduced, target_reduced, graph, port_crossing_penalty,
                )
                cache[key] = (score, Float64(constant_offset), valid)
                push!(generation_evaluations, _UnweightedEvaluation(
                    proposal, key, graph, boundary, score, valid, Float64(constant_offset),
                ))
                evaluated += 1
            end

            score, constant_offset, valid = cache[key]
            best_score = min(best_score, score)
            push!(ranked, (proposal.patch, score))
            if valid && !(key in gadget_keys)
                push!(gadget_keys, key)
                push!(gadgets, UnweightedGadget(
                    target_graph,
                    graph,
                    boundary,
                    constant_offset,
                    _lattice_symbol(lattice),
                    copy(proposal.patch.coordinates),
                    positions,
                    score[3],
                ))
                sort!(gadgets; by=gadget -> (
                    gadget.port_crossing_penalty,
                    nv(gadget.replacement_graph),
                    ne(gadget.replacement_graph),
                ))
                resize!(gadgets, min(length(gadgets), max_results))
            end
        end

        sort!(ranked; by=last)
        beam = _select_unweighted_beam(rng, ranked, beam_width, exploration_fraction)
        selected_keys = Set(_lattice_patch_key(lattice, patch) for patch in beam)
        for item in generation_evaluations
            patch = item.proposal.patch
            push!(trace, UnweightedSearchRecord(
                generations,
                item.key,
                _lattice_symbol(lattice),
                copy(patch.coordinates),
                copy(patch.pins),
                graph_to_g6(item.graph),
                copy(item.boundary),
                item.proposal.parent_key,
                item.proposal.action,
                nv(item.graph),
                ne(item.graph),
                item.score[1],
                item.score[2],
                item.score[3],
                item.valid,
                item.valid ? item.constant_offset : nothing,
                item.key in selected_keys,
            ))
        end
        evaluated == evaluated_before_generation && break
    end

    termination_reason = if !isempty(gadgets)
        :solution
    elseif evaluated == max_evaluations
        :budget
    else
        :search_space_exhausted
    end

    return UnweightedSearchResult(
        target_graph,
        copy(target_boundary),
        _lattice_symbol(lattice),
        gadgets,
        evaluated,
        generations,
        best_score[1],
        best_score[2],
        termination_reason,
        trace,
    )
end

function _select_unweighted_beam(
    rng::AbstractRNG,
    ranked::Vector{Tuple{_LatticePatch, Tuple{Int, Float64, Int, Int, Int}}},
    beam_width::Int,
    exploration_fraction::Float64,
)
    selected_count = min(beam_width, length(ranked))
    selected_count == 0 && return _LatticePatch[]
    exploration_count = min(floor(Int, selected_count * exploration_fraction), selected_count - 1)
    elite_count = selected_count - exploration_count
    selected = collect(Iterators.take(ranked, elite_count))
    if exploration_count > 0
        remaining = @view ranked[elite_count+1:end]
        chosen = randperm(rng, length(remaining))[1:exploration_count]
        append!(selected, remaining[chosen])
    end
    return [patch for (patch, _) in selected]
end

function _random_lattice_patch(
    rng::AbstractRNG,
    lattice::LatticeType,
    boundary_count::Int,
    min_vertices::Int,
    max_vertices::Int,
)
    target_size = rand(rng, min_vertices:max_vertices)
    occupied = Set{_LatticeCoordinate}([(0, 0)])
    while length(occupied) < target_size
        push!(occupied, rand(rng, _lattice_frontier(lattice, occupied)))
    end
    coordinates = sort!(collect(occupied))
    pins = coordinates[randperm(rng, length(coordinates))[1:boundary_count]]
    return _normalize_lattice_patch(lattice, coordinates, pins)
end

function _mutate_lattice_patch(
    rng::AbstractRNG,
    lattice::LatticeType,
    patch::_LatticePatch,
    min_vertices::Int,
    max_vertices::Int,
)
    operations = Symbol[:move_pin, :swap_pins]
    length(patch.coordinates) < max_vertices && push!(operations, :add_site)
    length(patch.coordinates) + 2 <= max_vertices && push!(operations, :extend_arm)
    removable = setdiff(patch.coordinates, patch.pins)
    if !isempty(removable)
        length(patch.coordinates) > min_vertices && push!(operations, :remove_site)
        push!(operations, :relocate_site)
        length(patch.coordinates) < max_vertices && push!(operations, :split_crowded_site)
    end

    for _ in 1:16
        action = rand(rng, operations)
        mutated = _apply_lattice_action(rng, lattice, patch, action, removable)
        mutated === nothing && continue
        normalized = _normalize_lattice_patch(lattice, mutated.coordinates, mutated.pins)
        _lattice_patch_key(lattice, normalized) != _lattice_patch_key(lattice, patch) &&
            return normalized, action
    end
    return patch, :rejected_edit
end

function _apply_lattice_action(
    rng::AbstractRNG,
    lattice::LatticeType,
    patch::_LatticePatch,
    action::Symbol,
    removable::Vector{_LatticeCoordinate},
)
    occupied = Set(patch.coordinates)
    if action == :add_site
        coordinates = [patch.coordinates; rand(rng, _lattice_frontier(lattice, occupied))]
        return _LatticePatch(coordinates, copy(patch.pins))
    elseif action == :remove_site
        removed = rand(rng, removable)
        coordinates = setdiff(patch.coordinates, [removed])
        return _connected_lattice_patch(lattice, coordinates) ? _LatticePatch(coordinates, copy(patch.pins)) : nothing
    elseif action == :relocate_site
        removed = rand(rng, removable)
        coordinates = setdiff(patch.coordinates, [removed])
        isempty(coordinates) && return nothing
        moved = rand(rng, _lattice_frontier(lattice, Set(coordinates)))
        relocated = [coordinates; moved]
        return _connected_lattice_patch(lattice, relocated) ? _LatticePatch(relocated, copy(patch.pins)) : nothing
    elseif action == :extend_arm
        base = rand(rng, patch.coordinates)
        direction = rand(rng, _lattice_directions(lattice))
        first = _lattice_step(lattice, base, direction, 1)
        second = _lattice_step(lattice, base, direction, 2)
        (first in occupied || second in occupied) && return nothing
        return _LatticePatch([patch.coordinates; first; second], copy(patch.pins))
    elseif action == :split_crowded_site
        graph, _, _ = _materialize_lattice_patch(lattice, patch)
        coordinate_index = Dict(coordinate => index for (index, coordinate) in enumerate(patch.coordinates))
        crowded = [coordinate for coordinate in removable if degree(graph, coordinate_index[coordinate]) >= 3]
        isempty(crowded) && return nothing
        removed = rand(rng, crowded)
        empty_neighbors = setdiff(_lattice_neighbors(lattice, removed), patch.coordinates)
        length(empty_neighbors) < 2 && return nothing
        chosen = empty_neighbors[randperm(rng, length(empty_neighbors))[1:2]]
        coordinates = [setdiff(patch.coordinates, [removed]); chosen]
        return _connected_lattice_patch(lattice, coordinates) ? _LatticePatch(coordinates, copy(patch.pins)) : nothing
    elseif action == :move_pin
        choices = setdiff(patch.coordinates, patch.pins)
        isempty(choices) && return nothing
        pins = copy(patch.pins)
        pins[rand(rng, eachindex(pins))] = rand(rng, choices)
        return _LatticePatch(copy(patch.coordinates), pins)
    else
        length(patch.pins) < 2 && return nothing
        first, second = randperm(rng, length(patch.pins))[1:2]
        pins = copy(patch.pins)
        pins[first], pins[second] = pins[second], pins[first]
        return _LatticePatch(copy(patch.coordinates), pins)
    end
end

function _materialize_lattice_patch(lattice::LatticeType, patch::_LatticePatch)
    positions = get_physical_positions(lattice, patch.coordinates)
    graph = unit_disk_graph(positions, get_radius(lattice))
    coordinate_index = Dict(coordinate => index for (index, coordinate) in enumerate(patch.coordinates))
    boundary = [coordinate_index[pin] for pin in patch.pins]
    return graph, boundary, positions
end

function _connected_lattice_patch(lattice::LatticeType, coordinates::Vector{_LatticeCoordinate})
    positions = get_physical_positions(lattice, sort(coordinates))
    return is_connected(unit_disk_graph(positions, get_radius(lattice)))
end

function _normalize_lattice_patch(
    ::Square,
    coordinates::Vector{_LatticeCoordinate},
    pins::Vector{_LatticeCoordinate},
)
    min_x = minimum(first, coordinates)
    min_y = minimum(last, coordinates)
    translate(point) = (point[1] - min_x, point[2] - min_y)
    return _LatticePatch(sort!(translate.(coordinates)), translate.(pins))
end

function _normalize_lattice_patch(
    ::Triangular,
    coordinates::Vector{_LatticeCoordinate},
    pins::Vector{_LatticeCoordinate},
)
    axial = [_offset_to_axial(point) for point in coordinates]
    pin_axial = [_offset_to_axial(point) for point in pins]
    min_q = minimum(first, axial)
    min_r = minimum(last, axial)
    translate(point) = (point[1] - min_q, point[2] - min_r)
    translated = _axial_to_offset.(translate.(axial))
    translated_pins = _axial_to_offset.(translate.(pin_axial))
    return _LatticePatch(sort!(translated), translated_pins)
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

_lattice_neighbors(lattice::LatticeType, point::_LatticeCoordinate) =
    [_lattice_step(lattice, point, direction, 1) for direction in _lattice_directions(lattice)]

function _lattice_frontier(lattice::LatticeType, occupied::Set{_LatticeCoordinate})
    frontier = Set{_LatticeCoordinate}()
    for point in occupied, neighbor in _lattice_neighbors(lattice, point)
        neighbor in occupied || push!(frontier, neighbor)
    end
    return collect(frontier)
end

_lattice_symbol(::Square) = :KSG
_lattice_symbol(::Triangular) = :triangular

function _lattice_patch_key(lattice::LatticeType, patch::_LatticePatch)
    coordinates = join(("$(x),$(y)" for (x, y) in patch.coordinates), ';')
    pins = join(("$(x),$(y)" for (x, y) in patch.pins), ';')
    return string(_lattice_symbol(lattice), ':', coordinates, '|', pins)
end

function _unweighted_tensor_distance(
    candidate::AbstractArray,
    target::AbstractArray,
    graph::SimpleGraph,
    port_crossing_penalty::Int,
)
    mask_mismatches = count(isinf(a) != isinf(b) for (a, b) in zip(candidate, target))
    differences = [a - b for (a, b) in zip(candidate, target) if isfinite(a) && isfinite(b)]
    offset_spread = Float64(maximum(differences) - minimum(differences))
    return mask_mismatches, offset_spread, port_crossing_penalty, nv(graph), ne(graph)
end

function _port_crossing_penalty(
    positions::Vector{Tuple{Float64, Float64}},
    boundary::Vector{Int},
)
    length(boundary) == 4 || return 0
    first_start, second_start, first_end, second_end = positions[boundary]
    orientation(a, b, c) =
        (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    first_side = orientation(first_start, first_end, second_start) *
        orientation(first_start, first_end, second_end)
    second_side = orientation(second_start, second_end, first_start) *
        orientation(second_start, second_end, first_end)
    return first_side < 0 && second_side < 0 ? 0 : 1
end

"""Write the self-contained lattice search trajectory as JSON Lines."""
function save_unweighted_trace(path::AbstractString, result::UnweightedSearchResult)
    target_graph6 = graph_to_g6(result.target_graph)
    open(path, "w") do io
        for record in result.trace
            JSON3.write(io, (
                target_graph6=target_graph6,
                target_boundary=result.target_boundary,
                lattice=String(record.lattice),
                generation=record.generation,
                key=record.key,
                lattice_coordinates=record.lattice_coordinates,
                pin_coordinates=record.pin_coordinates,
                graph6=record.graph6,
                boundary_vertices=record.boundary_vertices,
                parent_key=record.parent_key,
                action=String(record.action),
                vertices=record.vertices,
                edges=record.edges,
                mask_mismatches=record.mask_mismatches,
                offset_spread=record.offset_spread,
                port_crossing_penalty=record.port_crossing_penalty,
                is_solution=record.is_solution,
                constant_offset=record.constant_offset,
                selected=record.selected,
            ))
            write(io, '\n')
        end
    end
    return String(path)
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
