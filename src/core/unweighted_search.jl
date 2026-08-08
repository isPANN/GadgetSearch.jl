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
    key::String
    lattice::Symbol
    lattice_coordinates::Vector{_LatticeCoordinate}
    pin_coordinates::Vector{_LatticeCoordinate}
    pin_rays::Vector{_LatticeCoordinate}
    graph6::String
    boundary_vertices::Vector{Int}
    parent_key::Union{Nothing, String}
    action::Symbol
    vertices::Int
    edges::Int
    mask_mismatches::Int
    offset_spread::Float64
    frame_violations::Int
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
only logical acceptance criterion. Four-pin searches additionally require the
complete G1-G4 crossing-frame geometry.
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
                frame = _check_crossing_frame(lattice, proposal.patch)
                frame_violations = count(!, frame)
                score = _unweighted_tensor_distance(
                    candidate_reduced, target_reduced, graph, frame_violations,
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
            if valid && score[3] == 0 && !(key in gadget_keys)
                push!(gadget_keys, key)
                push!(gadgets, UnweightedGadget(
                    target_graph,
                    graph,
                    boundary,
                    constant_offset,
                    _lattice_symbol(lattice),
                    copy(proposal.patch.coordinates),
                    positions,
                    _patch_ray_directions(lattice, proposal.patch),
                ))
                sort!(gadgets; by=gadget -> (
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
                _patch_ray_directions(lattice, patch),
                graph_to_g6(item.graph),
                copy(item.boundary),
                item.proposal.parent_key,
                item.proposal.action,
                nv(item.graph),
                ne(item.graph),
                item.score[1],
                item.score[2],
                item.score[3],
                item.valid && item.score[3] == 0,
                item.valid && item.score[3] == 0 ? item.constant_offset : nothing,
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
    rays = rand(rng, eachindex(_lattice_directions(lattice)), boundary_count)
    return _normalize_lattice_patch(lattice, coordinates, pins, rays)
end

function _mutate_lattice_patch(
    rng::AbstractRNG,
    lattice::LatticeType,
    patch::_LatticePatch,
    min_vertices::Int,
    max_vertices::Int,
)
    operations = Symbol[:move_pin, :swap_pins, :change_ray]
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
        normalized = _normalize_lattice_patch(lattice, mutated.coordinates, mutated.pins, mutated.rays)
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
        return _LatticePatch(coordinates, copy(patch.pins), copy(patch.rays))
    elseif action == :remove_site
        removed = rand(rng, removable)
        coordinates = setdiff(patch.coordinates, [removed])
        return _connected_lattice_patch(lattice, coordinates) ? _LatticePatch(coordinates, copy(patch.pins), copy(patch.rays)) : nothing
    elseif action == :relocate_site
        removed = rand(rng, removable)
        coordinates = setdiff(patch.coordinates, [removed])
        isempty(coordinates) && return nothing
        moved = rand(rng, _lattice_frontier(lattice, Set(coordinates)))
        relocated = [coordinates; moved]
        return _connected_lattice_patch(lattice, relocated) ? _LatticePatch(relocated, copy(patch.pins), copy(patch.rays)) : nothing
    elseif action == :extend_arm
        base = rand(rng, patch.coordinates)
        direction = rand(rng, _lattice_directions(lattice))
        first = _lattice_step(lattice, base, direction, 1)
        second = _lattice_step(lattice, base, direction, 2)
        (first in occupied || second in occupied) && return nothing
        return _LatticePatch([patch.coordinates; first; second], copy(patch.pins), copy(patch.rays))
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
        return _connected_lattice_patch(lattice, coordinates) ? _LatticePatch(coordinates, copy(patch.pins), copy(patch.rays)) : nothing
    elseif action == :move_pin
        choices = setdiff(patch.coordinates, patch.pins)
        isempty(choices) && return nothing
        pins = copy(patch.pins)
        pins[rand(rng, eachindex(pins))] = rand(rng, choices)
        return _LatticePatch(copy(patch.coordinates), pins, copy(patch.rays))
    elseif action == :swap_pins
        length(patch.pins) < 2 && return nothing
        first, second = randperm(rng, length(patch.pins))[1:2]
        pins = copy(patch.pins)
        pins[first], pins[second] = pins[second], pins[first]
        rays = copy(patch.rays)
        rays[first], rays[second] = rays[second], rays[first]
        return _LatticePatch(copy(patch.coordinates), pins, rays)
    else
        rays = copy(patch.rays)
        slot = rand(rng, eachindex(rays))
        choices = setdiff(eachindex(_lattice_directions(lattice)), [rays[slot]])
        rays[slot] = rand(rng, choices)
        return _LatticePatch(copy(patch.coordinates), copy(patch.pins), rays)
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
    return string(_lattice_symbol(lattice), ':', coordinates, '|', pins, '|', join(patch.rays, ','))
end

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

function _unweighted_tensor_distance(
    candidate::AbstractArray,
    target::AbstractArray,
    graph::SimpleGraph,
    frame_violations::Int,
)
    mask_mismatches = count(isinf(a) != isinf(b) for (a, b) in zip(candidate, target))
    differences = [a - b for (a, b) in zip(candidate, target) if isfinite(a) && isfinite(b)]
    offset_spread = Float64(maximum(differences) - minimum(differences))
    return mask_mismatches, offset_spread, frame_violations, nv(graph), ne(graph)
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
                pin_rays=record.pin_rays,
                graph6=record.graph6,
                boundary_vertices=record.boundary_vertices,
                parent_key=record.parent_key,
                action=String(record.action),
                vertices=record.vertices,
                edges=record.edges,
                mask_mismatches=record.mask_mismatches,
                offset_spread=record.offset_spread,
                frame_violations=record.frame_violations,
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
