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

"""Outcome of a bounded direct SAT search on one concrete lattice."""
struct UnweightedSearchResult
    target_graph::SimpleGraph{Int}
    target_boundary::Vector{Int}
    lattice::Symbol
    gadgets::Vector{UnweightedGadget}
    evaluated::Int
    termination_reason::Symbol
end

"""One verifier-certified atom-reducing rewrite."""
struct UnweightedRewriteStep
    rule::Symbol
    before::UnweightedGadget
    after::UnweightedGadget
end

"""Result of rewriting and re-synthesizing one realized lattice gadget."""
struct UnweightedOptimizationResult
    gadget::UnweightedGadget
    steps::Vector{UnweightedRewriteStep}
    sat_evaluations::Int
    termination_reason::Symbol
end

struct _LatticePatch
    coordinates::Vector{_LatticeCoordinate}
    pins::Vector{_LatticeCoordinate}
    rays::Vector{Int}
end

"""
Search a four-pin unweighted gadget in a finite window of `lattice`.
Every returned gadget satisfies the reduced-alpha target up to a constant and
the four-direction crossing geometry.
"""
function search_unweighted_gadgets(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    lattice::LatticeType=Triangular();
    min_vertices::Int=length(target_boundary) + 1,
    max_vertices::Int=min_vertices + 8,
    max_evaluations::Int=2_000,
    max_frame_evaluations::Int=1_000_000,
    max_results::Int=1,
    window_side::Int=4,
    checkpoint_path::Union{Nothing, String}=nothing,
    checkpoint_interval::Int=10_000,
)
    boundary_count = length(target_boundary)
    boundary_count == 4 ||
        throw(ArgumentError("unweighted lattice search requires four boundary vertices"))
    min_vertices >= boundary_count ||
        throw(ArgumentError("min_vertices must be at least the number of boundary vertices"))
    max_vertices >= min_vertices ||
        throw(ArgumentError("max_vertices must be at least min_vertices"))
    max_evaluations > 0 ||
        throw(ArgumentError("max_evaluations must be positive"))
    max_frame_evaluations > 0 ||
        throw(ArgumentError("max_frame_evaluations must be positive"))
    max_results > 0 || throw(ArgumentError("max_results must be positive"))
    window_side >= 2 || throw(ArgumentError("window_side must be at least 2"))
    checkpoint_interval > 0 ||
        throw(ArgumentError("checkpoint_interval must be positive"))

    target_reduced = vec(calculate_reduced_alpha_tensor(target_graph, target_boundary))
    all(isinf, target_reduced) &&
        error("target graph has an entirely -Inf reduced alpha tensor")
    gadgets, evaluated, reason = _search_crossing_sat(
        target_graph, target_boundary, target_reduced, lattice;
        min_vertices, max_vertices, max_evaluations, max_frame_evaluations,
        max_results, window_side, checkpoint_path, checkpoint_interval,
    )
    return UnweightedSearchResult(
        target_graph, copy(target_boundary), _lattice_symbol(lattice),
        gadgets, evaluated, reason,
    )
end

"""
Search one joint occupancy-and-frame SAT instance in a fixed lattice window.

Unlike `search_unweighted_gadgets`, this formulation chooses the occupied sites,
four pins, and four outward rays in one CNF. `atom_count` and `offset` identify
the exact instance to solve. The unchanged gadget verifier checks every result.
"""
function search_unweighted_gadget_joint(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    lattice::LatticeType=Triangular();
    window_shape::Tuple{Int, Int},
    atom_count::Int,
    offset::Int,
    seconds::Int=600,
    kissat_executable::Union{Nothing, String}=nothing,
    seed::Int=1,
    canonical_shift::Tuple{Int, Int}=(0, 0),
    first_ray::Int=1,
)
    length(target_boundary) == 4 ||
        throw(ArgumentError("joint unweighted search requires four boundary vertices"))
    all(>=(2), window_shape) ||
        throw(ArgumentError("window dimensions must be at least 2"))
    length(target_boundary) <= atom_count <= prod(window_shape) ||
        throw(ArgumentError("atom_count must fit in the window"))
    seconds > 0 || throw(ArgumentError("seconds must be positive"))
    seed > 0 || throw(ArgumentError("seed must be positive"))
    directions = _lattice_directions(lattice)
    first_ray in eachindex(directions) ||
        throw(ArgumentError("first_ray is not a lattice direction index"))

    target_reduced = vec(calculate_reduced_alpha_tensor(
        target_graph, target_boundary,
    ))
    all(value -> isinf(value) || isinteger(value), target_reduced) ||
        error("the reduced alpha tensor must contain integers or -Inf")
    analysis = _solve_joint_crossing_sat(
        target_reduced, lattice, window_shape, atom_count, offset;
        seconds, kissat_executable, initial_seed=seed, canonical_shift,
        first_direction_index=first_ray,
    )
    analysis === nothing && return nothing
    verified, verified_offset = is_gadget_replacement(
        target_graph, analysis.graph, target_boundary, analysis.boundary,
    )
    verified || error("joint SAT result failed the gadget verifier")
    verified_offset == analysis.offset ||
        error("joint SAT result and gadget verifier disagree on the offset")
    return _unweighted_gadget(target_graph, lattice, analysis)
end

"""
    optimize_unweighted_gadget(gadget, target_boundary; kwargs...)

Reduce a verifier-accepted four-pin lattice gadget with certified rewrite rules.
The optimizer contracts even boundary tails, contracts opposite leaf pins in
pairs, and uses fixed-frame SAT to re-synthesize the interior after a one-step
frame rewrite. Direct rewrites are explored to a closure before re-synthesis,
so the optimizer does not commit to the first smaller direct result. It never
proposes arbitrary vertex deletion.
"""
function optimize_unweighted_gadget(
    gadget::UnweightedGadget,
    target_boundary::Vector{Int};
    min_vertices::Int=length(target_boundary),
    max_sat_evaluations::Int=256,
    host_radius::Int=1,
)
    length(target_boundary) == 4 ||
        throw(ArgumentError("unweighted rewrite optimization requires four target boundary vertices"))
    min_vertices >= length(target_boundary) ||
        throw(ArgumentError("min_vertices must include all boundary vertices"))
    min_vertices <= nv(gadget.replacement_graph) ||
        throw(ArgumentError("min_vertices exceeds the current gadget size"))
    max_sat_evaluations >= 0 ||
        throw(ArgumentError("max_sat_evaluations must be nonnegative"))
    host_radius >= 0 || throw(ArgumentError("host_radius must be nonnegative"))

    lattice = _gadget_lattice(gadget)
    target_reduced = vec(calculate_reduced_alpha_tensor(
        gadget.pattern_graph, target_boundary,
    ))
    valid, offset = is_gadget_replacement(
        gadget.pattern_graph, gadget.replacement_graph, target_boundary,
        gadget.boundary_vertices,
    )
    valid || throw(ArgumentError("the input gadget does not replace the requested target boundary"))
    offset == gadget.constant_offset ||
        throw(ArgumentError("the input gadget stores the wrong constant offset"))
    current = gadget
    steps = UnweightedRewriteStep[]
    sat_evaluations = 0
    while true
        nv(current.replacement_graph) == min_vertices &&
            return UnweightedOptimizationResult(current, steps, sat_evaluations, :minimum_vertices)

        current_key = _rewrite_state_key(current)
        direct_states = [(gadget=current, steps=steps)]
        seen_states = Set{Tuple}([current_key])
        for state in direct_states
            for candidate in _direct_unweighted_rewrites(
                state.gadget, target_boundary, target_reduced, lattice, min_vertices,
            )
                candidate_steps = [
                    state.steps;
                    UnweightedRewriteStep(
                        candidate.rule, state.gadget, candidate.gadget,
                    )
                ]
                nv(candidate.gadget.replacement_graph) == min_vertices &&
                    return UnweightedOptimizationResult(
                        candidate.gadget, candidate_steps, sat_evaluations,
                        :minimum_vertices,
                    )
                candidate_key = _rewrite_state_key(candidate.gadget)
                candidate_key in seen_states && continue
                push!(seen_states, candidate_key)
                push!(direct_states, (;
                    gadget=candidate.gadget, steps=candidate_steps,
                ))
            end
        end

        best = argmin(
            state -> nv(state.gadget.replacement_graph), direct_states,
        )
        sources = [
            sort(direct_states[2:end]; by=state -> -nv(state.gadget.replacement_graph));
            direct_states[1]
        ]
        resumed = false
        for source in sources
            candidate, used, status = _resynthesized_unweighted_rewrite(
                source.gadget, target_boundary, target_reduced, lattice, min_vertices,
                max_sat_evaluations - sat_evaluations, host_radius,
            )
            sat_evaluations += used
            if candidate !== nothing
                current = candidate.gadget
                steps = [
                    source.steps;
                    UnweightedRewriteStep(
                        candidate.rule, source.gadget, candidate.gadget,
                    )
                ]
                resumed = true
                break
            end
            status == :budget && return UnweightedOptimizationResult(
                best.gadget, best.steps, sat_evaluations, :sat_budget,
            )
        end
        resumed || return UnweightedOptimizationResult(
            best.gadget, best.steps, sat_evaluations, :rewrite_fixed_point,
        )
    end
end

_rewrite_state_key(gadget) = (
    Tuple(gadget.lattice_coordinates),
    Tuple(gadget.lattice_coordinates[gadget.boundary_vertices]),
    Tuple(gadget.pin_rays),
)

_gadget_lattice(gadget::UnweightedGadget) =
    gadget.lattice == :KSG ? Square() :
    gadget.lattice == :triangular ? Triangular() :
    error("unknown gadget lattice $(gadget.lattice)")

function _gadget_patch(gadget, lattice)
    pins = gadget.lattice_coordinates[gadget.boundary_vertices]
    directions = _lattice_directions(lattice)
    rays = [_lattice_direction_index(directions, ray) for ray in gadget.pin_rays]
    return _LatticePatch(copy(gadget.lattice_coordinates), pins, rays)
end

function _rewrite_gadget(
    gadget, target_boundary, target_reduced, lattice, sites, pins, rays,
)
    analysis = _analyze_crossing_candidate(
        target_reduced, lattice, collect(sites), collect(pins), collect(rays),
    )
    analysis.solved || return nothing
    return _certified_rewrite_gadget(gadget, target_boundary, lattice, analysis)
end

function _certified_rewrite_gadget(gadget, target_boundary, lattice, analysis)
    valid, offset = is_gadget_replacement(
        gadget.pattern_graph, analysis.graph, target_boundary, analysis.boundary,
    )
    valid || return nothing
    offset == analysis.offset || error("rewrite analysis and verifier disagree on the offset")
    return _unweighted_gadget(gadget.pattern_graph, lattice, analysis)
end

function _direct_unweighted_rewrites(
    gadget, target_boundary, target_reduced, lattice, min_vertices,
)
    patch = _gadget_patch(gadget, lattice)
    graph = gadget.replacement_graph
    boundary = gadget.boundary_vertices
    rewrites = NamedTuple[]

    nv(graph) - 2 >= min_vertices && for label in eachindex(boundary)
        pin = boundary[label]
        degree(graph, pin) == 1 || continue
        middle = only(neighbors(graph, pin))
        for endpoint in neighbors(graph, middle)
            endpoint == pin && continue
            endpoint in boundary && continue
            sites = setdiff(patch.coordinates, patch.coordinates[[pin, middle]])
            pins = copy(patch.pins)
            pins[label] = patch.coordinates[endpoint]
            rays = copy(patch.rays)
            rays[label] = _direction_index_between(
                lattice, pins[label], patch.coordinates[middle],
            )
            rewritten = _rewrite_gadget(
                gadget, target_boundary, target_reduced, lattice, sites, pins, rays,
            )
            rewritten === nothing || push!(rewrites, (;
                rule=:even_boundary_tail_contraction, gadget=rewritten,
            ))
        end
    end

    nv(graph) - 2 >= min_vertices && for (first_label, second_label) in ((1, 3), (2, 4))
        first_pin = boundary[first_label]
        second_pin = boundary[second_label]
        degree(graph, first_pin) == 1 || continue
        degree(graph, second_pin) == 1 || continue
        first_neighbor = only(neighbors(graph, first_pin))
        second_neighbor = only(neighbors(graph, second_pin))
        first_neighbor == second_neighbor && continue
        first_neighbor in boundary && continue
        second_neighbor in boundary && continue
        sites = setdiff(
            patch.coordinates, patch.coordinates[[first_pin, second_pin]],
        )
        pins = copy(patch.pins)
        pins[first_label] = patch.coordinates[first_neighbor]
        pins[second_label] = patch.coordinates[second_neighbor]
        rays = copy(patch.rays)
        rays[first_label] = _direction_index_between(
            lattice, pins[first_label], patch.coordinates[first_pin],
        )
        rays[second_label] = _direction_index_between(
            lattice, pins[second_label], patch.coordinates[second_pin],
        )
        rewritten = _rewrite_gadget(
            gadget, target_boundary, target_reduced, lattice, sites, pins, rays,
        )
        rewritten === nothing || push!(rewrites, (;
            rule=:opposite_leaf_pin_contraction, gadget=rewritten,
        ))
    end
    return rewrites
end

function _direction_index_between(lattice, source, destination)
    source_canonical = _canonical_coordinate(lattice, source)
    destination_canonical = _canonical_coordinate(lattice, destination)
    direction = (
        destination_canonical[1] - source_canonical[1],
        destination_canonical[2] - source_canonical[2],
    )
    return _lattice_direction_index(_lattice_directions(lattice), direction)
end

function _expanded_lattice_host(lattice, coordinates, radius)
    host = Set(coordinates)
    frontier = Set(coordinates)
    directions = _lattice_directions(lattice)
    for _ in 1:radius
        next_frontier = Set{_LatticeCoordinate}()
        for coordinate in frontier, direction in directions
            neighbor = _lattice_step(lattice, coordinate, direction, 1)
            neighbor in host || push!(next_frontier, neighbor)
        end
        union!(host, next_frontier)
        frontier = next_frontier
    end
    return sort!(collect(host))
end

function _one_step_frame_rewrites(lattice, patch)
    directions = _lattice_directions(lattice)
    frames = NamedTuple[]
    for label in eachindex(patch.pins), step in directions
        pins = copy(patch.pins)
        pins[label] = _lattice_step(lattice, pins[label], step, 1)
        length(unique(pins)) == 4 || continue
        for ray in eachindex(directions)
            rays = copy(patch.rays)
            rays[label] = ray
            _pin_rays_are_compatible(lattice, pins, rays) || continue
            _labels_alternate(lattice, pins, rays) || continue
            push!(frames, (; label, pins, rays))
        end
    end
    return frames
end

function _resynthesized_unweighted_rewrite(
    gadget, target_boundary, target_reduced, lattice, min_vertices,
    sat_budget, host_radius,
)
    iszero(sat_budget) && return nothing, 0, :budget
    patch = _gadget_patch(gadget, lattice)
    window = _expanded_lattice_host(lattice, patch.coordinates, host_radius)
    rewritten_frames = _one_step_frame_rewrites(lattice, patch)
    sort!(rewritten_frames; by=frame -> (
        degree(gadget.replacement_graph, gadget.boundary_vertices[frame.label]) == 1,
        frame.rays[frame.label] != patch.rays[frame.label],
        frame.label,
    ))
    frame_hosts = [begin
        allowed = _allowed_crossing_sites(
            lattice, window, frame.pins, frame.rays, min_vertices,
        )
        allowed === nothing ? nothing : (; frame, allowed)
    end for frame in rewritten_frames]
    filter!(!isnothing, frame_hosts)
    contexts = Dict{Int, _SatFrameContext}()
    completion = _target_completion(target_reduced)
    evaluated = 0
    for atom_count in nv(gadget.replacement_graph)-1:-1:min_vertices
        expected_offset = Int(gadget.constant_offset) -
            (nv(gadget.replacement_graph) - atom_count)
        offsets = _centered_offsets(
            -Int(minimum(completion)):
            atom_count-Int(maximum(completion)), expected_offset,
        )
        for offset in offsets
            for (frame_index, frame_host) in enumerate(frame_hosts)
                evaluated == sat_budget && return nothing, evaluated, :budget
                length(frame_host.allowed) >= atom_count || continue
                frame = (;
                    pins=frame_host.frame.pins,
                    rays=frame_host.frame.rays,
                    allowed=frame_host.allowed,
                )
                context = get!(contexts, frame_index) do
                    _prepare_sat_frame(lattice, frame)
                end
                evaluated += 1
                analysis = _solve_fixed_crossing_sat(
                    target_reduced, lattice, context, atom_count, offset,
                )
                analysis === nothing && continue
                rewritten = _certified_rewrite_gadget(
                    gadget, target_boundary, lattice, analysis,
                )
                rewritten === nothing &&
                    error("fixed-frame SAT result failed the gadget verifier")
                return (;
                    rule=:frame_rewrite_resynthesis,
                    gadget=rewritten,
                ), evaluated, :found
            end
        end
    end
    return nothing, evaluated, :exhausted
end

function _analyze_crossing_candidate(target_reduced, lattice, sites, pins, rays)
    patch = _LatticePatch(sort(copy(sites)), copy(pins), copy(rays))
    graph, boundary, positions = _materialize_lattice_patch(lattice, patch)
    reduced = vec(calculate_reduced_alpha_tensor(graph, boundary))
    tensor_valid, offset = is_diff_by_constant(reduced, target_reduced)
    geometry_valid = all(_check_crossing_frame(lattice, patch))
    solved = tensor_valid && geometry_valid && is_connected(graph)
    return (;
        graph, boundary, positions, patch, solved, offset=Float64(offset),
    )
end

function _unweighted_gadget(target_graph, lattice, analysis)
    return UnweightedGadget(
        target_graph, analysis.graph, analysis.boundary, analysis.offset,
        _lattice_symbol(lattice), copy(analysis.patch.coordinates),
        analysis.positions, _patch_ray_directions(lattice, analysis.patch),
    )
end

_rotate_lattice_coordinate(::Triangular, point) = (-point[2], point[1] + point[2])
_reflect_lattice_coordinate(::Triangular, point) = (point[2], point[1])
_rotate_lattice_coordinate(::Square, point) = (-point[2], point[1])
_reflect_lattice_coordinate(::Square, point) = (point[1], -point[2])

function _transform_lattice_coordinate(lattice, point, rotations, reflected)
    transformed = reflected ? _reflect_lattice_coordinate(lattice, point) : point
    for _ in 1:rotations
        transformed = _rotate_lattice_coordinate(lattice, transformed)
    end
    return transformed
end

function _canonical_crossing_frame_key(lattice, frame)
    allowed = _canonical_coordinate.(Ref(lattice), frame.allowed)
    pins = _canonical_coordinate.(Ref(lattice), frame.pins)
    rays = _lattice_directions(lattice)[frame.rays]
    return minimum((begin
        transformed_allowed = _transform_lattice_coordinate.(
            Ref(lattice), allowed, rotations, reflected,
        )
        transformed_pins = _transform_lattice_coordinate.(
            Ref(lattice), pins, rotations, reflected,
        )
        transformed_rays = _transform_lattice_coordinate.(
            Ref(lattice), rays, rotations, reflected,
        )
        minimum_first = minimum(first, transformed_allowed)
        minimum_last = minimum(last, transformed_allowed)
        normalize(point) = (point[1] - minimum_first, point[2] - minimum_last)
        (
            Tuple(sort(normalize.(transformed_allowed))),
            Tuple(normalize.(transformed_pins)),
            Tuple(transformed_rays),
        )
    end for reflected in (false, true) for rotations in
        0:(lattice isa Triangular ? 5 : 3)))
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

function _lattice_symbol(lattice::LatticeType)
    return lattice isa Square ? :KSG : :triangular
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
    length(unique(coordinates)) == length(coordinates) ||
        throw(ArgumentError("crossing-frame coordinates must be unique"))
    length(unique(pins)) == 4 ||
        throw(ArgumentError("crossing-frame pins must be distinct"))
    all(in(Set(coordinates)), pins) ||
        throw(ArgumentError("every crossing-frame pin must be present in coordinates"))
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
