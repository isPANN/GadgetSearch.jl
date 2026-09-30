mutable struct _SatCnf
    variables::Int
    clauses::Vector{Vector{Int}}
end

struct _SatFrameContext
    frame
    coordinates::Vector{_LatticeCoordinate}
    boundary::Vector{Int}
    edge_list::Vector{Tuple{Int, Int}}
    order::Vector{Int}
    layers::Vector{Vector{Tuple{Int, Int, Bool}}}
    adjacent::Vector{Vector{Int}}
end

_SatCnf() = _SatCnf(0, Vector{Int}[])

function _sat_variable!(cnf::_SatCnf)
    cnf.variables += 1
    return cnf.variables
end

_sat_clause!(cnf::_SatCnf, literals::Int...) = push!(cnf.clauses, collect(literals))

function _sat_at_most!(cnf::_SatCnf, literals::Vector{Int}, bound::Int)
    bound >= length(literals) && return
    bound >= 0 || return _sat_clause!(cnf)
    if bound == 0
        foreach(literal -> _sat_clause!(cnf, -literal), literals)
        return
    end
    counters = [
        [_sat_variable!(cnf) for _ in 1:bound]
        for _ in 1:length(literals)-1
    ]
    for row in 1:length(literals)-1
        _sat_clause!(cnf, -literals[row], counters[row][1])
    end
    for row in 2:length(literals)-1
        _sat_clause!(cnf, -counters[row-1][1], counters[row][1])
    end
    for column in 2:bound
        for row in column:length(literals)-1
            _sat_clause!(
                cnf, -literals[row], -counters[row-1][column-1],
                counters[row][column],
            )
        end
        for row in column+1:length(literals)-1
            _sat_clause!(cnf, -counters[row-1][column], counters[row][column])
        end
    end
    for row in bound+1:length(literals)
        _sat_clause!(cnf, -literals[row], -counters[row-1][bound])
    end
end

function _sat_exactly!(cnf::_SatCnf, literals::Vector{Int}, count::Int)
    _sat_at_most!(cnf, literals, count)
    _sat_at_most!(cnf, -literals, length(literals) - count)
end

function _sat_exactly_one!(cnf::_SatCnf, literals::Vector{Int})
    _sat_clause!(cnf, literals...)
    _sat_at_most!(cnf, literals, 1)
end

function _independent_set_bdd(coordinates, edge_list)
    adjacent = [Set{Int}() for _ in coordinates]
    for (first, second) in edge_list
        push!(adjacent[first], second)
        push!(adjacent[second], first)
    end
    canonical = _offset_to_axial.(coordinates)
    orderings = Vector{Vector{Int}}()
    for projection in (
        point -> (point[1], point[2]),
        point -> (point[2], point[1]),
        point -> (point[1] + point[2], point[1]),
    )
        order = sortperm(canonical; by=projection)
        push!(orderings, order, reverse(order))
    end
    function build_layers(order)
        position = zeros(Int, length(coordinates))
        for (step, vertex) in enumerate(order)
            position[vertex] = step
        end
        last_neighbor = [maximum([position[vertex]; position[collect(adjacent[vertex])]])
            for vertex in eachindex(coordinates)]
        layer = [Int[]]
        layers = Vector{Vector{Tuple{Int, Int, Bool}}}()
        layer_sizes = [1]
        for (step, vertex) in enumerate(order)
            next_layer = Vector{Vector{Int}}()
            next_index = Dict{Tuple{Vararg{Int}}, Int}()
            arcs = Tuple{Int, Int, Bool}[]
            for (source, occupied) in enumerate(layer), take in (false, true)
                take && any(neighbor -> neighbor in occupied, adjacent[vertex]) && continue
                next_occupied = take ? [occupied; vertex] : occupied
                frontier = Tuple(sort!(filter(v -> last_neighbor[v] > step, next_occupied)))
                destination = get!(next_index, frontier) do
                    push!(next_layer, collect(frontier))
                    length(next_layer)
                end
                push!(arcs, (source, destination, take))
            end
            push!(layers, arcs)
            push!(layer_sizes, length(next_layer))
            layer = next_layer
        end
        return sum(layer_sizes), order, layers
    end
    best = build_layers(first(orderings))
    for order in Iterators.drop(orderings, 1)
        candidate = build_layers(order)
        candidate[1] < best[1] && (best = candidate)
    end
    _, order, layers = best
    return order, layers
end

function _target_completion(target)
    return [maximum(target[subset+1] for subset in 0:state
        if subset & ~state == 0 && isfinite(target[subset+1]))
        for state in 0:length(target)-1]
end

function _essential_lower_states(target)
    states = Int[]
    for state in 0:length(target)-1
        value = target[state+1]
        isfinite(value) || continue
        implied = [target[superset+1] - count_ones(superset ⊻ state)
            for superset in state+1:length(target)-1
            if state & ~superset == 0 && isfinite(target[superset+1])]
        (isempty(implied) || value > maximum(implied)) && push!(states, state)
    end
    return states
end

function _essential_upper_states(completion)
    states = Int[]
    for state in 0:length(completion)-1
        implied = [completion[subset+1] + count_ones(state ⊻ subset)
            for subset in 0:state-1 if subset & ~state == 0]
        (isempty(implied) || completion[state+1] < minimum(implied)) &&
            push!(states, state)
    end
    return states
end

function _centered_offsets(offsets, center=(first(offsets) + last(offsets)) / 2)
    return sort!(collect(offsets); by=offset -> (abs(offset - center), offset))
end

function _add_independent_set_upper_bound!(
    cnf, selected, pins, order, layers, boundary_state, bound,
)
    root = _sat_variable!(cnf)
    _sat_clause!(cnf, root)
    reachable = Dict((1, 0) => root)
    pin_slots = Dict(pin => slot - 1 for (slot, pin) in enumerate(pins))
    for (step, (vertex, arcs)) in enumerate(zip(order, layers))
        next_reachable = Dict{Tuple{Int, Int}, Int}()
        required = haskey(pin_slots, vertex) ?
            !iszero(boundary_state & (1 << pin_slots[vertex])) : nothing
        for (source, destination, take) in arcs
            required !== nothing && take != required && continue
            for ((source_state, count), source_variable) in reachable
                source_state == source || continue
                next_count = min(bound + 1, count + Int(take))
                destination_variable = get!(next_reachable, (destination, next_count)) do
                    _sat_variable!(cnf)
                end
                take ? _sat_clause!(
                    cnf, -source_variable, -selected[vertex], destination_variable,
                ) : _sat_clause!(cnf, -source_variable, destination_variable)
            end
        end
        reachable = next_reachable
    end
    for ((_, count), variable) in reachable
        count == bound + 1 && _sat_clause!(cnf, -variable)
    end
end

function _add_selected_connectivity!(
    cnf, selected, adjacent, root, selected_count,
)
    reachable = [_sat_variable!(cnf) for _ in selected]
    for vertex in eachindex(reachable)
        _sat_clause!(cnf, vertex == root ? reachable[vertex] : -reachable[vertex])
    end
    for _ in 1:selected_count-1
        next_reachable = [_sat_variable!(cnf) for _ in selected]
        for vertex in eachindex(selected)
            _sat_clause!(cnf, -reachable[vertex], next_reachable[vertex])
            _sat_clause!(cnf, -next_reachable[vertex], selected[vertex])
            _sat_clause!(
                cnf, -next_reachable[vertex], reachable[vertex],
                (reachable[neighbor] for neighbor in adjacent[vertex])...,
            )
            for neighbor in adjacent[vertex]
                _sat_clause!(
                    cnf, -selected[vertex], -reachable[neighbor],
                    next_reachable[vertex],
                )
            end
        end
        reachable = next_reachable
    end
    for vertex in eachindex(selected)
        _sat_clause!(cnf, -selected[vertex], reachable[vertex])
    end
end

function _prepare_sat_frame(lattice, frame)
    coordinates = sort(copy(frame.allowed))
    patch = _LatticePatch(coordinates, copy(frame.pins), copy(frame.rays))
    host, boundary, _ = _materialize_lattice_patch(lattice, patch)
    edge_list = [(src(edge), dst(edge)) for edge in edges(host)]
    order, layers = _independent_set_bdd(coordinates, edge_list)
    adjacent = [collect(neighbors(host, vertex)) for vertex in vertices(host)]
    return _SatFrameContext(
        frame, coordinates, boundary, edge_list, order, layers, adjacent,
    )
end

function _solve_fixed_crossing_sat(target, lattice, frame, atom_count, offset)
    return _solve_fixed_crossing_sat(
        target, lattice, _prepare_sat_frame(lattice, frame), atom_count, offset,
    )
end

function _solve_fixed_crossing_sat(
    target, lattice, context::_SatFrameContext, atom_count, offset,
)
    solver, selected = _fixed_crossing_sat_problem(
        target, context, atom_count, offset,
    )
    return _solve_next_fixed_crossing_sat!(
        solver, selected, target, lattice, context,
    )
end

function _fixed_crossing_sat_problem(target, context, atom_count, offset)
    cnf, selected = _fixed_crossing_sat_cnf(
        target, context, atom_count, offset,
    )
    return _new_sat_solver(cnf), selected
end

function _fixed_crossing_sat_cnf(target, context, atom_count, offset)
    coordinates = context.coordinates
    pins = context.boundary
    cnf = _SatCnf()
    selected = [_sat_variable!(cnf) for _ in coordinates]
    foreach(pin -> _sat_clause!(cnf, selected[pin]), pins)
    _sat_exactly!(cnf, selected, atom_count)

    for state in _essential_lower_states(target)
        _add_lower_state_constraint!(
            cnf, selected, target, context, state, offset,
        )
    end

    completion = _target_completion(target)
    for state in _essential_upper_states(completion)
        _add_independent_set_upper_bound!(
            cnf, selected, pins, context.order, context.layers, state,
            Int(completion[state+1] + offset),
        )
    end
    _add_selected_connectivity!(
        cnf, selected, context.adjacent, pins[1], atom_count,
    )
    return cnf, selected
end

function _add_lower_state_constraint!(
    cnf, selected, target, context, state, offset,
)
    witness = [_sat_variable!(cnf) for _ in context.coordinates]
    for vertex in eachindex(context.coordinates)
        _sat_clause!(cnf, -witness[vertex], selected[vertex])
    end
    for (slot, pin) in enumerate(context.boundary)
        _sat_clause!(
            cnf, iszero(state & (1 << (slot - 1))) ?
            -witness[pin] : witness[pin],
        )
    end
    for (first, second) in context.edge_list
        _sat_clause!(cnf, -witness[first], -witness[second])
    end
    _sat_exactly!(cnf, witness, Int(target[state+1] + offset))
end

function _add_joint_lower_state_constraint!(
    cnf, selected, target, edge_list, frame_choices, state, offset,
)
    witness = [_sat_variable!(cnf) for _ in selected]
    for vertex in eachindex(selected)
        _sat_clause!(cnf, -witness[vertex], selected[vertex])
    end
    for label in eachindex(frame_choices),
        (vertex, _, port) in frame_choices[label]
            occupied = !iszero(state & (1 << (label - 1)))
            _sat_clause!(cnf, -port,
                occupied ? witness[vertex] : -witness[vertex])
    end
    for (first, second) in edge_list
        _sat_clause!(cnf, -witness[first], -witness[second])
    end
    _sat_exactly!(cnf, witness, Int(target[state+1] + offset))
end

function _add_joint_independent_set_upper_bound!(
    cnf, selected, choices_by_slot_vertex, order, layers, boundary_state, bound,
)
    maximum_remaining = [Dict{Int, Int}() for _ in 1:length(order)+1]
    maximum_remaining[end][1] = 0
    for step in length(order):-1:1
        for (source, destination, take) in layers[step]
            haskey(maximum_remaining[step+1], destination) || continue
            value = Int(take) + maximum_remaining[step+1][destination]
            maximum_remaining[step][source] = max(
                get(maximum_remaining[step], source, -1), value,
            )
        end
    end
    root = _sat_variable!(cnf)
    _sat_clause!(cnf, root)
    reachable = [(1, 0, root)]
    for (step, (vertex, arcs)) in enumerate(zip(order, layers))
        required_taken = Int[]
        required_skipped = Int[]
        for label in eachindex(choices_by_slot_vertex)
            destination = iszero(boundary_state & (1 << (label - 1))) ?
                required_skipped : required_taken
            append!(destination, choices_by_slot_vertex[label][vertex])
        end
        next_reachable = Tuple{Int, Int, Int}[]
        next_index = Dict{Tuple{Int, Int}, Int}()
        for (source, destination, take) in arcs
            for (source_state, count, source_variable) in reachable
                source_state == source || continue
                next_count = min(bound + 1, count + Int(take))
                next_count + get(maximum_remaining[step+1], destination, -1) <
                    bound + 1 && continue
                key = (destination, next_count)
                destination_variable = get(next_index, key, 0)
                if iszero(destination_variable)
                    destination_variable = _sat_variable!(cnf)
                    next_index[key] = destination_variable
                    push!(next_reachable, (destination, next_count, destination_variable))
                end
                if take
                    _sat_clause!(
                        cnf, -source_variable, -selected[vertex],
                        required_skipped..., destination_variable,
                    )
                else
                    _sat_clause!(
                        cnf, -source_variable, required_taken...,
                        destination_variable,
                    )
                end
            end
        end
        reachable = next_reachable
    end
    for (_, count, variable) in reachable
        count == bound + 1 && _sat_clause!(cnf, -variable)
    end
end

function _joint_frame_choice_is_outward(lattice, interface, direction, interfaces)
    pin = (interface[1] - direction[1], interface[2] - direction[2])
    geometry = _geometry_coordinate.(Ref(lattice), _from_canonical.(Ref(lattice), interfaces))
    pin_geometry = _geometry_coordinate(lattice, _from_canonical(lattice, pin))
    hull = _strict_convex_hull([geometry; pin_geometry])
    all(in(hull), geometry) || return false
    sum_q = sum(first, interfaces)
    sum_r = sum(last, interfaces)
    out_q = 4interface[1] - sum_q
    out_r = 4interface[2] - sum_r
    if lattice isa Square
        return out_q * direction[1] + out_r * direction[2] > 0
    end
    out_x = 2out_q + out_r
    direction_x = 2direction[1] + direction[2]
    return out_x * direction_x + 3out_r * direction[2] > 0
end

function _sat_group_indicator!(cnf, literals)
    indicator = _sat_variable!(cnf)
    for literal in literals
        _sat_clause!(cnf, -literal, indicator)
    end
    _sat_clause!(cnf, -indicator, literals...)
    return indicator
end

function _interfaces_form_alternating_quadrilateral(interfaces)
    length(unique(interfaces)) == 4 || return false
    return (
        _orientation(interfaces[1], interfaces[3], interfaces[2]) *
        _orientation(interfaces[1], interfaces[3], interfaces[4]) < 0 &&
        _orientation(interfaces[2], interfaces[4], interfaces[1]) *
        _orientation(interfaces[2], interfaces[4], interfaces[3]) < 0
    )
end

function _joint_rays_touch(start1, direction1, start2, direction2, directions)
    return any([_LatticeCoordinate[(0, 0)]; directions]) do offset
        _rays_touch(start1, direction1, start2, direction2, offset)
    end
end

function _joint_site_blocks_ray(site, pin, direction, directions)
    site == pin && return false
    start = (pin[1] + direction[1], pin[2] + direction[2])
    return any([_LatticeCoordinate[(0, 0)]; directions]) do offset
        _point_on_ray(
            (site[1] + offset[1], site[2] + offset[2]), start, direction,
        )
    end
end

function _add_joint_frame_choices!(
    cnf, selected, coordinates, lattice, first_direction_index,
)
    directions = _lattice_directions(lattice)
    canonical = _canonical_coordinate.(Ref(lattice), coordinates)
    origin = only(findall(==((0, 0)), canonical))
    frame_choices = [Tuple{Int, Int, Int}[] for _ in 1:4]
    push!(frame_choices[1], (
        origin, first_direction_index, _sat_variable!(cnf),
    ))
    fixed_start = directions[first_direction_index]
    for label in 2:4, vertex in eachindex(coordinates),
        direction_index in eachindex(directions)
        direction = directions[direction_index]
        pin = canonical[vertex]
        start = (pin[1] + direction[1], pin[2] + direction[2])
        (vertex == origin || _joint_rays_touch(
            fixed_start, directions[first_direction_index], start, direction,
            directions,
        )) && continue
        push!(frame_choices[label],
            (vertex, direction_index, _sat_variable!(cnf)))
    end

    choices_by_slot_vertex = [
        [Int[] for _ in coordinates] for _ in 1:4
    ]
    for label in eachindex(frame_choices)
        _sat_exactly_one!(cnf, last.(frame_choices[label]))
        for (vertex, _, variable) in frame_choices[label]
            push!(choices_by_slot_vertex[label][vertex], variable)
            _sat_clause!(cnf, -variable, selected[vertex])
        end
    end
    for first_label in 1:3, second_label in first_label+1:4,
        vertex in eachindex(coordinates),
        first in choices_by_slot_vertex[first_label][vertex],
        second in choices_by_slot_vertex[second_label][vertex]
        _sat_clause!(cnf, -first, -second)
    end

    interface_variables = [Dict{_LatticeCoordinate, Int}() for _ in 1:4]
    interface_order = [_LatticeCoordinate[] for _ in 1:4]
    records = [
        Dict{_LatticeCoordinate, Vector{Tuple{Int, Int}}}() for _ in 1:4
    ]
    for label in eachindex(frame_choices)
        for (vertex, direction_index, choice) in frame_choices[label]
            pin = canonical[vertex]
            direction = directions[direction_index]
            interface = (pin[1] + direction[1], pin[2] + direction[2])
            if !haskey(records[label], interface)
                records[label][interface] = Tuple{Int, Int}[]
                push!(interface_order[label], interface)
            end
            push!(records[label][interface], (direction_index, choice))
        end
        for interface in interface_order[label]
            choices = last.(records[label][interface])
            variable = _sat_group_indicator!(cnf, choices)
            interface_variables[label][interface] = variable
        end
    end

    first_interface = only(interface_order[1])
    adjacent_pairs = Tuple{_LatticeCoordinate, _LatticeCoordinate, Int}[]
    for adjacent_first in interface_order[2], adjacent_second in interface_order[4]
        adjacent_first < adjacent_second || continue
        pair = _sat_variable!(cnf)
        first_variable = interface_variables[2][adjacent_first]
        second_variable = interface_variables[4][adjacent_second]
        _sat_clause!(cnf, -pair, first_variable)
        _sat_clause!(cnf, -pair, second_variable)
        _sat_clause!(cnf, -first_variable, -second_variable, pair)
        push!(adjacent_pairs, (adjacent_first, adjacent_second, pair))
    end
    for opposite in interface_order[3]
        opposite_variable = interface_variables[3][opposite]
        allowed_pairs = Int[]
        for (adjacent_first, adjacent_second, pair) in adjacent_pairs
            interfaces = [
                first_interface, adjacent_first, opposite, adjacent_second,
            ]
            _interfaces_form_alternating_quadrilateral(interfaces) || continue
            allowed_choices = [
                [choice for (direction_index, choice) in records[label][interface]
                    if _joint_frame_choice_is_outward(
                        lattice, interface, directions[direction_index], interfaces,
                    )]
                for (label, interface) in enumerate(interfaces)
            ]
            any(isempty, allowed_choices) && continue
            push!(allowed_pairs, pair)
            for choices in allowed_choices[2:4]
                _sat_clause!(cnf, -opposite_variable, -pair, choices...)
            end
        end
        _sat_clause!(cnf, -opposite_variable, allowed_pairs...)
    end

    for choices in frame_choices, (vertex, direction_index, variable) in choices
        pin = canonical[vertex]
        direction = directions[direction_index]
        for (site_vertex, site) in enumerate(canonical)
            _joint_site_blocks_ray(site, pin, direction, directions) &&
                _sat_clause!(cnf, -variable, -selected[site_vertex])
        end
    end
    for first_label in 1:3, second_label in first_label+1:4,
        (first_vertex, first_direction, first_variable) in frame_choices[first_label],
        (second_vertex, second_direction, second_variable) in frame_choices[second_label]
        first_pin = canonical[first_vertex]
        second_pin = canonical[second_vertex]
        first_ray = directions[first_direction]
        second_ray = directions[second_direction]
        first_start = (first_pin[1] + first_ray[1], first_pin[2] + first_ray[2])
        second_start = (second_pin[1] + second_ray[1], second_pin[2] + second_ray[2])
        _joint_rays_touch(
            first_start, first_ray, second_start, second_ray, directions,
        ) && _sat_clause!(cnf, -first_variable, -second_variable)
    end
    return frame_choices, choices_by_slot_vertex
end

function _joint_crossing_sat_cnf(
    target, lattice, shape, atom_count, offset;
    canonical_shift=(0, 0), first_direction_index=1,
)
    columns, rows = shape
    coordinates = lattice isa Triangular ? sort!(_from_canonical.(
        Ref(lattice), _LatticeCoordinate[
            (q + canonical_shift[1], r + canonical_shift[2])
            for q in -1:columns-2 for r in 1-rows:0
        ],
    )) : _LatticeCoordinate[
        (column + canonical_shift[1], row + canonical_shift[2])
        for column in 0:columns-1 for row in 0:rows-1
    ]
    canonical = _canonical_coordinate.(Ref(lattice), coordinates)
    coordinate_index = Dict(point => vertex for (vertex, point) in enumerate(canonical))
    directions = _lattice_directions(lattice)
    edge_list = Tuple{Int, Int}[]
    for (vertex, point) in enumerate(canonical), direction in directions
        other = get(coordinate_index,
            (point[1] + direction[1], point[2] + direction[2]), 0)
        vertex < other && push!(edge_list, (vertex, other))
    end
    adjacent = [Int[] for _ in coordinates]
    for (first, second) in edge_list
        push!(adjacent[first], second)
        push!(adjacent[second], first)
    end
    order, layers = _independent_set_bdd(coordinates, edge_list)

    cnf = _SatCnf()
    selected = [_sat_variable!(cnf) for _ in coordinates]
    _sat_exactly!(cnf, selected, atom_count)
    origin = only(findall(==((0, 0)), canonical))
    _add_selected_connectivity!(cnf, selected, adjacent, origin, atom_count)
    frame_choices, choices_by_slot_vertex = _add_joint_frame_choices!(
        cnf, selected, coordinates, lattice, first_direction_index,
    )
    for state in _essential_lower_states(target)
        _add_joint_lower_state_constraint!(
            cnf, selected, target, edge_list, frame_choices, state, offset,
        )
    end
    completion = _target_completion(target)
    for state in _essential_upper_states(completion)
        _add_joint_independent_set_upper_bound!(
            cnf, selected, choices_by_slot_vertex, order, layers, state,
            Int(completion[state+1] + offset),
        )
    end
    return cnf, selected, frame_choices, coordinates
end

function _solve_joint_crossing_sat(
    target, lattice, shape, atom_count, offset;
    verbose=false, seconds=600, kissat_executable=nothing, initial_seed=1,
    canonical_shift=(0, 0), first_direction_index=1,
)
    cnf, selected, frame_choices, coordinates = _joint_crossing_sat_cnf(
        target, lattice, shape, atom_count, offset;
        canonical_shift, first_direction_index,
    )
    solver = _new_sat_solver(cnf)
    choice_variables = [choice[3] for choices in frame_choices for choice in choices]
    projected = [selected; choice_variables]
    model_index = 0
    deadline = time() + seconds
    while time() < deadline
        model_index += 1
        status, assignment = _next_joint_assignment!(
            solver, projected; seed=initial_seed + model_index - 1,
            seconds=max(1, floor(Int, deadline - time())),
            kissat_executable,
        )
        status == :unknown && return nothing
        status == :unsat && return nothing
        selected_assignment = assignment[1:length(selected)]
        choice_assignment = assignment[length(selected)+1:end]
        pins = _LatticeCoordinate[]
        rays = Int[]
        chosen_choices = Tuple{Int, Int, Int}[]
        cursor = 0
        for choices in frame_choices
            chosen = only(choice for choice in choices if choice_assignment[cursor += 1])
            push!(chosen_choices, chosen)
            push!(pins, coordinates[chosen[1]])
            push!(rays, chosen[2])
        end
        selected_indices = findall(identity, selected_assignment)
        sites = coordinates[selected_indices]
        checks = _check_gadget_geometry(lattice, _LatticePatch(sites, pins, rays))
        if all(checks)
            analysis = _analyze_crossing_candidate(
                target, lattice, sites, pins, rays,
            )
            analysis.solved && return analysis
            error("joint SAT candidate failed a non-geometric encoded constraint")
        end
        chosen_ports = last.(chosen_choices)
        intrinsic = _check_gadget_geometry(
            lattice, _LatticePatch(pins, pins, rays),
        )
        clause = all(intrinsic) ?
            [-selected[selected_indices]; -chosen_ports] :
            -chosen_ports
        _add_solver_clause!(solver, clause)
        verbose && println((; model_index, frame_checks=checks,
            intrinsic_frame_checks=intrinsic, blocker_length=length(clause)))
        verbose && flush(stdout)
    end
    return nothing
end

function _solve_next_fixed_crossing_sat!(
    solver, selected, target, lattice, context,
)
    assignment = _next_selected_assignment!(solver, selected)
    assignment === nothing && return nothing
    chosen = findall(identity, assignment)
    sites = context.coordinates[chosen]
    frame = context.frame
    analysis = _analyze_crossing_candidate(
        target, lattice, sites, frame.pins, frame.rays,
    )
    analysis.solved || error("SAT candidate failed its encoded constraints")
    return analysis
end

mutable struct _KissatEnumerator
    variables::Int
    clauses::Vector{Vector{Int}}
end

_new_sat_solver(cnf::_SatCnf) =
    _KissatEnumerator(cnf.variables, copy(cnf.clauses))
_add_solver_clause!(solver::_KissatEnumerator, clause) =
    push!(solver.clauses, collect(clause))
_kissat_init() = ccall((:kissat_init, Kissat_jll.libkissat), Ptr{Cvoid}, ())
_kissat_add(solver, literal) = ccall((:kissat_add, Kissat_jll.libkissat), Cvoid, (Ptr{Cvoid}, Cint), solver, literal)
_kissat_solve(solver) = ccall((:kissat_solve, Kissat_jll.libkissat), Cint, (Ptr{Cvoid},), solver)
_kissat_value(solver, variable) = ccall((:kissat_value, Kissat_jll.libkissat), Cint, (Ptr{Cvoid}, Cint), solver, variable)
_kissat_release(solver) = ccall((:kissat_release, Kissat_jll.libkissat), Cvoid, (Ptr{Cvoid},), solver)
_kissat_quiet(solver) = ccall((:kissat_set_option, Kissat_jll.libkissat), Cint, (Ptr{Cvoid}, Cstring, Cint), solver, "quiet", 1)
_kissat_seed(solver, seed) = ccall(
    (:kissat_set_option, Kissat_jll.libkissat), Cint,
    (Ptr{Cvoid}, Cstring, Cint), solver, "seed", seed,
)
_kissat_walkinitially(solver) = ccall(
    (:kissat_set_option, Kissat_jll.libkissat), Cint,
    (Ptr{Cvoid}, Cstring, Cint), solver, "walkinitially", 1,
)
_kissat_set_conflict_limit(solver, limit) = ccall(
    (:kissat_set_conflict_limit, Kissat_jll.libkissat), Cint,
    (Ptr{Cvoid}, Cuint), solver, limit,
)

function _next_selected_assignment!(enumerator::_KissatEnumerator, selected)
    solver = _kissat_init()
    _kissat_quiet(solver)
    for clause in enumerator.clauses
        foreach(literal -> _kissat_add(solver, literal), clause)
        _kissat_add(solver, 0)
    end
    status = _kissat_solve(solver)
    assignment = status == 10 ?
        Bool[_kissat_value(solver, variable) > 0 for variable in selected] : nothing
    _kissat_release(solver)
    status == 20 && return nothing
    status == 10 || error("SAT solver returned an undefined result")
    push!(enumerator.clauses, [
        assignment[index] ? -variable : variable
        for (index, variable) in enumerate(selected)
    ])
    return assignment
end

function _next_joint_assignment!(
    enumerator::_KissatEnumerator, variables;
    seed, seconds, kissat_executable,
)
    path, stream = mktemp()
    try
        println(stream, "p cnf $(enumerator.variables) $(length(enumerator.clauses))")
        for clause in enumerator.clauses
            println(stream, join(clause, ' '), " 0")
        end
        close(stream)
        output = IOBuffer()
        executable = isnothing(kissat_executable) ?
            Kissat_jll.kissat() : kissat_executable
        command = `$executable --sat --walkinitially --seed=$seed -q --time=$seconds $path`
        process = run(pipeline(ignorestatus(command), stdout=output, stderr=stderr))
        status = process.exitcode
        status == 20 && return :unsat, nothing
        status == 0 && return :unknown, nothing
        status == 10 || error("Kissat exited with status $status")
        positive = Set{Int}()
        for line in eachline(seekstart(output))
            startswith(line, "v ") || continue
            for literal in split(line)[2:end]
                value = parse(Int, literal)
                value > 0 && push!(positive, value)
            end
        end
        return :sat, Bool[variable in positive for variable in variables]
    finally
        isopen(stream) && close(stream)
        rm(path)
    end
end


function _next_selected_assignment_limited!(
    enumerator::_KissatEnumerator, selected, conflict_limit; seed=1,
)
    solver = _kissat_init()
    _kissat_quiet(solver)
    _kissat_seed(solver, seed)
    for clause in enumerator.clauses
        foreach(literal -> _kissat_add(solver, literal), clause)
        _kissat_add(solver, 0)
    end
    _kissat_set_conflict_limit(solver, conflict_limit)
    status = _kissat_solve(solver)
    assignment = status == 10 ?
        Bool[_kissat_value(solver, variable) > 0 for variable in selected] : nothing
    _kissat_release(solver)
    status == 0 && return :unknown, nothing
    status == 20 && return :unsat, nothing
    status == 10 || error("SAT solver returned an undefined result")
    push!(enumerator.clauses, [
        assignment[index] ? -variable : variable
        for (index, variable) in enumerate(selected)
    ])
    return :sat, assignment
end

function _pin_rays_are_compatible(lattice, pins, rays)
    directions = _lattice_directions(lattice)[collect(rays)]
    interfaces = [
        _lattice_step(lattice, pin, direction, 1)
        for (pin, direction) in zip(pins, directions)
    ]
    starts = _canonical_coordinate.(Ref(lattice), interfaces)
    adjacency_offsets = [_LatticeCoordinate[(0, 0)]; _lattice_directions(lattice)]
    for first in 1:3, second in first+1:4, offset in adjacency_offsets
        _rays_touch(
            starts[first], directions[first], starts[second], directions[second],
            offset,
        ) && return false
    end
    return true
end

function _labels_alternate(lattice, pins, rays)
    directions = _lattice_directions(lattice)[rays]
    interfaces = [
        _lattice_step(lattice, pin, direction, 1)
        for (pin, direction) in zip(pins, directions)
    ]
    interface_geometry = _geometry_coordinate.(Ref(lattice), interfaces)
    occupied_geometry = _geometry_coordinate.(Ref(lattice), pins)
    hull = _strict_convex_hull([occupied_geometry; interface_geometry])
    return _interfaces_alternate(hull, interface_geometry)
end

function _allowed_crossing_sites(lattice, window, pins, rays, min_allowed)
    directions = _lattice_directions(lattice)[rays]
    interfaces = [
        _lattice_step(lattice, pin, direction, 1)
        for (pin, direction) in zip(pins, directions)
    ]
    interface_geometry = _geometry_coordinate.(Ref(lattice), interfaces)
    pin_geometry = _geometry_coordinate.(Ref(lattice), pins)
    starts = _canonical_coordinate.(Ref(lattice), interfaces)
    adjacency_offsets = [_LatticeCoordinate[(0, 0)]; _lattice_directions(lattice)]
    allowed = copy(pins)
    candidates = [site for site in window if site ∉ pins]
    for (candidate_index, site) in enumerate(candidates)
        canonical_site = _canonical_coordinate(lattice, site)
        corridor_clear = all(eachindex(starts)) do ray_index
            all(adjacency_offsets) do offset
                !_point_on_ray(
                    (canonical_site[1] + offset[1], canonical_site[2] + offset[2]),
                    starts[ray_index], directions[ray_index],
                )
            end
        end
        if corridor_clear
            site_geometry = _geometry_coordinate(lattice, site)
            hull = _strict_convex_hull([
                pin_geometry; site_geometry; interface_geometry
            ])
            all(in(hull), interface_geometry) && push!(allowed, site)
        end
        remaining = length(candidates) - candidate_index
        length(allowed) + remaining >= min_allowed || return nothing
    end
    sort!(allowed)
    return allowed
end
