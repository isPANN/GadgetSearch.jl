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

function _independent_set_bdd(coordinates, edge_list)
    adjacent = [Set{Int}() for _ in coordinates]
    for (first, second) in edge_list
        push!(adjacent[first], second)
        push!(adjacent[second], first)
    end
    order = sortperm(coordinates)
    position = zeros(Int, length(coordinates))
    for (step, vertex) in enumerate(order)
        position[vertex] = step
    end
    last_neighbor = [maximum([position[vertex]; position[collect(adjacent[vertex])]])
        for vertex in eachindex(coordinates)]
    layer = [Int[]]
    layers = Vector{Vector{Tuple{Int, Int, Bool}}}()
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
        layer = next_layer
    end
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

function _add_independent_set_upper_bound!(
    cnf, selected, pins, order, layers, boundary_state, bound,
)
    root = _sat_variable!(cnf)
    _sat_clause!(cnf, root)
    reachable = Dict((1, 0) => root)
    pin_slots = Dict(pin => slot - 1 for (slot, pin) in enumerate(pins))
    for (vertex, arcs) in zip(order, layers)
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

function _add_selected_connectivity!(cnf, selected, adjacent, root, selected_count)
    reachable = [_sat_variable!(cnf) for _ in selected]
    for vertex in eachindex(reachable)
        _sat_clause!(cnf, vertex == root ? reachable[vertex] : -reachable[vertex])
    end
    for _ in 1:selected_count-1
        next_reachable = [_sat_variable!(cnf) for _ in selected]
        for vertex in eachindex(selected)
            arrivals = Int[]
            for neighbor in adjacent[vertex]
                arrival = _sat_variable!(cnf)
                push!(arrivals, arrival)
                _sat_clause!(cnf, -arrival, selected[vertex])
                _sat_clause!(cnf, -arrival, reachable[neighbor])
                _sat_clause!(cnf, -selected[vertex], -reachable[neighbor], arrival)
                _sat_clause!(cnf, -arrival, next_reachable[vertex])
            end
            _sat_clause!(cnf, -reachable[vertex], next_reachable[vertex])
            _sat_clause!(cnf, -next_reachable[vertex], reachable[vertex], arrivals...)
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
    coordinates = context.coordinates
    pins = context.boundary
    edge_list = context.edge_list
    order = context.order
    layers = context.layers
    cnf = _SatCnf()
    selected = [_sat_variable!(cnf) for _ in coordinates]
    foreach(pin -> _sat_clause!(cnf, selected[pin]), pins)
    _sat_exactly!(cnf, selected, atom_count)

    for state in _essential_lower_states(target)
        witness = [_sat_variable!(cnf) for _ in coordinates]
        for vertex in eachindex(coordinates)
            _sat_clause!(cnf, -witness[vertex], selected[vertex])
        end
        for (slot, pin) in enumerate(pins)
            _sat_clause!(
                cnf,
                iszero(state & (1 << (slot - 1))) ? -witness[pin] : witness[pin],
            )
        end
        for (first, second) in edge_list
            _sat_clause!(cnf, -witness[first], -witness[second])
        end
        _sat_exactly!(cnf, witness, Int(target[state+1] + offset))
    end

    completion = _target_completion(target)
    for state in _essential_upper_states(completion)
        _add_independent_set_upper_bound!(
            cnf, selected, pins, order, layers, state,
            Int(completion[state+1] + offset),
        )
    end
    _add_selected_connectivity!(
        cnf, selected, context.adjacent, pins[1], atom_count,
    )
    return _new_sat_solver(cnf), selected
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

@static if Sys.iswindows()
    function _new_sat_solver(cnf::_SatCnf)
        solver = CryptoMiniSat.CMS(cnf.variables; num_threads=1)
        foreach(clause -> CryptoMiniSat.add_clause(solver, clause), cnf.clauses)
        return solver
    end

    function _next_selected_assignment!(solver, selected)
        status = CryptoMiniSat.solve(solver)
        status === false && return nothing
        status === true || error("SAT solver returned an undefined result")
        model = CryptoMiniSat.get_model(solver)
        assignment = Bool[model[variable] for variable in selected]
        blocking_clause = [
            assignment[index] ? -variable : variable
            for (index, variable) in enumerate(selected)
        ]
        CryptoMiniSat.add_clause(solver, blocking_clause)
        return assignment
    end
else
    mutable struct _KissatEnumerator
        clauses::Vector{Vector{Int}}
    end

    _new_sat_solver(cnf::_SatCnf) = _KissatEnumerator(copy(cnf.clauses))
    _kissat_init() = ccall((:kissat_init, Kissat_jll.libkissat), Ptr{Cvoid}, ())
    _kissat_add(solver, literal) = ccall((:kissat_add, Kissat_jll.libkissat), Cvoid, (Ptr{Cvoid}, Cint), solver, literal)
    _kissat_solve(solver) = ccall((:kissat_solve, Kissat_jll.libkissat), Cint, (Ptr{Cvoid},), solver)
    _kissat_value(solver, variable) = ccall((:kissat_value, Kissat_jll.libkissat), Cint, (Ptr{Cvoid}, Cint), solver, variable)
    _kissat_release(solver) = ccall((:kissat_release, Kissat_jll.libkissat), Cvoid, (Ptr{Cvoid},), solver)
    _kissat_quiet(solver) = ccall((:kissat_set_option, Kissat_jll.libkissat), Cint, (Ptr{Cvoid}, Cstring, Cint), solver, "quiet", 1)

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
end

function _search_crossing_sat(
    target_graph, target_boundary, target_reduced, lattice;
    min_vertices, max_vertices, max_evaluations, max_frame_evaluations,
    max_results, window_side,
)
    all(value -> isinf(value) || isinteger(value), target_reduced) ||
        error("the reduced alpha tensor must contain integers or -Inf")
    completion = _target_completion(target_reduced)
    gadgets = UnweightedGadget[]
    seen_gadgets = Set{Tuple}()
    evaluated = 0
    frame_evaluated = Ref(0)
    for atom_count in min_vertices:max_vertices
        atom_count > window_side^2 && break
        stopped = _foreach_crossing_frame(
            lattice, window_side, frame_evaluated, max_frame_evaluations,
        ) do frame
            atom_count <= length(frame.allowed) || return nothing
            offsets = (-Int(minimum(completion))):(atom_count - Int(maximum(completion)))
            isempty(offsets) && return nothing
            evaluated == max_evaluations && return :budget
            context = _prepare_sat_frame(lattice, frame)
            for offset in offsets
                evaluated == max_evaluations && return :budget
                solver, selected = _fixed_crossing_sat_problem(
                    target_reduced, context, atom_count, offset,
                )
                while true
                    evaluated == max_evaluations && return :budget
                    evaluated += 1
                    analysis = _solve_next_fixed_crossing_sat!(
                        solver, selected, target_reduced, lattice, context,
                    )
                    analysis === nothing && break
                    valid, verified_offset = is_gadget_replacement(
                        target_graph, analysis.graph, target_boundary,
                        analysis.boundary,
                    )
                    valid || error("SAT candidate failed the fixed verifier")
                    analysis.offset == verified_offset ||
                        error("SAT candidate offset mismatch")
                    key = _canonical_crossing_frame_key(lattice, (
                        allowed=analysis.patch.coordinates,
                        pins=analysis.patch.pins,
                        rays=analysis.patch.rays,
                    ))
                    key in seen_gadgets && continue
                    push!(seen_gadgets, key)
                    push!(
                        gadgets,
                        _unweighted_gadget(target_graph, lattice, analysis),
                    )
                    length(gadgets) == max_results && return :solution
                end
            end
            return nothing
        end
        stopped === nothing || return gadgets, evaluated, stopped
    end
    return gadgets, evaluated, :search_space_exhausted
end

function _foreach_crossing_frame(visit, lattice, side, evaluated, limit)
    window = _LatticeCoordinate[
        (column, row) for column in 0:side-1 for row in 0:side-1
    ]
    ray_options = eachindex(_lattice_directions(lattice))
    seen = Set{Tuple}()
    for pin_set in combinations(window, 4)
        for pins in permutations(pin_set)
            for rays in Iterators.product(ntuple(_ -> ray_options, 4)...)
                evaluated[] == limit && return :frame_budget
                evaluated[] += 1
                ordered_pins = collect(pins)
                ray_indices = collect(rays)
                patch = _LatticePatch(ordered_pins, ordered_pins, ray_indices)
                all(_check_crossing_frame(lattice, patch)) || continue
                allowed = [site for site in window if site in ordered_pins || all(
                    _check_crossing_frame(
                        lattice,
                        _LatticePatch(
                            [ordered_pins; site], ordered_pins, ray_indices,
                        ),
                    ),
                )]
                frame = (
                    pins=ordered_pins, rays=ray_indices, allowed=allowed,
                )
                key = _canonical_crossing_frame_key(lattice, frame)
                key in seen && continue
                push!(seen, key)
                result = visit(frame)
                result === nothing || return result
            end
        end
    end
    return nothing
end
