using GadgetSearch
using Combinatorics
using Graphs
using Test

function cross_graph()
    graph = SimpleGraph(4)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 2, 4)
    return graph
end

function cross_edge_graph()
    graph = SimpleGraph(6)
    for edge in ((1,2), (1,5), (2,6), (3,5), (4,6))
        add_edge!(graph, edge...)
    end
    return graph
end
@testset "Unweighted search" begin
    @testset "selected connectivity matches enumeration" begin
        lattice_host, _, _ = GadgetSearch._materialize_lattice_patch(
            Triangular(), GadgetSearch._LatticePatch(
                [(0,0),(1,0),(0,1),(1,1),(2,1)], [(0,0)], [1],
            ),
        )
        for host in (path_graph(5), cycle_graph(5), lattice_host)
            adjacent = [collect(neighbors(host, vertex)) for vertex in vertices(host)]
            for mask in 0:(1 << nv(host))-1
                !iszero(mask & 1) || continue
                chosen = [vertex for vertex in vertices(host)
                    if !iszero(mask & (1 << (vertex - 1)))]
                cnf = GadgetSearch._SatCnf()
                selected = [GadgetSearch._sat_variable!(cnf) for _ in vertices(host)]
                GadgetSearch._add_selected_connectivity!(
                    cnf, selected, adjacent, 1, length(chosen),
                )
                for vertex in vertices(host)
                    GadgetSearch._sat_clause!(
                        cnf, vertex in chosen ? selected[vertex] : -selected[vertex],
                    )
                end
                actual = GadgetSearch._next_selected_assignment!(
                    GadgetSearch._new_sat_solver(cnf), selected,
                ) !== nothing
                @test actual == is_connected(induced_subgraph(host, chosen)[1])
            end
        end
    end

    @testset "fixed verifier" begin
        reduced = calculate_reduced_alpha_tensor(cross_graph(), [1, 2, 3, 4])
        @test GadgetSearch.inf_mask(reduced) == BigInt(60576)
        @test is_diff_by_constant(reduced .+ 3, reduced) == (true, 3.0)
    end

    @testset "crossing frame" begin
        square = [(0, 0), (-1, 0), (0, 1), (1, 0), (0, -1)]
        pins = square[2:5]
        rays = [(-1, 0), (0, 1), (1, 0), (0, -1)]
        @test all(check_crossing_frame(Square(), square, pins, rays))
        @test !check_crossing_frame(Square(), [square; (-2, 1)], pins, rays).G4
        @test_throws ArgumentError check_crossing_frame(
            Square(), [(0, 0)], pins, rays,
        )
        @test_throws ArgumentError check_crossing_frame(
            Square(), [square; square[1]], pins, rays,
        )
        @test_throws ArgumentError check_crossing_frame(
            Square(), square, [pins[1], pins[1], pins[3], pins[4]], rays,
        )
        triangular = [(0, 0), (-1, 0), (0, 1), (1, 0), (-1, -1)]
        triangular_pins = triangular[2:5]
        @test all(check_crossing_frame(
            Triangular(), triangular, triangular_pins, rays,
        ))
    end

    @testset "pruned frame generation preserves geometric candidates" begin
        lattice = Triangular()
        pins = [(2,0),(2,3),(0,1),(1,0)]
        ray_indices = eachindex(GadgetSearch._lattice_directions(lattice))
        brute_force = Set{Tuple}()
        for order in permutations(1:4)
            ordered_pins = pins[collect(order)]
            for rays in Iterators.product(ntuple(_ -> ray_indices, 4)...)
                all(GadgetSearch._check_crossing_frame(
                    lattice,
                    GadgetSearch._LatticePatch(
                        ordered_pins, ordered_pins, collect(rays),
                    ),
                )) || continue
                push!(brute_force, (Tuple(ordered_pins), rays))
            end
        end

        pruned = Set{Tuple}()
        ray_options = [
            [ray for ray in ray_indices if GadgetSearch._single_pin_corridor_clear(
                lattice, pins, index, ray,
            )] for index in eachindex(pins)
        ]
        for rays_tuple in Iterators.product(ray_options...)
            GadgetSearch._pin_rays_are_compatible(
                lattice, pins, rays_tuple,
            ) || continue
            physical_rays = collect(rays_tuple)
            checks = GadgetSearch._check_crossing_frame(
                lattice, GadgetSearch._LatticePatch(pins, pins, physical_rays),
            )
            checks[1] && checks[3] && checks[4] || continue
            for order_tuple in permutations(1:4)
                order = collect(order_tuple)
                ordered_pins = pins[order]
                ordered_rays = physical_rays[order]
                GadgetSearch._labels_alternate(
                    lattice, ordered_pins, ordered_rays,
                ) || continue
                push!(pruned, (Tuple(ordered_pins), Tuple(ordered_rays)))
            end
        end
        @test pruned == brute_force

        clique_frames = Set{Tuple}()
        expected_prefixed_frames = Set{Tuple}()
        clique_evaluated = Ref(0)
        candidates = GadgetSearch._crossing_port_candidates(
            lattice, [(column, row) for column in 0:2 for row in 0:2],
        )
        GadgetSearch._foreach_crossing_frame_clique(
            lattice, (3,3), clique_evaluated, typemax(Int); min_allowed=4,
        ) do frame, _, _
            key = GadgetSearch._canonical_crossing_frame_key(lattice, frame)
            push!(clique_frames, key)
            indices = sort([findfirst(==(port), candidates)
                for port in zip(frame.pins, frame.rays)])
            indices[1:2] == [1, 9] && push!(expected_prefixed_frames, key)
            return nothing
        end
        sharded_frames = Set{Tuple}()
        for shard_index in 0:2
            evaluated = Ref(0)
            GadgetSearch._foreach_crossing_frame_clique(
                lattice, (3,3), evaluated, typemax(Int);
                min_allowed=4, shard_index, shard_count=3,
            ) do frame, _, _
                push!(sharded_frames,
                    GadgetSearch._canonical_crossing_frame_key(lattice, frame))
                return nothing
            end
        end
        @test sharded_frames == clique_frames

        prefixed_frames = Set{Tuple}()
        prefixed_evaluated = Ref(0)
        GadgetSearch._foreach_crossing_frame_clique(
            lattice, (3,3), prefixed_evaluated, typemax(Int);
            min_allowed=4, port_prefix=(1, 9),
        ) do frame, _, _
            push!(prefixed_frames,
                GadgetSearch._canonical_crossing_frame_key(lattice, frame))
            return nothing
        end
        @test prefixed_frames == expected_prefixed_frames
        cross23_window = [(column, row)
            for column in 0:7 for row in 0:7]
        @test GadgetSearch._crossing_port_prefix(
            Triangular(), cross23_window, 5117,
        ) == (21, 177)

        for (candidate_pins, candidate_rays, side) in (
            (pins, [6,2,4,5], 4),
            ([(5,1),(0,3),(3,5),(6,6)], [6,3,3,2], 8),
        )
            window = [(column, row)
                for column in 0:side-1 for row in 0:side-1]
            expected = sort([site for site in window
                if site in candidate_pins || begin
                    checks = GadgetSearch._check_crossing_frame(
                        lattice, GadgetSearch._LatticePatch(
                            [candidate_pins; site], candidate_pins, candidate_rays,
                        ),
                    )
                    checks[1] && checks[3] && checks[4]
                end])
            @test GadgetSearch._allowed_crossing_sites(
                lattice, window, candidate_pins, candidate_rays, 0,
            ) == expected
        end

    end

    @testset "fixed SAT positive control" begin
        coordinates = [
            (0,3),(6,6),(1,3),(2,1),(2,2),(2,3),(3,1),(3,3),
            (3,4),(3,5),(4,1),(4,2),(4,4),(4,5),(5,1),(5,2),
            (5,3),(5,4),(5,5),(6,2),(6,3),(6,5),(7,4),
        ]
        frame = (
            pins=[(5,1),(0,3),(3,5),(6,6)],
            rays=[6,3,3,2],
            allowed=coordinates,
        )
        target_reduced = vec(calculate_reduced_alpha_tensor(
            cross_graph(), [1,2,3,4],
        ))
        analysis = GadgetSearch._solve_fixed_crossing_sat(
            target_reduced, Triangular(), frame, 23, 7,
        )
        @test analysis !== nothing
        @test analysis.offset == 7
        context = GadgetSearch._prepare_sat_frame(Triangular(), frame)
        filter_solver, filter_selected = GadgetSearch._lower_state_filter_problem(
            target_reduced, context, 23, 7,
            first(GadgetSearch._essential_lower_states(target_reduced)),
        )
        @test GadgetSearch._next_selected_assignment!(
            filter_solver, filter_selected,
        ) !== nothing
        @test nv(analysis.graph) == 23
        @test is_gadget_replacement(
            cross_graph(), analysis.graph, [1,2,3,4], analysis.boundary,
        ) == (true, 7.0)
        @test all(check_crossing_frame(
            Triangular(), analysis.patch.coordinates, analysis.patch.pins,
            GadgetSearch._patch_ray_directions(Triangular(), analysis.patch),
        ))

        extended_sites = [coordinates; (6,0); (3,6)]
        extended_pins = [(6,0), (0,3), (3,6), (6,6)]
        extended_analysis = GadgetSearch._analyze_crossing_candidate(
            target_reduced, Triangular(), extended_sites, extended_pins,
            frame.rays,
        )
        @test extended_analysis.solved
        @test extended_analysis.offset == 8
        extended_gadget = GadgetSearch._unweighted_gadget(
            cross_graph(), Triangular(), extended_analysis,
        )
        optimized = optimize_unweighted_gadget(
            extended_gadget, [1,2,3,4];
            min_vertices=23, max_sat_evaluations=0,
        )
        @test nv(optimized.gadget.replacement_graph) == 23
        @test optimized.gadget.constant_offset == 7
        @test only(optimized.steps).rule == :opposite_leaf_pin_contraction
        @test nv(only(optimized.steps).before.replacement_graph) == 25
        @test nv(only(optimized.steps).after.replacement_graph) == 23
        @test optimized.termination_reason == :minimum_vertices
        @test is_gadget_replacement(
            cross_graph(), optimized.gadget.replacement_graph, [1,2,3,4],
            optimized.gadget.boundary_vertices,
        ) == (true, 7.0)

        directions = GadgetSearch._lattice_directions(Triangular())
        tail_middle = GadgetSearch._lattice_step(
            Triangular(), frame.pins[1], directions[frame.rays[1]], 1,
        )
        tail_pin = GadgetSearch._lattice_step(
            Triangular(), frame.pins[1], directions[frame.rays[1]], 2,
        )
        tail_analysis = GadgetSearch._analyze_crossing_candidate(
            target_reduced, Triangular(),
            [coordinates; tail_middle; tail_pin],
            [tail_pin; frame.pins[2:4]], frame.rays,
        )
        @test tail_analysis.solved
        tail_gadget = GadgetSearch._unweighted_gadget(
            cross_graph(), Triangular(), tail_analysis,
        )
        tail_optimized = optimize_unweighted_gadget(
            tail_gadget, [1,2,3,4];
            min_vertices=23, max_sat_evaluations=0,
        )
        @test only(tail_optimized.steps).rule ==
            :even_boundary_tail_contraction
        @test nv(tail_optimized.gadget.replacement_graph) == 23

        rewrite_start_coordinates = [
            (4,7),(2,0),(6,0),(7,3),(1,3),(2,4),(3,6),(3,7),
            (2,2),(3,4),(3,5),(2,1),(3,3),(4,4),(4,5),(3,0),
            (4,2),(5,4),(5,5),(4,0),(4,1),(5,2),(5,3),(6,5),
            (5,1),(6,2),(6,3),(7,4),
        ]
        rewrite_start_analysis = GadgetSearch._analyze_crossing_candidate(
            target_reduced, Triangular(), rewrite_start_coordinates,
            [(4,7),(2,0),(6,0),(7,3)], [2,5,5,1],
        )
        @test rewrite_start_analysis.solved
        @test rewrite_start_analysis.offset == 10
        rewrite_start = GadgetSearch._unweighted_gadget(
            cross_graph(), Triangular(), rewrite_start_analysis,
        )
        rewritten = optimize_unweighted_gadget(
            rewrite_start, [1,2,3,4];
            min_vertices=23, max_sat_evaluations=8, host_radius=1,
        )
        @test rewritten.termination_reason == :minimum_vertices
        @test rewritten.sat_evaluations <= 8
        @test [step.rule for step in rewritten.steps] == [
            :even_boundary_tail_contraction,
            :frame_rewrite_resynthesis,
            :opposite_leaf_pin_contraction,
        ]
        @test [
            nv(step.before.replacement_graph) => nv(step.after.replacement_graph)
            for step in rewritten.steps
        ] == [28 => 26, 26 => 25, 25 => 23]
        @test is_gadget_replacement(
            cross_graph(), rewritten.gadget.replacement_graph, [1,2,3,4],
            rewritten.gadget.boundary_vertices,
        ) == (true, 7.0)

        rewrite_budget = optimize_unweighted_gadget(
            GadgetSearch._unweighted_gadget(
                cross_graph(), Triangular(), analysis,
            ),
            [1,2,3,4];
            min_vertices=22, max_sat_evaluations=1, host_radius=0,
        )
        @test isempty(rewrite_budget.steps)
        @test rewrite_budget.sat_evaluations == 1
        @test rewrite_budget.termination_reason == :sat_budget
        if !Sys.iswindows()
            canonical_coordinates = [
                (0,-5),(-1,-4),(-1,-3),(1,-5),(0,-3),(3,-7),(2,-6),
                (2,-5),(1,-4),(1,-3),(0,-2),(0,-1),(3,-6),(2,-4),
                (1,-2),(0,0),(5,-7),(4,-6),(3,-4),(3,-3),(2,-2),
                (6,-7),(4,-4),(4,-3),(6,-6),(6,-5),(5,-4),(4,-2),
            ]
            joint_frame = [
                ((0,0),(1,0)), ((-1,-3),(-1,1)),
                ((3,-7),(0,-1)), ((4,-2),(1,0)),
            ]
            joint_solution = GadgetSearch._from_canonical.(
                Ref(Triangular()), canonical_coordinates,
            )
            joint_pins = GadgetSearch._from_canonical.(
                Ref(Triangular()), first.(joint_frame),
            )
            directions = GadgetSearch._lattice_directions(Triangular())
            joint_rays = [findfirst(==(ray), directions) for ray in last.(joint_frame)]
            joint_analysis = GadgetSearch._analyze_crossing_candidate(
                target_reduced, Triangular(), joint_solution, joint_pins, joint_rays,
            )
            @test joint_analysis.solved
            @test joint_analysis.offset == 10
            joint_cnf, joint_selected, joint_choices, joint_coordinates =
                GadgetSearch._joint_crossing_sat_cnf(
                    target_reduced, Triangular(), (8,8), 28, 10,
                )
            selected_coordinates = Set(joint_solution)
            for (vertex, coordinate) in enumerate(joint_coordinates)
                GadgetSearch._sat_clause!(
                    joint_cnf, coordinate in selected_coordinates ?
                    joint_selected[vertex] : -joint_selected[vertex],
                )
            end
            for label in 1:4,
                (vertex, direction, variable) in joint_choices[label]
                chosen = joint_coordinates[vertex] == joint_pins[label] &&
                    direction == joint_rays[label]
                GadgetSearch._sat_clause!(
                    joint_cnf, chosen ? variable : -variable,
                )
            end
            status, _ = GadgetSearch._next_selected_assignment_limited!(
                GadgetSearch._new_sat_solver(joint_cnf), joint_selected, 100_000,
            )
            @test status == :sat

            cross23_sites = [
                (0,0),(3,-7),(0,-1),(-2,-1),(-1,-1),(0,-2),(-2,-2),
                (0,-3),(1,-3),(2,-4),(-2,-3),(-1,-3),(1,-4),(2,-5),
                (-2,-4),(-1,-4),(0,-5),(1,-5),(2,-6),(-1,-5),(0,-6),
                (2,-7),(1,-7),
            ]
            cross23_pins = [(0,0),(-2,-4),(3,-7),(2,-4)]
            cross23_rays = [1,4,6,1]
            cross23_solution = GadgetSearch._from_canonical.(
                Ref(Triangular()), cross23_sites,
            )
            cross23_boundary = GadgetSearch._from_canonical.(
                Ref(Triangular()), cross23_pins,
            )
            cross23_cnf, cross23_selected, cross23_choices, cross23_coordinates =
                GadgetSearch._joint_crossing_sat_cnf(
                    target_reduced, Triangular(), (8,8), 23, 7;
                    canonical_shift=(-1,0),
                )
            selected_coordinates = Set(cross23_solution)
            for (vertex, coordinate) in enumerate(cross23_coordinates)
                GadgetSearch._sat_clause!(
                    cross23_cnf, coordinate in selected_coordinates ?
                    cross23_selected[vertex] : -cross23_selected[vertex],
                )
            end
            for label in 1:4,
                (vertex, direction, variable) in cross23_choices[label]
                chosen = cross23_coordinates[vertex] == cross23_boundary[label] &&
                    direction == cross23_rays[label]
                GadgetSearch._sat_clause!(
                    cross23_cnf, chosen ? variable : -variable,
                )
            end
            status, _ = GadgetSearch._next_selected_assignment_limited!(
                GadgetSearch._new_sat_solver(cross23_cnf), cross23_selected,
                100_000,
            )
            @test status == :sat
            cross23_analysis = GadgetSearch._analyze_crossing_candidate(
                target_reduced, Triangular(), cross23_solution,
                cross23_boundary, cross23_rays,
            )
            @test cross23_analysis.solved
            @test cross23_analysis.offset == 7
        end
        shifted_target = map(value -> isfinite(value) ? value + 8 : value, target_reduced)
        negative_offset = GadgetSearch._solve_fixed_crossing_sat(
            shifted_target, Triangular(), frame, 23, -1,
        )
        @test negative_offset !== nothing
        @test negative_offset.offset == -1

        ksg_coordinates = [(0, 0), (-1, 0), (0, 1), (1, 0), (0, -1)]
        ksg_frame = (
            pins=ksg_coordinates[2:5],
            rays=[4, 7, 5, 2],
            allowed=ksg_coordinates,
        )
        ksg_patch = GadgetSearch._LatticePatch(
            ksg_coordinates, ksg_frame.pins, ksg_frame.rays,
        )
        ksg_graph, ksg_boundary, _ = GadgetSearch._materialize_lattice_patch(
            Square(), ksg_patch,
        )
        ksg_target = vec(calculate_reduced_alpha_tensor(ksg_graph, ksg_boundary))
        ksg_analysis = GadgetSearch._solve_fixed_crossing_sat(
            ksg_target, Square(), ksg_frame, 5, 0,
        )
        @test ksg_analysis !== nothing
        @test nv(ksg_analysis.graph) == 5
        @test all(check_crossing_frame(
            Square(), ksg_analysis.patch.coordinates, ksg_analysis.patch.pins,
            GadgetSearch._patch_ray_directions(Square(), ksg_analysis.patch),
        ))
        @test isnothing(GadgetSearch._solve_fixed_crossing_sat(
            ksg_target, Square(), ksg_frame, 6, 0,
        ))
    end

    @testset "joint frame and occupancy SAT" begin
        target = complete_graph(4)
        target_reduced = vec(calculate_reduced_alpha_tensor(
            target, collect(1:4),
        ))
        public_joint = search_unweighted_gadget_joint(
            target, collect(1:4), Square();
            window_shape=(2,2), atom_count=4, offset=0, seconds=30,
        )
        @test public_joint !== nothing
        @test nv(public_joint.replacement_graph) == 4
        @test first(is_gadget_replacement(
            target, public_joint.replacement_graph, collect(1:4),
            public_joint.boundary_vertices,
        ))
        @test all(check_crossing_frame(
            Square(), public_joint.lattice_coordinates,
            public_joint.lattice_coordinates[public_joint.boundary_vertices],
            public_joint.pin_rays,
        ))
        cross_edge_target = cross_edge_graph()
        cross_edge = search_unweighted_gadget_joint(
            cross_edge_target, collect(1:4), Triangular();
            window_shape=(3,4), atom_count=9, offset=1, seconds=30,
            canonical_shift=(-1,1), first_ray=1,
        )
        @test cross_edge !== nothing
        @test first(is_gadget_replacement(
            cross_edge_target, cross_edge.replacement_graph, collect(1:4),
            cross_edge.boundary_vertices,
        ))
        optimized_cross_edge = optimize_unweighted_gadget(
            cross_edge, collect(1:4);
            min_vertices=4, max_sat_evaluations=256, host_radius=1,
        )
        @test nv(optimized_cross_edge.gadget.replacement_graph) == 9
        @test isempty(optimized_cross_edge.steps)
        @test optimized_cross_edge.sat_evaluations == 256
        @test optimized_cross_edge.termination_reason == :sat_budget
        rotated_cross_edge = search_unweighted_gadget_joint(
            cross_edge_target, collect(1:4), Triangular();
            window_shape=(3,4), atom_count=9, offset=1, seconds=30,
            canonical_shift=(0,0), first_ray=2,
        )
        @test rotated_cross_edge !== nothing
        @test first(is_gadget_replacement(
            cross_edge_target, rotated_cross_edge.replacement_graph,
            collect(1:4), rotated_cross_edge.boundary_vertices,
        ))
        @test isnothing(GadgetSearch._solve_joint_crossing_sat(
            target_reduced, Square(), (2,2), 4, 1,
        ))

        crossing_points = [(0,0), (0,2), (2,0), (2,2)]
        for chosen in Iterators.product(ntuple(_ -> 1:4, 4)...)
            points = crossing_points[collect(chosen)]
            hull = GadgetSearch._strict_convex_hull(points)
            expected = length(unique(points)) == 4 && length(hull) == 4 &&
                GadgetSearch._interfaces_alternate(hull, points)
            actual = GadgetSearch._interfaces_form_alternating_quadrilateral(
                points,
            )
            @test actual == expected
        end
    end

    @testset "SAT model enumeration" begin
        cnf = GadgetSearch._SatCnf()
        foreach(_ -> GadgetSearch._sat_variable!(cnf), 1:3)
        solver = GadgetSearch._new_sat_solver(cnf)
        assignments = Set{Tuple{Bool, Bool}}()
        while true
            assignment = GadgetSearch._next_selected_assignment!(solver, [1, 2])
            assignment === nothing && break
            push!(assignments, Tuple(assignment))
        end
        @test assignments == Set([
            (false, false), (false, true), (true, false), (true, true),
        ])
        if !Sys.iswindows()
            sat_cnf = GadgetSearch._SatCnf()
            variable = GadgetSearch._sat_variable!(sat_cnf)
            GadgetSearch._sat_clause!(sat_cnf, variable)
            status, assignment = GadgetSearch._next_selected_assignment_limited!(
                GadgetSearch._new_sat_solver(sat_cnf), [variable], 10,
            )
            @test status == :sat
            @test assignment == [true]

            unsat_cnf = GadgetSearch._SatCnf()
            variable = GadgetSearch._sat_variable!(unsat_cnf)
            GadgetSearch._sat_clause!(unsat_cnf, variable)
            GadgetSearch._sat_clause!(unsat_cnf, -variable)
            status, assignment = GadgetSearch._next_selected_assignment_limited!(
                GadgetSearch._new_sat_solver(unsat_cnf), [variable], 10,
            )
            @test status == :unsat
            @test assignment === nothing
        end
    end

    @testset "public bounded search" begin
        @test GadgetSearch._crossing_window_shapes(8, 23) == [
            (8,8), (8,7), (7,8), (8,6), (6,8),
            (8,5), (5,8), (8,4), (4,8), (8,3), (3,8),
        ]
        @test_throws ArgumentError search_unweighted_gadgets(
            path_graph(3), [1,2,3], Triangular(),
        )
        @test_throws ArgumentError search_unweighted_gadgets(
            path_graph(4), [1,2,3,4]; window_side=1,
        )
        @test_throws ArgumentError search_unweighted_gadgets(
            path_graph(4), [1,2,3,4]; checkpoint_interval=0,
        )
        ksg = search_unweighted_gadgets(
            cross_graph(), [1,2,3,4], Square();
            min_vertices=4, max_vertices=4, max_evaluations=1,
        )
        @test ksg.lattice == :KSG
        @test ksg.evaluated <= 1
        @test ksg.termination_reason in (:budget, :frame_budget)

        exhausted = search_unweighted_gadgets(
            cross_graph(), [1,2,3,4], Triangular();
            min_vertices=4, max_vertices=4, max_evaluations=1,
        )
        @test exhausted.evaluated == 1
        @test isempty(exhausted.gadgets)
        @test exhausted.termination_reason in (:budget, :frame_budget)

        ksg_target = complete_graph(4)
        solved = search_unweighted_gadgets(
            ksg_target, collect(1:4), Square();
            min_vertices=4, max_vertices=4, max_evaluations=50, max_results=1,
            max_frame_evaluations=100_000, window_side=2,
        )
        @test solved.termination_reason == :solution
        @test length(solved.gadgets) == 1
        @test is_gadget_replacement(
            ksg_target, solved.gadgets[1].replacement_graph,
            collect(1:4), solved.gadgets[1].boundary_vertices,
        )[1]

        fully_enumerated = search_unweighted_gadgets(
            ksg_target, collect(1:4), Square();
            min_vertices=4, max_vertices=4, max_evaluations=50, max_results=2,
            max_frame_evaluations=100_000, window_side=2,
        )
        @test fully_enumerated.termination_reason == :search_space_exhausted
        @test length(fully_enumerated.gadgets) == 1
        @test fully_enumerated.evaluated > solved.evaluated

        mktempdir() do directory
            checkpoint_path = joinpath(directory, "search.checkpoint")
            interrupted = search_unweighted_gadgets(
                ksg_target, collect(1:4), Square();
                min_vertices=4, max_vertices=4, max_evaluations=3,
                max_results=2, max_frame_evaluations=100_000,
                window_side=2, checkpoint_path, checkpoint_interval=1,
            )
            @test interrupted.termination_reason == :budget
            @test isfile(checkpoint_path)
            checkpoint = read_unweighted_search_checkpoint(checkpoint_path)
            @test checkpoint.window_shape == (2, 2)
            @test checkpoint.frame_cursor >= 1
            @test checkpoint.order_cursor >= 1
            @test checkpoint.offset_cursor >= 0
            resumed = search_unweighted_gadgets(
                ksg_target, collect(1:4), Square();
                min_vertices=4, max_vertices=4, max_evaluations=20,
                max_results=2, max_frame_evaluations=100_000,
                window_side=2, checkpoint_path, checkpoint_interval=1,
            )
            @test resumed.termination_reason == :search_space_exhausted
            @test length(resumed.gadgets) == 1
        end
        point = (3, 4)
        @test GadgetSearch._from_canonical(Square(), point) == point
        canonical = GadgetSearch._canonical_coordinate(Triangular(), point)
        @test GadgetSearch._from_canonical(Triangular(), canonical) == point
        @test GadgetSearch._lattice_distance(Square(), (0, 0), (2, -1)) == 2
        @test GadgetSearch._lattice_distance(Triangular(), (0, 0), (2, -1)) == 2
    end
end
