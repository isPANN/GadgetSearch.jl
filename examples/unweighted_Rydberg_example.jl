# Run from the repository root:
# julia --project=. -e 'using Pkg; Pkg.instantiate()'
# julia --project=. examples/unweighted_Rydberg_example.jl

using GadgetSearch, Graphs, JSON3

# 1. Define both targets. Pins are ordered P1, P2, P3, P4.
# CROSS+EDGE has internal vertices 5 and 6. Every atom has unit weight.
cross_edge = SimpleGraph(6)
for edge in ((1, 2), (1, 5), (2, 6), (3, 5), (4, 6))
    add_edge!(cross_edge, edge...)
end
cross = SimpleGraph(4)
add_edge!(cross, 1, 3)
add_edge!(cross, 2, 4)
target_pins = [1, 2, 3, 4]
lattice = Triangular()

# 2. Set one blank-window search instance per target.
# Require reduced_alpha(replacement) = reduced_alpha(target) + offset.
# These settings admit known realizations but supply no atom placement.
# CROSS is a harder search, so it gets a larger SAT-solving budget.
searches = [
    ("CROSS+EDGE", "CROSS_EDGE", cross_edge, (;
        window_shape=(3, 4), atom_count=9, offset=1,
        seconds=30, canonical_shift=(-1, 1), first_ray=1,
    )),
    ("CROSS", "CROSS", cross, (;
        window_shape=(8, 6), atom_count=23, offset=7,
        seconds=120, canonical_shift=(-6, 3), first_ray=2,
    )),
]

for (label, name, target, settings) in searches
    println("\nSearching $label: $(settings.atom_count) atoms, offset $(settings.offset)")
    gadget = search_unweighted_gadgets(
        target, target_pins, lattice; settings..., seed=1,
    )
    isnothing(gadget) && error(
        "$label: no gadget returned; this instance was unsatisfiable or timed out. " *
        "Increase seconds or try another instance.",
    )

    # 3. Recompute the tensors and check all four crossing-frame conditions.
    @assert nv(gadget.replacement_graph) == settings.atom_count
    @assert is_connected(gadget.replacement_graph)
    @assert is_gadget_replacement(
        target, gadget.replacement_graph, target_pins, gadget.boundary_vertices,
    ) == (true, Float64(settings.offset))
    geometry = check_gadget_geometry(
        lattice, gadget.lattice_coordinates,
        gadget.lattice_coordinates[gadget.boundary_vertices], gadget.pin_rays,
    )
    @assert all(geometry)

    println("Atoms: ", nv(gadget.replacement_graph))
    println("Pins P1..P4: ", gadget.boundary_vertices)
    println("Reduced-alpha offset: ", gadget.constant_offset)
    println("Geometry: ", geometry)
    for (vertex, coordinate) in enumerate(gadget.lattice_coordinates)
        println("Vertex $vertex: lattice $coordinate, position $(gadget.pos[vertex])")
    end

    # 4. Save separate embeddings and plots for CROSS_EDGE and CROSS.
    output = joinpath(@__DIR__, "unweighted_Rydberg_$name")
    open(output * ".json", "w") do io
        JSON3.write(io, (;
            lattice=gadget.lattice,
            edges=[(src(edge), dst(edge)) for edge in edges(gadget.replacement_graph)],
            pins=gadget.boundary_vertices, pin_rays=gadget.pin_rays,
            lattice_coordinates=gadget.lattice_coordinates, pos=gadget.pos,
            constant_offset=gadget.constant_offset,
        ))
    end
    GadgetSearch.plot_graph(gadget.replacement_graph, output * ".png"; pos=gadget.pos)
    println("Embedding saved as ", output * ".json")
end
