# Crossing Gadget — Pin Selection Walkthrough
#
# A crossing gadget has 4 pins representing two independent wires:
#   Wire 1: pin1 → pin3 (copy), must be on opposite boundary sides
#   Wire 2: pin2 → pin4 (copy), must be on opposite boundary sides

using GadgetSearch
using Graphs

outdir = joinpath(pkgdir(GadgetSearch), "examples", "walkthrough_output")
mkpath(outdir)

# %% Step 1: crossing constraint truth table

crossing_constraint = TruthTableConstraint(Bool[
    0 0 0 0
    0 1 0 1
    1 0 1 0
    1 1 1 1
])

println("Truth table:")
display(crossing_constraint.truth_table)

# %% Step 2: a cross-shaped subgraph on KingsGraph

lattice = Triangular()
positions = Tuple{Int,Int}[(1,1), (2,1), (3,1), (1,2), (2,2), (3,2), (2,3)]

physical = GadgetSearch.get_physical_positions(lattice, positions)
g = GadgetSearch.unit_disk_graph(physical, get_radius(lattice))
println("\n$(nv(g)) vertices, $(ne(g)) edges")

# %% Step 3: plot open vertices

plot_open_vertices(positions, lattice, joinpath(outdir, "step3_open_vertices.png"))

# %% Step 4: boundary and direction classification

boundary = GadgetSearch.classify_boundary(positions)
open_boundary = GadgetSearch._get_open_boundary(positions, lattice)

println("\nBoundary (★ = open):")
for i in sort(collect(keys(boundary)))
    marker = haskey(open_boundary, i) ? " ★" : ""
    println("  v$i $(positions[i]) → $(boundary[i])$marker")
end

dir_verts = Dict{Symbol, Vector{Int}}()
for (v, dirs) in open_boundary
    for d in dirs
        push!(get!(dir_verts, d, Int[]), v)
    end
end
println("\nDirection groups:")
for d in [:up, :down, :left, :right]
    vs = sort(get(dir_verts, d, Int[]))
    isempty(vs) && continue
    println("  $d → [$(join(["v$v$(positions[v])" for v in vs], ", "))]")
end

# %% Step 5: pin candidates — with and without opposite constraint

opposite_pairs = [(1, 3), (2, 4)]

candidates_all = compute_open_pin_candidates(positions, "KSG", 4)
candidates = compute_open_pin_candidates(positions, "KSG", 4; opposite_pairs=opposite_pairs)

println("\nWithout constraint: $(length(candidates_all)) combinations")
println("With opposite_pairs=$opposite_pairs: $(length(candidates)) combinations\n")

for (ci, c) in enumerate(candidates)
    ds = [open_boundary[c[j]] for j in 1:4]
    pair1 = "pin1=v$(c[1])$(ds[1]) ↔ pin3=v$(c[3])$(ds[3])"
    pair2 = "pin2=v$(c[2])$(ds[2]) ↔ pin4=v$(c[4])$(ds[4])"
    println("  #$ci  $pair1,  $pair2")
end

# %% Step 6: plot pin assignments

for (ci, c) in enumerate(candidates)
    outpath = joinpath(outdir, "step6_pins_$ci.png")
    plot_pin_assignment(positions, lattice, c, outpath; opposite_pairs=opposite_pairs)
end
