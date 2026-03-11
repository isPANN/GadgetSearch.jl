abstract type LatticeType end
struct Square <: LatticeType end
struct Triangular <: LatticeType end

get_radius(::Square) = 1.5
get_radius(::Triangular) = 1.1

get_physical_positions(::Square, pos::Vector{Tuple{Int, Int}}) = Vector{Tuple{Float64, Float64}}(pos)
function get_physical_positions(::Triangular, pos::Vector{Tuple{Int, Int}})
    h = sqrt(3) / 2
    return Vector{Tuple{Float64, Float64}}([(x + (y - 1) * 0.5, y * h) for (x, y) in pos])
end

"""
    get_pin_positions(lattice::LatticeType, nx::Int, ny::Int)

Return candidate pin positions on the four boundary sides of an `(nx+2)×(ny+2)`
padded grid.

Coordinates use **matrix convention** `(row, col)` where row 1 is the top and
column 1 is the left.  The physical mapping to plot coordinates is handled
separately by `get_physical_positions`.

# Returns
`(top, bottom, left, right)` — each a vector of physical `(Float64, Float64)`
positions.
"""
function get_pin_positions(lattice::LatticeType, nx::Int, ny::Int)
    top_candidates    = get_physical_positions(lattice, Tuple{Int, Int}[(1, y)     for y in 2:ny+1])
    bottom_candidates = get_physical_positions(lattice, Tuple{Int, Int}[(nx+2, y)  for y in 2:ny+1])
    left_candidates   = get_physical_positions(lattice, Tuple{Int, Int}[(x, 1)     for x in 2:nx+1])
    right_candidates  = get_physical_positions(lattice, Tuple{Int, Int}[(x, ny+2)  for x in 2:nx+1])
    return top_candidates, bottom_candidates, left_candidates, right_candidates
end

function get_inner_points(lattice::LatticeType, nx::Int, ny::Int)
    original_points = Tuple{Int, Int}[(x, y) for x in 2:nx+1, y in 2:ny+1]
    return get_physical_positions(lattice, vec(original_points))
end

"""
    triangular_adjacency(i1, j1, i2, j2) -> Bool

Return `true` if integer grid positions `(i1, j1)` and `(i2, j2)` are
nearest neighbours on the triangular lattice.

The triangular lattice is a rectangular grid with **one diagonal** added per
cell (the `(+1,−1)` / `(−1,+1)` direction). Under the parallelogram physical
mapping `(x, y) → (x + (y−1)·0.5,  y·√3/2)` all six resulting neighbours lie
at Euclidean distance exactly 1.0.

```
neighbours of (i,j):
  (i±1, j)   – horizontal
  (i, j±1)   – vertical
  (i+1, j−1) – diagonal ↘
  (i−1, j+1) – diagonal ↖
```

No floating-point computation is required.
"""
function triangular_adjacency(i1::Int, j1::Int, i2::Int, j2::Int)
    di, dj = i2 - i1, j2 - j1
    return (abs(di) == 1 && dj == 0)  ||  # horizontal
           (di == 0 && abs(dj) == 1)  ||  # vertical
           (di ==  1 && dj == -1)     ||  # diagonal ↘
           (di == -1 && dj ==  1)         # diagonal ↖
end

"""
    triangular_lattice_graph(nx, ny) -> (SimpleGraph{Int}, Vector{Tuple{Float64,Float64}})

Build the full `nx × ny` triangular-lattice graph directly from integer `(i, j)`
indices, **without** floating-point distance computation.

The connectivity rule is `triangular_adjacency`: each vertex connects to its
horizontal, vertical, and one diagonal neighbour (six neighbours in the interior).
Physical positions follow the parallelogram layout used by
`get_physical_positions(Triangular(), ...)` and are returned as the second
element of the tuple for use in visualisation or UDG construction.

# Arguments
- `nx`: Number of columns (x direction)
- `ny`: Number of rows    (y direction)

# Returns
- `(SimpleGraph{Int}, Vector{Tuple{Float64,Float64}})`: graph and vertex positions

# Example
```julia
g, pos = triangular_lattice_graph(3, 3)   # 9-vertex triangular lattice
```
"""
function triangular_lattice_graph(nx::Int, ny::Int)
    coords = Tuple{Int,Int}[(i, j) for i in 1:nx for j in 1:ny]
    n      = length(coords)
    g      = SimpleGraph(n)
    idx    = Dict{Tuple{Int,Int}, Int}(coords[k] => k for k in 1:n)

    for k in 1:n
        i, j = coords[k]
        for (di, dj) in ((1, 0), (0, 1), (1, -1))   # 3 undirected half-edges
            nb = (i + di, j + dj)
            if haskey(idx, nb)
                add_edge!(g, k, idx[nb])
            end
        end
    end

    pos = get_physical_positions(Triangular(), coords)
    return g, pos
end
