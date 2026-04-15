# visualize.jl Refactor Design

## Goal

Simplify `src/utils/visualize.jl` by removing duplicate code, eliminating the x/y coordinate swap hack, and adding shape-aware lattice grid backgrounds. Reduce from 14 internal functions + 5 public functions to 6 internal + 3 public.

## Public API (3 functions)

### `plot_gadget(gadget, path; shape=nothing, ...)`

The single Gadget plotting function. Replaces both old `plot_single_gadget` and `plot_gadget`.

**Parameters:**
- `gadget::Gadget` — the gadget to plot
- `save_path::String` — output file path (.pdf/.svg/.png)
- `shape::Union{Nothing, String}=nothing` — lattice type ("grid", "KSG", "TLSG") for background grid
- `plot_size=400` — canvas size in pixels
- `margin=30` — margin around the plot
- `preserve_aspect_ratio=true`
- `background_grid=false` — draw lattice grid when shape is provided
- `show_weights=true` — show vertex weight numbers
- `show_edge_weights=true` — show edge weight labels (QUBO)
- `round_weights=false` — round weights to integers for display
- `discrete_color_scheme=ColorSchemes.seaborn_bright`
- `continuous_color_scheme=ColorSchemes.viridis`

**Behavior:**
1. Filter weight=0 vertices via `_filter_zero_weights(gadget)`
2. Compute layout via `_positions_to_layout(pos, plot_size, margin)`
3. If `background_grid && shape !== nothing`, draw lattice grid via `_draw_lattice_grid(...)`
4. Draw edges (with optional edge weight labels for QUBO)
5. Draw vertices with colors and optional weight labels
6. Draw PIN labels for pin vertices

**Coordinate handling:** Use `gadget.pos` directly as `(x, y)`. No x/y swap. The new `_compute_physical_positions` already outputs correct `(x, y)`.

### `plot_graph(graph, path; pos=nothing, ...)`

Simple graph plot. Unchanged API.

### `plot_single_gadget_mis(gadget, path; ...)`

MIS coloring visualization. Unchanged — uses spring layout, no coordinates.

## Internal Functions (6 functions)

### `_positions_to_layout(pos, plot_size, margin; preserve_aspect_ratio=true) -> Vector{Point}`

Replaces `_map_to_symmetric_range`, `_scale_values_with_margin`, `_layout_from_positions`, `_layout_from_points`, and `_map_and_scale`. One function that:
1. Extracts x and y from `pos::Vector{Tuple{Float64, Float64}}`
2. Centers around origin
3. Scales to fit `plot_size` with `margin`
4. Returns `Vector{Point}`

No x/y swap.

### `_filter_zero_weights(gadget) -> (new_graph, valid_vertices, old_to_new, new_weights, new_pins, new_pos, new_edge_weights, new_edge_list)`

Extracts the weight=0 filtering logic from `plot_gadget`. Returns all the filtered data needed for plotting. If no zero-weight vertices exist, returns the original data unchanged (fast path).

### `_generate_vertex_color(weights, discrete_scheme, continuous_scheme) -> Vector{Color}`

Unchanged.

### `_draw_lattice_grid(shape, positions, scale_x, scale_y, x_range, y_range)`

New. Draws lattice background lines based on shape:
- `"grid"`: horizontal + vertical lines only
- `"KSG"`: horizontal + vertical + both diagonal lines
- `"TLSG"`: horizontal + vertical + one diagonal direction (the connected one: `(+1,+1)/(-1,-1)`)

Uses the centered/scaled coordinate system. Lines drawn in light gray with thin stroke.

### `_with_drawing(fn, path, w, h) -> path`

Unchanged.

### `_suffix_path(path, suffix) -> String`

Unchanged.

## Deleted Functions

| Function | Reason |
|----------|--------|
| `plot_single_gadget` | Merged into `plot_gadget` |
| `plot_single_gadget_new` | Deprecated alias, delete |
| `_map_and_scale` | Contains x/y swap hack, replaced by `_positions_to_layout` |
| `_map_to_symmetric_range` | Absorbed into `_positions_to_layout` |
| `_scale_values_with_margin` | Absorbed into `_positions_to_layout` |
| `_layout_from_positions` | Absorbed into `_positions_to_layout` |
| `_layout_from_points` | Duplicate, absorbed into `_positions_to_layout` |
| `_generate_mask` | Only used by deleted `plot_single_gadget` |
| `_generate_pin_mask` | Only used by deleted `plot_single_gadget` |

## Test Changes

Update `test/utils/visualize.jl`:
- Remove `plot_single_gadget` test calls (function deleted)
- Remove `_layout_from_points` test (function deleted)
- Add test for `_positions_to_layout`
- Add tests for `plot_gadget` with `shape` and `background_grid=true`
- Keep `plot_single_gadget_mis` and `plot_graph` tests unchanged

## Export Changes

In `src/GadgetSearch.jl`:
- `plot_gadget` already exported — no change needed
- `plot_single_gadget` was not exported — no change needed

## Constraints

- `Gadget` struct is NOT changed (no shape field added)
- Shape info is passed as a keyword argument to `plot_gadget`
- All visual output stays Luxor-based (.pdf/.svg/.png by extension)
