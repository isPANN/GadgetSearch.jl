function _map_to_symmetric_range(x_min, x_max, y_min, y_max)
    x_center = (x_min + x_max) / 2
    y_center = (y_min + y_max) / 2

    map_xcoordinates = x -> x - x_center
    map_ycoordinates = y -> y - y_center

    return map_xcoordinates, map_ycoordinates
end

function _scale_values_with_margin(xvals, yvals, xmax_abs_val, ymax_abs_val, margin, preserve_aspect_ratio)
    allowed_x_max = xmax_abs_val / 2 - margin
    allowed_y_max = ymax_abs_val / 2 - margin
    
    x_denom = max(maximum(abs.(xvals)), eps(Float64))
    y_denom = max(maximum(abs.(yvals)), eps(Float64))
    scale_factor_x = allowed_x_max / x_denom
    scale_factor_y = allowed_y_max / y_denom

    if preserve_aspect_ratio
        scale_factor = min(scale_factor_x, scale_factor_y)
        return scale_factor, scale_factor
    end
    return scale_factor_x, scale_factor_y
end

function _map_and_scale(gadget::Gadget, x_plot_size, y_plot_size, margin, preserve_aspect_ratio)
    #Note: In Luxor, the X-axis represents the horizontal direction, and the Y-axis represents the vertical direction, so we need to adjust the coordinate order accordingly.
    y_vals = [pos[1] for pos in gadget.pos]
    x_vals = [pos[2] for pos in gadget.pos]
    
    x_min, x_max = minimum(x_vals), maximum(x_vals)
    y_min, y_max = minimum(y_vals), maximum(y_vals)

    map_xcoordinates, map_ycoordinates = _map_to_symmetric_range(x_min, x_max, y_min, y_max)

    x_new_vals = map_xcoordinates.(x_vals)
    y_new_vals = map_ycoordinates.(y_vals)
    
    x_scale_factor, y_scale_factor = _scale_values_with_margin(x_new_vals, y_new_vals, x_plot_size, y_plot_size, margin, preserve_aspect_ratio)

    return x_new_vals, y_new_vals, x_scale_factor, y_scale_factor
end

function _generate_vertex_color(weights::AbstractVector, discrete_color_scheme, continuous_color_scheme)
    if any(w < 1 || w > length(discrete_color_scheme) || typeof(w) != Int for w in weights)
        # @info "Using continuous color scheme."
        min_w, max_w = minimum(weights), maximum(weights)
        # generate colormap
        cmap = get(continuous_color_scheme, range(0, 1, length=length(weights)))
        # handle equal weights to avoid divide-by-zero
        if max_w == min_w
            mid = clamp(round(Int, length(cmap)/2), 1, length(cmap))
            node_colors = fill(cmap[mid], length(weights))
        else
            # assign colors based on normalized weights
            node_colors = [cmap[clamp(floor(Int, (w - min_w) / (max_w - min_w) * (length(cmap) - 1) + 1), 1, length(cmap))] for w in weights]
        end
    else
        # @info "Using discrete color scheme."
        node_colors = [discrete_color_scheme[w] for w in weights]
    end
    return node_colors
    
end

function _generate_mask(index_map, n)
    mask = fill("", n)  # generate an empty mask
    for (i, idx) in enumerate(index_map)
        mask[idx] = string(i)  # use the index_map to fill the mask
    end
    return mask
end


function _generate_pin_mask(special_index, n; normal_fill::String="", special_fill::String="")
    mask = fill(normal_fill, n)  # generate an empty mask
    for (i, idx) in enumerate(special_index)
        if length(special_fill) == 0
            mask[idx] = string(i)  # use the special_index to fill the mask
        else
            mask[idx] = special_fill  # use the special fill value
        end
    end
    return mask
end

# ----- Drawing helpers -----
function _with_drawing(drawfn::Function, save_path::AbstractString, w::Integer, h::Integer)
    # Let Luxor auto-select backend by file extension (.pdf/.svg/.png)
    Drawing(w, h, save_path)
    try
        origin()
        drawfn()
    finally
        finish()
    end
    return save_path
end

function _suffix_path(path::AbstractString, suffix::AbstractString)
    base, ext = splitext(path)
    return string(base, suffix, ext)
end

function _layout_from_positions(pos::Vector{Tuple{Float64, Float64}}, plot_size::Int, margin::Int; preserve_aspect_ratio::Bool=true)
    x_vals = [p[1] for p in pos]
    y_vals = [p[2] for p in pos]
    x_min, x_max = minimum(x_vals), maximum(x_vals)
    y_min, y_max = minimum(y_vals), maximum(y_vals)
    map_xcoordinates, map_ycoordinates = _map_to_symmetric_range(x_min, x_max, y_min, y_max)
    x_new_vals = map_xcoordinates.(x_vals)
    y_new_vals = map_ycoordinates.(y_vals)
    x_scale_factor, y_scale_factor = _scale_values_with_margin(x_new_vals, y_new_vals, plot_size, plot_size, margin, preserve_aspect_ratio)
    return [Point(x*x_scale_factor, y*y_scale_factor) for (x,y) in zip(x_new_vals, y_new_vals)]
end

function plot_single_gadget(gadget::Gadget, save_path::String; plot_size=400, margin=20, preserve_aspect_ratio=true, discrete_color_scheme=ColorSchemes.seaborn_bright, continuous_color_scheme=ColorSchemes.viridis)
    x_new_vals, y_new_vals, x_scale_factor, y_scale_factor = _map_and_scale(gadget, plot_size, plot_size, margin, preserve_aspect_ratio)
    pts = [Point(x*x_scale_factor, y*y_scale_factor) for (x,y) in zip(x_new_vals, y_new_vals)]
    node_colors = _generate_vertex_color(gadget.weights, discrete_color_scheme, continuous_color_scheme)
    mask = _generate_mask(gadget.pins, nv(gadget.graph))

    _with_drawing(save_path, plot_size, plot_size) do
        background("white")
        drawgraph(gadget.graph,
                  vertexlabels=mask,
                  vertexshapesizes=10,
                  vertexlabelfontsizes=15,
                  layout=pts,
                  edgestrokeweights=1,
                  vertexfillcolors=node_colors)
    end
    println("Drawing saved as $save_path")
end


function plot_gadget(
                            gadget::Gadget, save_path::String; 
                            plot_size=400, margin=30, 
                            preserve_aspect_ratio=true, 
                            background_grid=false,
                            show_weights=true,
                            show_edge_weights=true,
                            round_weights=false,
                            discrete_color_scheme=ColorSchemes.seaborn_bright, 
                            continuous_color_scheme=ColorSchemes.viridis
                            )

    # Create a new graph with only non-zero weight vertices
    valid_vertices = findall(w -> w != 0, gadget.weights)
    new_graph = SimpleGraph(length(valid_vertices))
    
    # Create mapping from old vertex indices to new ones
    old_to_new = Dict(v => i for (i, v) in enumerate(valid_vertices))
    
    # Track edge weights for QUBO gadgets
    has_edge_weights = !isempty(gadget.edge_weights)
    edge_weight_map = Dict{Tuple{Int,Int}, Float64}()
    
    # Add edges between valid vertices and map edge weights
    if has_edge_weights
        for (idx, (src, dst)) in enumerate(gadget.edge_list)
            if src in valid_vertices && dst in valid_vertices
                new_src = old_to_new[src]
                new_dst = old_to_new[dst]
                add_edge!(new_graph, new_src, new_dst)
                # Store edge weight (normalized order)
                edge_key = new_src < new_dst ? (new_src, new_dst) : (new_dst, new_src)
                edge_weight_map[edge_key] = gadget.edge_weights[idx]
            end
        end
    else
        for e in edges(gadget.graph)
            src, dst = Graphs.src(e), Graphs.dst(e)
            if src in valid_vertices && dst in valid_vertices
                add_edge!(new_graph, old_to_new[src], old_to_new[dst])
            end
        end
    end
    
    # Create new weights and pins arrays
    new_weights = gadget.weights[valid_vertices]
    new_pins = [old_to_new[p] for p in gadget.pins if p in valid_vertices]
    
    # Create new positions array
    new_pos = gadget.pos[valid_vertices]
    
    # Create new gadget with filtered data (for compatibility)
    new_gadget = Gadget(gadget.ground_states, new_graph, new_pins, new_weights, new_pos)

    x_new_vals, y_new_vals, x_scale_factor, y_scale_factor = _map_and_scale(new_gadget, plot_size, plot_size, margin, preserve_aspect_ratio)
    x_max, y_max = maximum(x_new_vals), maximum(y_new_vals)
    x_min, y_min = minimum(x_new_vals), minimum(y_new_vals)

    _with_drawing(save_path, plot_size, plot_size) do
        background("white")

        if background_grid
            x_list = collect(x_min:x_max) .* x_scale_factor
            y_list = collect(y_min:y_max) .* y_scale_factor

            sethue("gray")
            setline(0.3)
            for x in x_list
                line(Point(x, y_min*y_scale_factor), Point(x, y_max*y_scale_factor), :stroke)
            end
            for y in y_list
                line(Point(x_min*x_scale_factor, y), Point(x_max*x_scale_factor, y), :stroke)
            end
        end

        sethue("black")
        pts = [Point(x*x_scale_factor, y*y_scale_factor) for (x,y) in zip(x_new_vals, y_new_vals)]
        node_colors = _generate_vertex_color(new_gadget.weights, discrete_color_scheme, continuous_color_scheme)

        # Draw edges first (so edge labels appear below vertices)
        if has_edge_weights && show_edge_weights
            for e in edges(new_graph)
                src, dst = Graphs.src(e), Graphs.dst(e)
                edge_key = src < dst ? (src, dst) : (dst, src)
                
                # Draw edge
                sethue("black")
                setline(2)
                line(pts[src], pts[dst], :stroke)
                
                # Draw edge weight label offset from edge
                if haskey(edge_weight_map, edge_key)
                    weight = edge_weight_map[edge_key]
                    
                    # Calculate midpoint
                    mid_x = (pts[src].x + pts[dst].x) / 2
                    mid_y = (pts[src].y + pts[dst].y) / 2
                    
                    # Calculate perpendicular offset direction
                    dx = pts[dst].x - pts[src].x
                    dy = pts[dst].y - pts[src].y
                    length = sqrt(dx^2 + dy^2)
                    
                    # Perpendicular vector (rotated 90 degrees)
                    perp_x = -dy / length
                    perp_y = dx / length
                    
                    # Offset distance (adjust based on edge direction to avoid overlap)
                    offset = 12
                    label_point = Point(mid_x + perp_x * offset, mid_y + perp_y * offset)
                    
                    # Background for text
                    sethue("white")
                    fontsize(11)
                    weight_text = round_weights ? string(round(Int, weight)) : string(weight)
                    text_width = textextents(weight_text)[3]
                    text_height = textextents(weight_text)[4]
                    box(label_point, text_width + 4, text_height + 4, :fill)
                    
                    # Draw text
                    sethue("blue")
                    fontsize(11)
                    text(weight_text, label_point, halign=:center, valign=:middle)
                end
            end
        end

        drawgraph(new_gadget.graph,
                  vertexfunction=(v, c) -> begin
                      index = findfirst(==(v), new_gadget.pins)
                      sethue(node_colors[v])
                      translate(pts[v])
                      circle(O, 10, :fill)
                      sethue("white")
                      fontsize(15)
                      if show_weights
                          if round_weights
                              text("$(round(Int, new_gadget.weights[v]))", O, halign=:center, valign=:middle)
                          else
                              text("$(new_gadget.weights[v])", O, halign=:center, valign=:middle)
                          end
                      end
                      sethue("red")
                      if !isnothing(index)
                          # Place PIN label with a consistent vertical offset from the vertex.
                          # Using relative coordinates avoids accumulating transforms.
                          # If y is exactly 0, ensure a nonzero offset by preferring upward placement.
                          y = pts[v].y
                          dy = y >= 0 ? 20 : -20
                          text("PIN $index", Point(0, dy), halign=:center, valign=:middle)
                      end
                  end,
                  edgestrokeweights=has_edge_weights && show_edge_weights ? 0 : 2,  # Hide edges if we drew them manually
                  layout=pts,
                 )
    end
    println("Drawing saved as $save_path")
end

function plot_single_gadget_new(args...; kwargs...)
    @warn "plot_single_gadget_new is deprecated; use plot_gadget instead."
    return plot_gadget(args...; kwargs...)
end
function plot_single_gadget_mis(
                            gadget::Gadget, save_path::String; 
                            plot_size=400, margin=30
                            )
    masks, n = find_maximal_independent_sets(gadget.graph)
    for i in 1:n
        save_single_path = _suffix_path(save_path, "_" * string(i))
        _with_drawing(save_single_path, plot_size, plot_size) do
            # Build color vector from bit mask for each vertex
            m = masks[i]
            colors = [(((m >> (v - 1)) & 0x1) == 1) ? RGB(1, 0, 0) : RGB(0, 0, 0) for v in 1:nv(gadget.graph)]
            drawgraph(gadget.graph,
                      vertexlabels=collect(1:nv(gadget.graph)),
                      vertexshapesizes=18,
                      vertexlabelfontsizes=22,
                      vertexfillcolors=colors,
                      edgestrokeweights=3,
                      layout=spring)
        end
        println("Drawing saved as $save_single_path")
    end
end

function plot_graph(g::SimpleGraph, save_path::String; 
                   pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing,
                   plot_size=400, 
                   margin=30,
                   vertex_size=10,
                   vertex_label_size=15,
                   edge_width=2)
    _with_drawing(save_path, plot_size, plot_size) do
        background("white")

        layout = pos === nothing ? stress : _layout_from_positions(pos, plot_size, margin; preserve_aspect_ratio=true)
        sethue("black")
        drawgraph(g,
                  vertexlabels=collect(1:nv(g)),
                  vertexshapesizes=vertex_size,
                  vertexlabelfontsizes=vertex_label_size,
                  edgestrokeweights=edge_width,
                  layout=layout)
    end
    println("Drawing saved as $save_path")
end

function _layout_from_points(points, plot_size::Int, margin::Int; preserve_aspect_ratio::Bool=true)
    x_vals = [p[1] for p in points]
    y_vals = [p[2] for p in points]
    x_min, x_max = minimum(x_vals), maximum(x_vals)
    y_min, y_max = minimum(y_vals), maximum(y_vals)
    map_xcoordinates, map_ycoordinates = _map_to_symmetric_range(x_min, x_max, y_min, y_max)
    x_new_vals = map_xcoordinates.(x_vals)
    y_new_vals = map_ycoordinates.(y_vals)
    x_scale_factor, y_scale_factor = _scale_values_with_margin(x_new_vals, y_new_vals, plot_size, plot_size, margin, preserve_aspect_ratio)
    return [Point(x*x_scale_factor, y*y_scale_factor) for (x, y) in zip(x_new_vals, y_new_vals)]
end

function _format_tensor_entry(value::Real)
    if isinf(value) && value < 0
        return "-Inf"
    end

    rounded = round(float(value); digits=2)
    rounded_int = round(Int, rounded)
    if isapprox(rounded, rounded_int; atol=1e-9)
        return string(rounded_int)
    end

    return string(rounded)
end

function _format_reduced_tensor_lines(reduced::AbstractVector{<:Real}; entries_per_line::Int=6)
    entries_per_line >= 1 || throw(ArgumentError("entries_per_line must be positive"))
    formatted = _format_tensor_entry.(reduced)
    lines = String[]

    for start_idx in 1:entries_per_line:length(formatted)
        end_idx = min(start_idx + entries_per_line - 1, length(formatted))
        prefix = start_idx == 1 ? "alpha = " : "        "
        push!(lines, prefix * "[" * join(formatted[start_idx:end_idx], ", ") * "]")
    end

    return lines
end

function _equivalent_representation_lines(
    base_graph::SimpleGraph{Int},
    repr_graph::SimpleGraph{Int},
    repr_boundary::Vector{Int};
    show_reduced_tensor::Bool=true,
    tensor_entries_per_line::Int=6,
)
    reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(repr_graph, repr_boundary))))
    lines = String[
        "|V|=$(nv(repr_graph)), |E|=$(ne(repr_graph)), extra=$(nv(repr_graph) - nv(base_graph))",
        "boundary = $(repr_boundary)",
        "inf mask = 0x$(string(inf_mask(reduced); base=16))",
    ]

    if show_reduced_tensor
        append!(lines, _format_reduced_tensor_lines(reduced; entries_per_line=tensor_entries_per_line))
    end

    return lines
end

function _draw_equivalent_representation_graph(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    layout_points::AbstractVector;
    vertex_size::Int=12,
    vertex_label_size::Int=11,
    edge_width::Real=1.5,
)
    boundary_order = Dict(vertex => idx for (idx, vertex) in enumerate(boundary))
    vertex_labels = [
        haskey(boundary_order, vertex) ? "$(vertex)/b$(boundary_order[vertex])" : string(vertex)
        for vertex in 1:nv(graph)
    ]
    vertex_colors = [
        haskey(boundary_order, vertex) ? RGB(0.95, 0.69, 0.20) : RGB(0.66, 0.82, 0.96)
        for vertex in 1:nv(graph)
    ]

    sethue("black")
    drawgraph(
        graph,
        vertexlabels=vertex_labels,
        vertexshapesizes=vertex_size,
        vertexlabelfontsizes=vertex_label_size,
        vertexfillcolors=vertex_colors,
        edgestrokeweights=edge_width,
        layout=layout_points,
    )
end

"""
    plot_equivalent_representations(graph, boundary, save_path; max_added_vertices=0, preserve_boundary_roles=true, columns=3, panel_graph_size=220)

Create a gallery-style visualization of all distinct outputs of
[`equivalent_representations`](@ref). Each panel shows one representation with
boundary-order annotations and a summary of the reduced alpha tensor used for
deduplication.

# Arguments
- `graph::SimpleGraph{Int}`: Input graph
- `boundary::Vector{Int}`: Boundary vertices of the input graph
- `save_path::String`: Output image path (`.svg`, `.png`, `.pdf`, ...)

# Keywords
- `max_added_vertices::Int=0`: Forwarded to [`equivalent_representations`](@ref)
- `preserve_boundary_roles::Bool=true`: Forwarded to
  [`equivalent_representations`](@ref)
- `columns::Int=3`: Number of panels per row
- `panel_graph_size::Int=220`: Size of each graph drawing region
- `panel_margin::Int=18`: Inner graph margin within each panel
- `panel_padding::Int=18`: Outer padding inside each panel
- `panel_gap::Int=18`: Gap between panels
- `title_height::Int=56`: Reserved height for the figure title
- `show_reduced_tensor::Bool=true`: Show flattened reduced alpha tensor lines
- `tensor_entries_per_line::Int=6`: How many tensor entries to show per line
- `pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing`: Optional
  positions for the original graph. These are reused when the representation has
  the same number of vertices as the input graph.
"""
function plot_equivalent_representations(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    save_path::String;
    max_added_vertices::Int=0,
    preserve_boundary_roles::Bool=true,
    columns::Int=3,
    panel_graph_size::Int=220,
    panel_margin::Int=18,
    panel_padding::Int=18,
    panel_gap::Int=18,
    title_height::Int=56,
    show_reduced_tensor::Bool=true,
    tensor_entries_per_line::Int=6,
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing,
)
    columns >= 1 || throw(ArgumentError("columns must be positive"))
    panel_graph_size > 0 || throw(ArgumentError("panel_graph_size must be positive"))
    max_added_vertices >= 0 || throw(ArgumentError("max_added_vertices must be non-negative"))

    reprs = equivalent_representations(
        graph,
        boundary;
        max_added_vertices=max_added_vertices,
        preserve_boundary_roles=preserve_boundary_roles,
    )
    annotation_lines = [
        _equivalent_representation_lines(graph, repr_graph, repr_boundary;
            show_reduced_tensor=show_reduced_tensor,
            tensor_entries_per_line=tensor_entries_per_line,
        )
        for (repr_graph, repr_boundary) in reprs
    ]

    line_height = 16
    annotation_height = max(90, maximum(length.(annotation_lines)) * line_height + 20)
    panel_width = panel_graph_size + 2 * panel_padding
    panel_height = panel_graph_size + annotation_height + 2 * panel_padding
    rows = cld(length(reprs), columns)
    canvas_width = columns * panel_width + (columns + 1) * panel_gap
    canvas_height = title_height + rows * panel_height + (rows + 1) * panel_gap

    _with_drawing(save_path, canvas_width, canvas_height) do
        background("white")

        sethue("black")
        fontsize(20)
        text(
            "Equivalent Representations",
            Point(-canvas_width / 2 + panel_gap, -canvas_height / 2 + 24);
            halign=:left,
            valign=:middle,
        )

        fontsize(11)
        text(
            "Input boundary = $(boundary), max_added_vertices = $(max_added_vertices), preserve_boundary_roles = $(preserve_boundary_roles), orange nodes are boundary vertices labeled as vertex/boundary-order.",
            Point(-canvas_width / 2 + panel_gap, -canvas_height / 2 + 44);
            halign=:left,
            valign=:middle,
        )

        for (idx, (repr_graph, repr_boundary)) in enumerate(reprs)
            row = fld(idx - 1, columns)
            col = mod(idx - 1, columns)
            left = -canvas_width / 2 + panel_gap + col * (panel_width + panel_gap)
            top = -canvas_height / 2 + title_height + panel_gap + row * (panel_height + panel_gap)

            sethue(RGB(0.97, 0.97, 0.98))
            box(Point(left + panel_width / 2, top + panel_height / 2), panel_width, panel_height, :fill)
            sethue(RGB(0.82, 0.84, 0.88))
            box(Point(left + panel_width / 2, top + panel_height / 2), panel_width, panel_height, :stroke)

            sethue("black")
            fontsize(14)
            text(
                "repr $(idx)",
                Point(left + panel_padding, top + 18);
                halign=:left,
                valign=:middle,
            )

            layout_points =
                if !isnothing(pos) && length(pos) == nv(repr_graph)
                    _layout_from_positions(pos, panel_graph_size, panel_margin; preserve_aspect_ratio=true)
                else
                    _layout_from_points(stress(repr_graph), panel_graph_size, panel_margin; preserve_aspect_ratio=true)
                end

            gsave()
            translate(Point(left + panel_width / 2, top + panel_padding + panel_graph_size / 2))
            _draw_equivalent_representation_graph(repr_graph, repr_boundary, layout_points)
            grestore()

            fontsize(10)
            for (line_idx, line) in enumerate(annotation_lines[idx])
                text(
                    line,
                    Point(left + panel_padding, top + panel_padding + panel_graph_size + 12 + (line_idx - 1) * line_height);
                    halign=:left,
                    valign=:middle,
                )
            end
        end
    end

    println("Drawing saved as $save_path")
    return save_path
end

function _simple_path_between(graph::SimpleGraph{Int}, start::Int, goal::Int)
    start == goal && return [start]

    parents = fill(0, nv(graph))
    visited = falses(nv(graph))
    queue = [start]
    visited[start] = true
    head = 1

    while head <= length(queue)
        current = queue[head]
        head += 1

        for neighbor in neighbors(graph, current)
            visited[neighbor] && continue
            visited[neighbor] = true
            parents[neighbor] = current
            neighbor == goal && break
            push!(queue, neighbor)
        end

        visited[goal] && break
    end

    visited[goal] || throw(ArgumentError("no path exists between $start and $goal"))

    path = Int[goal]
    while path[end] != start
        push!(path, parents[path[end]])
    end
    reverse!(path)
    return path
end

function _crossing_arm_paths(graph::SimpleGraph{Int}, boundary::Vector{Int})
    length(boundary) == 4 || throw(ArgumentError("crossing gallery expects exactly 4 boundary vertices"))

    horizontal_path = _simple_path_between(graph, boundary[1], boundary[3])
    vertical_path = _simple_path_between(graph, boundary[2], boundary[4])
    covered_vertices = Set(vcat(horizontal_path, vertical_path))
    covered_vertices == Set(1:nv(graph)) || throw(ArgumentError(
        "crossing gallery expects a graph formed only by the 1-3 and 2-4 arms"))

    return horizontal_path, vertical_path
end

function _crossing_signature(graph::SimpleGraph{Int}, boundary::Vector{Int})
    horizontal_path, vertical_path = _crossing_arm_paths(graph, boundary)
    return (
        horizontal_subdivisions=length(horizontal_path) - 2,
        vertical_subdivisions=length(vertical_path) - 2,
        horizontal_path=horizontal_path,
        vertical_path=vertical_path,
    )
end

function _can_draw_with_crossing_geometry(graph::SimpleGraph{Int}, boundary::Vector{Int})
    try
        _crossing_signature(graph, boundary)
        return true
    catch err
        if err isa ArgumentError
            return false
        end
        rethrow()
    end
end

function _crossing_layout_map(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    horizontal_path::Vector{Int},
    vertical_path::Vector{Int},
)
    axis_extent = 0.86
    crossing_gap = 0.16
    positions = Dict{Int, Point}()

    function _axis_coordinate(t::Real)
        if t < 0.5
            return -axis_extent + (t / 0.5) * (axis_extent - crossing_gap)
        elseif t > 0.5
            return crossing_gap + ((t - 0.5) / 0.5) * (axis_extent - crossing_gap)
        end

        # Keep the geometric crossing empty: a path midpoint is nudged off-center
        # instead of becoming an artificial shared graph vertex.
        return crossing_gap
    end

    for (idx, vertex) in enumerate(horizontal_path)
        t = length(horizontal_path) == 1 ? 0.5 : (idx - 1) / (length(horizontal_path) - 1)
        positions[vertex] = Point(_axis_coordinate(t), 0.0)
    end

    for (idx, vertex) in enumerate(vertical_path)
        t = length(vertical_path) == 1 ? 0.5 : (idx - 1) / (length(vertical_path) - 1)
        positions[vertex] = Point(0.0, _axis_coordinate(t))
    end

    length(positions) == nv(graph) || throw(ArgumentError(
        "crossing gallery could not assign positions to every vertex"))

    return positions
end

function _draw_crossing_path(
    path_vertices::Vector{Int},
    positions::Dict{Int, Point},
    scale::Real;
    gap_at_center::Bool=false,
    center_gap_radius::Real=0.0,
)
    for (u, v) in zip(path_vertices[1:end-1], path_vertices[2:end])
        start_point = positions[u] * scale
        end_point = positions[v] * scale

        if gap_at_center
            if isapprox(start_point.y, 0.0; atol=1e-9) &&
               isapprox(end_point.y, 0.0; atol=1e-9) &&
               start_point.x * end_point.x < 0
                line(start_point, Point(-center_gap_radius, 0), :stroke)
                line(Point(center_gap_radius, 0), end_point, :stroke)
                continue
            elseif isapprox(start_point.x, 0.0; atol=1e-9) &&
                   isapprox(end_point.x, 0.0; atol=1e-9) &&
                   start_point.y * end_point.y < 0
                line(start_point, Point(0, -center_gap_radius), :stroke)
                line(Point(0, center_gap_radius), end_point, :stroke)
                continue
            end
        end

        line(start_point, end_point, :stroke)
    end
end

function _draw_crossing_variant_graph(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    panel_graph_size::Int,
)
    signature = _crossing_signature(graph, boundary)
    positions = _crossing_layout_map(
        graph,
        boundary,
        signature.horizontal_path,
        signature.vertical_path,
    )
    scale = panel_graph_size / 2
    boundary_set = Set(boundary)
    horizontal_center_vertex = findfirst(v -> positions[v] == Point(0.0, 0.0), signature.horizontal_path)
    vertical_center_vertex = findfirst(v -> positions[v] == Point(0.0, 0.0), signature.vertical_path)
    horizontal_owns_center = !isnothing(horizontal_center_vertex)
    vertical_owns_center = !isnothing(vertical_center_vertex)
    center_gap_radius = 0.12 * scale

    setline(3.2)
    sethue(RGB(0.22, 0.45, 0.82))
    _draw_crossing_path(
        signature.horizontal_path,
        positions,
        scale;
        gap_at_center=vertical_owns_center && !horizontal_owns_center,
        center_gap_radius=center_gap_radius,
    )

    sethue(RGB(0.84, 0.31, 0.26))
    _draw_crossing_path(
        signature.vertical_path,
        positions,
        scale;
        gap_at_center=horizontal_owns_center && !vertical_owns_center,
        center_gap_radius=center_gap_radius,
    )

    for vertex in 1:nv(graph)
        point = positions[vertex] * scale
        if vertex in boundary_set
            sethue(RGB(0.95, 0.69, 0.20))
            circle(point, 13, :fill)
        else
            sethue(RGB(0.78, 0.88, 0.98))
            circle(point, 11, :fill)
        end
        sethue("black")
        circle(point, vertex in boundary_set ? 13 : 11, :stroke)
        fontsize(12)
        text(string(vertex), point; halign=:center, valign=:middle)
    end

    return signature
end

function _crossing_panel_title(horizontal_subdivisions::Int, vertical_subdivisions::Int)
    if horizontal_subdivisions == 0 && vertical_subdivisions == 0
        return "canonical crossing"
    elseif horizontal_subdivisions > 0 && vertical_subdivisions == 0
        return "stretch 1-3 arm"
    elseif horizontal_subdivisions == 0 && vertical_subdivisions > 0
        return "stretch 2-4 arm"
    end

    return "stretch both arms"
end

"""
    plot_crossing_equivalence_gallery(graph, boundary, save_path; ...)

Create a report-friendly gallery for a canonical 4-pin crossing target with
fixed boundary roles. The figure uses a fixed geometric layout:
`1 = left`, `2 = top`, `3 = right`, `4 = bottom`.

The gallery highlights structural examples from
[`equivalent_representations`](@ref) using the fixed crossing geometry. Logical
flips remain tensor-level variants and are therefore described in text rather
than drawn as separate graph rewrites here.
"""
function plot_crossing_equivalence_gallery(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    save_path::String;
    max_added_vertices::Int=4,
    preserve_boundary_roles::Bool=true,
    example_signatures::Vector{Tuple{Int, Int}}=[(0, 0), (1, 0), (0, 1), (2, 2)],
    columns::Int=2,
    panel_graph_size::Int=220,
    panel_padding::Int=18,
    panel_gap::Int=18,
    title_height::Int=78,
)
    preserve_boundary_roles || throw(ArgumentError(
        "crossing gallery requires preserve_boundary_roles=true"))
    length(boundary) == 4 || throw(ArgumentError("crossing gallery expects exactly 4 boundary vertices"))
    max_added_vertices >= 0 || throw(ArgumentError("max_added_vertices must be non-negative"))
    columns >= 1 || throw(ArgumentError("columns must be positive"))

    reprs = equivalent_representations(
        graph,
        boundary;
        max_added_vertices=max_added_vertices,
        preserve_boundary_roles=preserve_boundary_roles,
    )

    repr_lookup = Dict{Tuple{Int, Int}, NamedTuple{(:graph, :boundary, :signature), Tuple{SimpleGraph{Int}, Vector{Int}, NamedTuple}}}()
    for (repr_graph, repr_boundary) in reprs
        signature = _crossing_signature(repr_graph, repr_boundary)
        key = (signature.horizontal_subdivisions, signature.vertical_subdivisions)
        haskey(repr_lookup, key) && continue
        repr_lookup[key] = (graph=repr_graph, boundary=repr_boundary, signature=signature)
    end

    selected = [
        merge((key=signature_key,), repr_lookup[signature_key])
        for signature_key in example_signatures
        if haskey(repr_lookup, signature_key)
    ]
    isempty(selected) && throw(ArgumentError("no requested crossing examples were found"))

    line_height = 16
    annotation_height = 78
    panel_width = panel_graph_size + 2 * panel_padding
    panel_height = panel_graph_size + annotation_height + 2 * panel_padding
    rows = cld(length(selected), columns)
    canvas_width = columns * panel_width + (columns + 1) * panel_gap
    canvas_height = title_height + rows * panel_height + (rows + 1) * panel_gap

    _with_drawing(save_path, canvas_width, canvas_height) do
        background("white")

        sethue("black")
        fontsize(20)
        text(
            "Canonical Crossing Equivalence Gallery",
            Point(-canvas_width / 2 + panel_gap, -canvas_height / 2 + 24);
            halign=:left,
            valign=:middle,
        )

        fontsize(11)
        text(
            "Fixed roles: 1 = left, 2 = top, 3 = right, 4 = bottom. Blue shows the 1-3 arm; red shows the 2-4 arm.",
            Point(-canvas_width / 2 + panel_gap, -canvas_height / 2 + 46);
            halign=:left,
            valign=:middle,
        )

        text(
            "Panels show structural examples from equivalent_representations(...); logical flips remain tensor-level variants.",
            Point(-canvas_width / 2 + panel_gap, -canvas_height / 2 + 64);
            halign=:left,
            valign=:middle,
        )

        for (idx, example) in enumerate(selected)
            row = fld(idx - 1, columns)
            col = mod(idx - 1, columns)
            left = -canvas_width / 2 + panel_gap + col * (panel_width + panel_gap)
            top = -canvas_height / 2 + title_height + panel_gap + row * (panel_height + panel_gap)

            sethue(RGB(0.97, 0.97, 0.98))
            box(Point(left + panel_width / 2, top + panel_height / 2), panel_width, panel_height, :fill)
            sethue(RGB(0.82, 0.84, 0.88))
            box(Point(left + panel_width / 2, top + panel_height / 2), panel_width, panel_height, :stroke)

            sethue("black")
            fontsize(14)
            text(
                _crossing_panel_title(example.signature.horizontal_subdivisions, example.signature.vertical_subdivisions),
                Point(left + panel_padding, top + 18);
                halign=:left,
                valign=:middle,
            )

            gsave()
            translate(Point(left + panel_width / 2, top + panel_padding + panel_graph_size / 2))
            _draw_crossing_variant_graph(example.graph, example.boundary, panel_graph_size - 22)
            grestore()

            fontsize(10)
            annotation_lines = [
                "1-3 arm: $(join(example.signature.horizontal_path, "-"))",
                "2-4 arm: $(join(example.signature.vertical_path, "-"))",
                "extra vertices: $(nv(example.graph) - nv(graph))",
            ]
            for (line_idx, line) in enumerate(annotation_lines)
                text(
                    line,
                    Point(left + panel_padding, top + panel_padding + panel_graph_size + 12 + (line_idx - 1) * line_height);
                    halign=:left,
                    valign=:middle,
                )
            end
        end
    end

    println("Drawing saved as $save_path")
    return save_path
end

function _draw_report_box(
    center::Point,
    width::Real,
    height::Real,
    title::AbstractString;
    fill::Colorant=RGB(0.98, 0.99, 1.0),
    border::Colorant=RGB(0.76, 0.82, 0.92),
    header_fill::Colorant=RGB(0.90, 0.94, 0.99),
)
    sethue(fill)
    box(center, width, height, :fill)
    sethue(border)
    setline(1.6)
    box(center, width, height, :stroke)

    header_height = 32
    header_center = Point(center.x, center.y - height / 2 + header_height / 2)
    sethue(header_fill)
    box(header_center, width, header_height, :fill)
    sethue(border)
    box(header_center, width, header_height, :stroke)

    sethue("black")
    fontsize(16)
    text(
        title,
        Point(center.x - width / 2 + 16, center.y - height / 2 + 20);
        halign=:left,
        valign=:middle,
    )
end

function _draw_report_lines(
    top_left::Point,
    lines::AbstractVector{<:AbstractString};
    font_size::Real=11,
    line_height::Real=15,
)
    sethue("black")
    fontsize(font_size)
    for (idx, line) in enumerate(lines)
        text(
            line,
            Point(top_left.x, top_left.y + (idx - 1) * line_height);
            halign=:left,
            valign=:middle,
        )
    end
end

function _draw_report_arrow(start_point::Point, end_point::Point; label::Union{Nothing, AbstractString}=nothing)
    sethue(RGB(0.42, 0.47, 0.55))
    setline(2.4)
    line(start_point, end_point, :stroke)

    dx = end_point.x - start_point.x
    dy = end_point.y - start_point.y
    length = sqrt(dx^2 + dy^2)
    length <= eps(Float64) && return

    ux = dx / length
    uy = dy / length
    arrow_size = 10.0

    p1 = Point(
        end_point.x - ux * arrow_size - uy * arrow_size / 2,
        end_point.y - uy * arrow_size + ux * arrow_size / 2,
    )
    p2 = Point(
        end_point.x - ux * arrow_size + uy * arrow_size / 2,
        end_point.y - uy * arrow_size - ux * arrow_size / 2,
    )
    line(end_point, p1, :stroke)
    line(end_point, p2, :stroke)

    if !isnothing(label)
        sethue(RGB(0.33, 0.37, 0.44))
        fontsize(10)
        text(
            label,
            Point((start_point.x + end_point.x) / 2 + 10, (start_point.y + end_point.y) / 2);
            halign=:left,
            valign=:middle,
        )
    end
end

function _draw_graph_thumbnail(
    graph::SimpleGraph{Int},
    boundary::Vector{Int},
    center::Point,
    plot_size::Int;
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing,
    show_boundary_order::Bool=false,
    margin::Int=14,
    crossing_geometry::Bool=false,
)
    if crossing_geometry
        gsave()
        translate(center)
        _draw_crossing_variant_graph(graph, boundary, plot_size - 2 * margin)
        grestore()
        return
    end

    boundary_order = Dict(vertex => idx for (idx, vertex) in enumerate(boundary))
    vertex_labels = [
        if show_boundary_order && haskey(boundary_order, vertex)
            "$(vertex)/b$(boundary_order[vertex])"
        else
            string(vertex)
        end
        for vertex in 1:nv(graph)
    ]
    vertex_colors = [
        haskey(boundary_order, vertex) ? RGB(0.95, 0.69, 0.20) : RGB(0.66, 0.82, 0.96)
        for vertex in 1:nv(graph)
    ]

    layout_points =
        if !isnothing(pos) && length(pos) == nv(graph)
            _layout_from_positions(pos, plot_size, margin; preserve_aspect_ratio=true)
        else
            _layout_from_points(stress(graph), plot_size, margin; preserve_aspect_ratio=true)
        end

    gsave()
    translate(center)
    sethue("black")
    drawgraph(
        graph,
        vertexlabels=vertex_labels,
        vertexshapesizes=11,
        vertexlabelfontsizes=10,
        vertexfillcolors=vertex_colors,
        edgestrokeweights=1.6,
        layout=layout_points,
    )
    grestore()
end

"""
    plot_unweighted_search_report(target_graph, target_boundary, save_path; ...)

Create a report-style overview of the `search_unweighted_gadgets` pipeline used
in PR1. The figure explains, step by step, how the target graph is expanded into
equivalent representations, deduplicated, cached, compared against candidate
graphs, and finally returned as `UnweightedGadget` results.
"""
function plot_unweighted_search_report(
    target_graph::SimpleGraph{Int},
    target_boundary::Vector{Int},
    save_path::String;
    max_added_vertices::Int=1,
    preserve_boundary_roles::Bool=true,
    prefilter::Bool=true,
    sample_candidate::Union{Nothing, SimpleGraph{Int}}=nothing,
    sample_candidate_boundary::Union{Nothing, Vector{Int}}=nothing,
    target_pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing,
    sample_candidate_pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}=nothing,
)
    max_added_vertices >= 0 || throw(ArgumentError("max_added_vertices must be non-negative"))

    reprs = equivalent_representations(
        target_graph,
        target_boundary;
        max_added_vertices=max_added_vertices,
        preserve_boundary_roles=preserve_boundary_roles,
    )
    structural_variants = _subdivision_graph_variants(target_graph, max_added_vertices)
    seed_count = sum(
        length(
            _boundary_equivalent_representations(
                g,
                target_boundary;
                preserve_boundary_roles=preserve_boundary_roles,
            ),
        ) for g in structural_variants
    )
    apply_prefilter = prefilter && !_allows_unpinned_components(reprs)
    target_data = GadgetSearch._build_target_data(reprs; include_logical_flips=true)
    flipped_target_count = count(td -> !isempty(td.flip_mask), target_data)

    pattern_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(target_graph, target_boundary))))
    pattern_mask = inf_mask(pattern_reduced)

    sample_graph =
        if !isnothing(sample_candidate)
            sample_candidate
        elseif max_added_vertices > 0
            first(filter(g -> nv(g) > nv(target_graph), structural_variants))
        else
            copy(target_graph)
        end
    sample_boundary = something(sample_candidate_boundary, copy(target_boundary))
    sample_valid, sample_offset = is_gadget_replacement(
        target_graph,
        sample_graph,
        target_boundary,
        sample_boundary,
    )
    sample_reduced = vec(Float64.(content.(calculate_reduced_alpha_tensor(sample_graph, sample_boundary))))
    sample_mask = inf_mask(sample_reduced)

    canvas_width = 1300
    title_height = 84
    box_width = 1120
    box_height = 205
    box_gap = 26
    box_count = 6
    canvas_height = title_height + box_count * box_height + (box_count + 1) * box_gap

    step_titles = [
        "1. Input target graph and boundary",
        "2. Expand into equivalent seeds",
        "3. Deduplicate by reduced alpha tensor",
        "4. Build the reusable filter",
        "5. Scan candidate graphs from the loader",
        "6. Return an UnweightedGadget match",
    ]

    step_lines = [
        [
            "Call site: search_unweighted_gadgets(target_graph, target_boundary, loader)",
            "Input target R has |V|=$(nv(target_graph)), |E|=$(ne(target_graph)), boundary size k=$(length(target_boundary)).",
            "The search API is now single-target: the caller passes one graph and one boundary.",
            "Boundary vertices are highlighted in orange in the thumbnail on the right; preserve_boundary_roles = $(preserve_boundary_roles).",
        ],
        [
            "Internally call equivalent_representations(target_graph, target_boundary; max_added_vertices=...).",
            "Distribute the added-vertex budget across the target edges and subdivide those edges into longer paths.",
            "Example here: $(seed_count) raw seeds are generated before any deduplication.",
            preserve_boundary_roles ?
            "So a crossing arm like 1-3 may become 1-5-3 while the boundary roles stay fixed." :
            "So a crossing arm like 1-3 may become 1-5-3, with optional boundary renumbering also explored.",
        ],
        [
            "For every seed, compute calculate_reduced_alpha_tensor(graph, boundary).",
            "Seeds with the same flattened reduced tensor are considered equivalent and only one is kept.",
            "Example here: $(seed_count) seeds collapse to $(length(reprs)) distinct representations.",
            "Pattern inf-mask example: 0x$(string(pattern_mask; base=16)).",
        ],
        [
            "_make_unweighted_filter(reprs) precomputes per-representation cache entries.",
            "Each entry stores: reduced tensor, inf mask, offset_from_pattern, and logical-flip metadata.",
            "Logical flips are applied at the tensor level; this example generated $(flipped_target_count) flipped target tensors.",
            "Prefilter active in this example: $(apply_prefilter).",
            "If no representation allows an unpinned component, the pins_prefilter can reject candidates early.",
        ],
        [
            "For each candidate graph from the GraphLoader, choose vertex_pool = pinset or 1:nv(candidate).",
            "If the prefilter passes, iterate over all k-subsets of vertex_pool as boundary choices.",
            "For each candidate boundary, compare inf mask first, then test is_diff_by_constant(candidate, target_repr).",
            "A candidate matches as soon as one structural or flip-aware target representation agrees up to a constant offset.",
        ],
        [
            sample_valid ?
            "Sample candidate shown here is a valid replacement: constant_offset = $(round(sample_offset; digits=2))." :
            "Sample candidate shown here does not match the target under the provided boundary.",
            "Returned fields: pattern_graph, replacement_graph, boundary_vertices, constant_offset, pos.",
            "Sample candidate mask: 0x$(string(sample_mask; base=16)), boundary = $(sample_boundary).",
            "This is the object returned to the caller and later used for inspection or plotting.",
        ],
    ]

    _with_drawing(save_path, canvas_width, canvas_height) do
        background("white")

        sethue("black")
        fontsize(24)
        text(
            "PR1 Unweighted Search Flow Report",
            Point(-canvas_width / 2 + 34, -canvas_height / 2 + 28);
            halign=:left,
            valign=:middle,
        )

        fontsize(12)
        text(
            "This figure explains how search_unweighted_gadgets expands a canonical crossing target into subdivided and logical-flip-aware representations with fixed boundary roles by default.",
            Point(-canvas_width / 2 + 34, -canvas_height / 2 + 54);
            halign=:left,
            valign=:middle,
        )

        centers = Point[]
        for step_idx in 1:box_count
            center_y = -canvas_height / 2 + title_height + box_gap + box_height / 2 + (step_idx - 1) * (box_height + box_gap)
            push!(centers, Point(0, center_y))
        end

        for (idx, center) in enumerate(centers)
            _draw_report_box(center, box_width, box_height, step_titles[idx])
            left = center.x - box_width / 2 + 18
            top = center.y - box_height / 2 + 54
            _draw_report_lines(Point(left, top), step_lines[idx]; font_size=11, line_height=16)
        end

        for idx in 1:(length(centers) - 1)
            start_point = Point(0, centers[idx].y + box_height / 2)
            end_point = Point(0, centers[idx + 1].y - box_height / 2)
            _draw_report_arrow(start_point, end_point)
        end

        right_x = canvas_width / 2 - 205
        small_left_x = canvas_width / 2 - 280
        small_right_x = canvas_width / 2 - 120

        _draw_graph_thumbnail(
            target_graph,
            target_boundary,
            Point(right_x, centers[1].y + 12),
            185;
            crossing_geometry=preserve_boundary_roles &&
                _can_draw_with_crossing_geometry(target_graph, target_boundary),
        )

        _draw_graph_thumbnail(
            target_graph,
            target_boundary,
            Point(small_left_x, centers[2].y + 12),
            120;
            crossing_geometry=preserve_boundary_roles &&
                _can_draw_with_crossing_geometry(target_graph, target_boundary),
        )
        if max_added_vertices > 0
            expanded_graph = first(filter(g -> nv(g) > nv(target_graph), structural_variants))
            _draw_graph_thumbnail(
                expanded_graph,
                target_boundary,
                Point(small_right_x, centers[2].y + 12),
                120;
                crossing_geometry=preserve_boundary_roles &&
                    _can_draw_with_crossing_geometry(expanded_graph, target_boundary),
            )
        end
        fontsize(10)
        sethue(RGB(0.33, 0.37, 0.44))
        text("original seed", Point(small_left_x, centers[2].y + 82); halign=:center, valign=:middle)
        if max_added_vertices > 0
            text("seed with one subdivided edge", Point(small_right_x, centers[2].y + 82); halign=:center, valign=:middle)
        end

        repr_preview = reprs[1]
        _draw_graph_thumbnail(
            repr_preview[1],
            repr_preview[2],
            Point(right_x, centers[3].y + 8),
            150;
            crossing_geometry=preserve_boundary_roles &&
                _can_draw_with_crossing_geometry(repr_preview[1], repr_preview[2]),
        )

        sethue(RGB(0.95, 0.96, 0.99))
        box(Point(right_x, centers[4].y + 8), 210, 92, :fill)
        sethue(RGB(0.76, 0.82, 0.92))
        box(Point(right_x, centers[4].y + 8), 210, 92, :stroke)
        _draw_report_lines(
            Point(right_x - 92, centers[4].y - 18),
            [
                "cached tuple",
                "reduced = vector(...)",
                "mask = 0x$(string(pattern_mask; base=16))",
                "flip metadata included",
            ];
            font_size=10,
            line_height=14,
        )

        _draw_graph_thumbnail(
            sample_graph,
            sample_boundary,
            Point(right_x, centers[5].y + 8),
            170;
            pos=sample_candidate_pos,
            crossing_geometry=preserve_boundary_roles &&
                _can_draw_with_crossing_geometry(sample_graph, sample_boundary),
        )

        result_border = sample_valid ? RGB(0.53, 0.75, 0.56) : RGB(0.86, 0.55, 0.55)
        result_fill = sample_valid ? RGB(0.94, 0.99, 0.94) : RGB(1.0, 0.96, 0.96)
        sethue(result_fill)
        box(Point(right_x, centers[6].y + 8), 220, 110, :fill)
        sethue(result_border)
        box(Point(right_x, centers[6].y + 8), 220, 110, :stroke)
        _draw_report_lines(
            Point(right_x - 92, centers[6].y - 26),
            [
                sample_valid ? "status = MATCH" : "status = NO MATCH",
                "boundary_vertices = $(sample_boundary)",
                "constant_offset = $(round(sample_offset; digits=2))",
                "replacement |V| = $(nv(sample_graph))",
            ];
            font_size=10,
            line_height=15,
        )
    end

    println("Drawing saved as $save_path")
    return save_path
end
