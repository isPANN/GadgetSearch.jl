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
    
    scale_factor_x = allowed_x_max / maximum(abs.(xvals))
    scale_factor_y = allowed_y_max / maximum(abs.(yvals))

    if preserve_aspect_ratio
        scale_factor = min(scale_factor_x, scale_factor_y)
        return scale_factor, scale_factor
    end
    return scale_factor_x, scale_factor_y
end

function _map_and_scale(gadget::Gadget{T}, x_plot_size, y_plot_size, margin, preserve_aspect_ratio) where T
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

function plot_single_gadget(gadget::Gadget{T}, save_path::String; plot_size=400, margin=20, preserve_aspect_ratio=true, discrete_color_scheme=ColorSchemes.seaborn_bright, continuous_color_scheme=ColorSchemes.viridis) where T
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
                            gadget::Gadget{T}, save_path::String; 
                            plot_size=400, margin=30, 
                            preserve_aspect_ratio=true, 
                            background_grid=false,
                            show_weights=true,
                            round_weights=false,
                            discrete_color_scheme=ColorSchemes.seaborn_bright, 
                            continuous_color_scheme=ColorSchemes.viridis
                            ) where T

    # Create a new graph with only non-zero weight vertices
    valid_vertices = findall(w -> w != 0, gadget.weights)
    new_graph = SimpleGraph(length(valid_vertices))
    
    # Create mapping from old vertex indices to new ones
    old_to_new = Dict(v => i for (i, v) in enumerate(valid_vertices))
    
    # Add edges between valid vertices
    for e in edges(gadget.graph)
        src, dst = Graphs.src(e), Graphs.dst(e)
        if src in valid_vertices && dst in valid_vertices
            add_edge!(new_graph, old_to_new[src], old_to_new[dst])
        end
    end
    
    # Create new weights and pins arrays
    new_weights = gadget.weights[valid_vertices]
    new_pins = [old_to_new[p] for p in gadget.pins if p in valid_vertices]
    
    # Create new positions array
    new_pos = gadget.pos[valid_vertices]
    
    # Create new gadget with filtered data
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
                  edgestrokeweights=2,
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
