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

function _map_and_scale(gadget::Union{grid_gadget, unweighted_grid_gadget}, x_plot_size, y_plot_size, margin, preserve_aspect_ratio)
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
        @info "Using continuous color scheme."
        min_w, max_w = minimum(weights), maximum(weights)
        # generate colormap
        cmap = get(continuous_color_scheme, range(0, 1, length=length(weights)))
        # assign colors based on weights
        node_colors = [cmap[floor(Int, (w - min_w) / (max_w - min_w) * (length(cmap) - 1) + 1)] for w in weights]
    else
        @info "Using discrete color scheme."
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

#TODO: use Luxor to add colorbar
function plot_single_gadget(gadget::grid_gadget, save_path::String; plot_size=400, margin=20, preserve_aspect_ratio=true, discrete_color_scheme=ColorSchemes.seaborn_bright, continuous_color_scheme=ColorSchemes.viridis)
    x_new_vals, y_new_vals, x_scale_factor, y_scale_factor = _map_and_scale(gadget, plot_size, plot_size, margin, preserve_aspect_ratio)
    
    @pdf begin
        background("white")
        pts = [Point(x*x_scale_factor, y*y_scale_factor) for (x,y) in zip(x_new_vals, y_new_vals)]
        node_colors = _generate_vertex_color(gadget.weights, discrete_color_scheme, continuous_color_scheme)
        mask = _generate_mask(gadget.pins, nv(gadget.graph))

        drawgraph(gadget.graph, 
                  vertexlabels=mask,  
                  vertexshapesizes=10, 
                  vertexlabelfontsizes=15, 
                  layout=pts, 
                  edgestrokeweights=1,
                  vertexfillcolors=node_colors) 
        
    end plot_size plot_size save_path

    finish()
    println("Drawing saved as $save_path")
end