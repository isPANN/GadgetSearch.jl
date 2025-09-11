struct Gadget{T<:Real}
    ground_states::BitMatrix
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::Vector{T}
    pos::Union{Nothing, Vector{Tuple{Float64, Float64}}}
end


function save_results_to_json(results::Vector{Gadget}, file_path::String)
    # Convert each result to a serializable dictionary
    json_results = map(results) do res
        base_dict = Dict(
            "ground_states" => [vec(res.ground_states[i, :]) for i in 1:size(res.ground_states, 1)],
            "pins" => res.pins,
            "graph" => Dict(
                "nodes" => [Dict("id" => i, "weight" => res.weights[i]) for i in 1:length(res.weights)],
                "edges" => [Dict("source" => src(e), "target" => dst(e)) for e in Graphs.edges(res.graph)],
            )
        )

        # Add position data if it's a grid_gadget
        if !isnothing(res.pos)
            base_dict["graph"]["positions"] = [Dict("id" => i, "position" => [pos[1], pos[2]]) for (i, pos) in enumerate(res.pos)]
        end

        return base_dict
    end

    # Write to JSON file with pretty formatting
    open(file_path, "w") do io
        write(io, JSON3.write(json_results; pretty=true))
    end

    return file_path
end

"""
    check_gadget(gadget::Gadget; _return_info::Bool=false)

Validate a `Gadget` by enumerating its maximal independent sets (MIS) and
reporting the maximum energy among them and the corresponding pin assignments.

# Arguments
- `gadget::Gadget`: The gadget to check.

# Keyword Arguments
- `_return_info::Bool=false`: If `true`, return a string; otherwise log with `@info`.

# Returns
- When `_return_info` is `true`, returns a formatted `String`; otherwise returns `nothing`.
"""
function check_gadget(gadget::Gadget; _return_info::Bool=false)
    g = gadget.graph
    weights = gadget.weights
    pins = gadget.pins

    # Basic validations
    num_vertices = nv(g)
    length(weights) == num_vertices || throw(ArgumentError("length(weights)=$(length(weights)) must equal nv(g)=$num_vertices"))
    all(1 .<= pins .<= num_vertices) || throw(ArgumentError("pins must be within 1:nv(g)"))

    # Enumerate MIS
    mis_result, mis_count = find_maximal_independent_sets(g)
    mis_count == length(mis_result) || (mis_count = length(mis_result))
    mis_count > 0 || return _return_info ? "No maximal independent sets found." : (@info "No maximal independent sets found."; nothing)

    # Helper to accumulate weights of set bits in config
    @inline function _energy_of_config(config::Unsigned)
        total = zero(float(eltype(weights)))
        v = 1
        while v <= num_vertices
            if ((config >> (v - 1)) & 0x1) == 1
                total += weights[v]
            end
            v += 1
        end
        return total
    end

    # Compute energies
    WeightFloat = float(eltype(weights))
    energy_values = Vector{WeightFloat}(undef, mis_count)
    @inbounds for i in 1:mis_count
        energy_values[i] = _energy_of_config(mis_result[i])
    end

    # Find maxima
    max_energy = maximum(energy_values)
    max_indices = findall(==(max_energy), energy_values)

    # Format report
    lines = ["Max energy value: $(max_energy)"]
    for idx in max_indices
        config = mis_result[idx]
        pin_values = @inbounds Int[((config >> (p - 1)) & 0x1) for p in pins]
        push!(lines, "MaximalIS index=$(idx), pins=$(pin_values)")
    end

    msg = join(lines, "\n")
    return _return_info ? msg : (@info msg; nothing)
end