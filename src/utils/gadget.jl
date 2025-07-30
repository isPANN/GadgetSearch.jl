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

function check_gadget(gadget::Gadget)
    g = gadget.graph
    weights = gadget.weights
    pins = gadget.pins

    # Get all maximal independent sets (MIS), assume the result type is a Vector{UInt}
    mis_result, mis_num = find_maximal_independent_sets(g)

    # Calculate energy for each MIS
    energy_value = Vector{Float64}(undef, mis_num)
    for i in 1:mis_num
        config = mis_result[i]
        total = 0.0
        for v in 1:nv(g)
            if ((config >> (v - 1)) & 0x1) == 1
                total += weights[v]
            end
        end
        energy_value[i] = total
    end

    # Check the energy value of each MIS
    @show energy_value
    max_energy = maximum(energy_value)
    max_indices = findall(x -> x == max_energy, energy_value)

    @info "Max energy value: $max_energy"
    for idx in max_indices
        config = mis_result[idx]
        pin_values = [((config >> (p - 1)) & 0x1) for p in pins]
        @info "MIS index: $idx, pins = $pin_values"
    end

    return nothing
end