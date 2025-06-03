function generic_rule(rule_id::Int, bit_num::Tuple{Int, Int}; show_info::Bool=false)::BitMatrix
    @assert length(bit_num) == 2
    input_bits = bit_num[1]
    output_bits = bit_num[2]
    total_bits = input_bits + output_bits
    num_inputs = 2^input_bits
    num_outputs = 2^output_bits
    max_gateid = num_outputs^num_inputs

    if rule_id < 0 || rule_id >= max_gateid
        error("Gate ID must be between 0 and $(max_gateid - 1)")
    end

    result = BitMatrix(undef, num_inputs, total_bits)
    mask_output = num_outputs - 1
    shift_amount = output_bits

    for input in 0:(num_inputs - 1)
        output = (rule_id >> (input * shift_amount)) & mask_output
        combined = (input << shift_amount) | output
        for b in 1:total_bits
            result[input + 1, b] = ((combined >> (total_bits - b)) & 1) == 1
        end
    end

    show_info && show_rule_info(rule_id, bit_num)
    return result
end


function show_rule_info(rule_id::Int, bit_num::Tuple{Int, Int})
    @assert length(bit_num) == 2
    input_bits = bit_num[1]; output_bits = bit_num[2]
    num_inputs = 2 ^ input_bits
    num_outputs = 2 ^ output_bits
    max_rule_id = num_outputs ^ num_inputs

    if rule_id < 0 || rule_id > max_rule_id
        @error("Gate ID must be between 0 and $max_rule_id")
    end
    @info "==== Gate ID: $rule_id ===="

    # Create the header for the table
    header = "Input (Binary) | Output (Binary)"
    separator = "-" ^ length(header)
    println(header)
    println(separator)

    for input in 0:(num_inputs - 1)
        # Calculate output using bitwise operations
        output = (rule_id >> (input * output_bits)) & (num_outputs - 1)

        # Format input and output to fixed-width binary strings
        input_bin = string(input, base=2, pad=input_bits)
        output_bin = string(output, base=2, pad=output_bits)

        # Print the formatted result
        println("$input_bin | $output_bin")
    end
end