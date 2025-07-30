# Type definitions for GadgetSearch

"""
    SearchParameters(;
        pin_set::Vector{Int} = Int[],
        max_file_size_mb::Int = 30,
        split_size::Int = 700_000,
        start_idx::Int = 0,
        end_idx::Int = 0,
        greedy::Bool = false,
        threshold::Int = 0,
        max_samples::Int = 0
    )

Contains parameters for graph search operations

# Fields
- `pin_set::Vector{Int}`: Set of pins to search for
- `max_file_size_mb::Int`: Maximum file size in MB to process before splitting
- `split_size::Int`: Maximum number of rows in each split file
- `start_idx::Int`: Starting index for graph search
- `end_idx::Int`: Ending index for graph search
- `greedy::Bool`: Whether to use greedy approach
- `threshold::Int`: Maximum number of MIS combinations to consider
- `max_samples::Int`: Maximum number of samples to generate
"""
struct SearchParameters
    pin_set::Vector{Int}
    max_file_size_mb::Int
    split_size::Int
    start_idx::Int
    end_idx::Int
    greedy::Bool
    threshold::Int
    max_samples::Int

    # Constructor with default values and validation
    function SearchParameters(;
        pin_set::Vector{Int} = Int[],
        max_file_size_mb::Int = 30,
        split_size::Int = 700_000,
        start_idx::Int = 0,
        end_idx::Int = 0,
        greedy::Bool = false,
        threshold::Int = 0,
        max_samples::Int = 0
    )
        # Parameter validation
        max_file_size_mb >= 0 || throw(ArgumentError("max_file_size_mb must be positive"))
        split_size > 0 || throw(ArgumentError("split_size must be positive"))
        start_idx >= 0 || throw(ArgumentError("start_idx must be non-negative"))
        end_idx >= 0 || throw(ArgumentError("end_idx must be non-negative"))
        threshold >= 0 || throw(ArgumentError("threshold must be non-negative"))
        max_samples >= 0 || throw(ArgumentError("max_samples must be non-negative"))

        new(pin_set, max_file_size_mb, split_size, start_idx, end_idx, greedy, threshold, max_samples)
    end
end
