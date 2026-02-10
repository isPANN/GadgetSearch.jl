# Feature: Add `RydbergUnweightedModel` for Unweighted MIS Gadget Search

## Background

Currently GadgetSearch.jl supports two energy models:
- `RydbergModel` — weighted MIS, optimizes vertex weights h_i via ILP
- `QUBOModel` — weighted QUBO, optimizes both vertex and edge weights

Both require an optimization solver (e.g., HiGHS).

Based on the paper:

> Liu et al., *"Computer-assisted gadget design and problem reduction of unweighted maximum independent set"*

We propose adding `RydbergUnweightedModel`, which searches for gadgets where **all vertex weights are uniformly 1**. Ground states are the Maximum Independent Sets (largest cardinality MIS). No optimizer is needed — gadgets are found purely by graph structure.

## Key Difference from `RydbergModel`

| | `RydbergModel` (Weighted) | `RydbergUnweightedModel` |
|---|---|---|
| **Weights** | Optimized h_i ≥ 1 (integer) | Fixed h_i = 1 for all i |
| **Energy** | E(σ) = -Σ h_i σ_i | E(σ) = -Σ σ_i = -\|S\| |
| **Ground states** | Depend on weight assignment | Maximum Independent Sets |
| **Optimizer** | Required (HiGHS, etc.) | Not needed |
| **Search method** | ILP optimization | Feasibility check (popcount comparison) |

## Implementation Plan

### Core changes in `src/core/search.jl`:
- [x] Define `RydbergUnweightedModel <: EnergyModel`
- [x] Implement `get_state_space(::Type{RydbergUnweightedModel}, graph)` — reuses MIS enumeration
- [x] Implement `_find_weights(::Type{RydbergUnweightedModel}, ...)` — checks if uniform weights satisfy constraint
- [x] Update `solve_weights` to allow `optimizer=nothing` for this model
- [x] Update `make_filter` to handle `RydbergUnweightedModel` in Gadget construction
- [x] Update `search_gadgets` — optimizer is now optional (defaults to `nothing`)
- [x] Add `Gadget` constructor for `RydbergUnweightedModel`

### Module exports in `src/GadgetSearch.jl`:
- [x] Export `RydbergUnweightedModel`
- [x] Export `check_gadget_unweighted`

### Validation in `src/utils/gadget.jl`:
- [x] Add `check_gadget_unweighted` convenience function
- [x] Update `check_gadget` model name display

### Tests in `test/core/search_unweighted.jl`:
- [x] Type system tests
- [x] Gadget constructor tests
- [x] `_find_weights` feasibility check (pass and fail cases)
- [x] `solve_weights` without optimizer
- [x] `search_gadgets` integration test
- [x] `check_gadget_unweighted` validation
- [x] Error handling (other models still require optimizer)

### Example in `examples/triangular_unweighted_example.jl`:
- [x] Complete working example

## Usage Example

```julia
using GadgetSearch

# No optimizer needed!
results, failed = search_gadgets(
    RydbergUnweightedModel,
    loader, 
    constraints;
    pin_candidates=[[1,2,3], [1,2,4]],
    max_result_num=5
)
```

## Notes

- Unweighted gadgets are **much harder to find** than weighted ones, because the constraint is stricter (no weight optimization freedom)
- Larger lattice sizes may be needed to find valid gadgets
- The implementation is backward compatible — existing `RydbergModel` and `QUBOModel` usage is unchanged
- The `optimizer` keyword in `search_gadgets` now defaults to `nothing` but is still validated for models that require it

## References

- Liu, J.-G., Wurtz, J., Nguyen, M.-T., Lukin, M.D., Pichler, H., Wang, S.-T., *"Computer-assisted gadget design and problem reduction of unweighted maximum independent set"*
- Pan, X.-W., Zhou, H.-H., Lu, Y.-M., Liu, J.-G., *"Encoding computationally hard problems in triangular Rydberg atom arrays"*, arXiv:2510.25249
