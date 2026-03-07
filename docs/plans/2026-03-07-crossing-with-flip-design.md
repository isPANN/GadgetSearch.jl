# Crossing Gadget Search with Logical Flip

**Date:** 2026-03-07
**Goal:** Find crossing gadgets by searching for "crossing + logical flip" combinations

## Background

When searching for crossing gadgets on triangular lattice UDGs, we can expand the search space by looking for gadgets that implement "crossing + NOT gate on one or more boundary pins". If we find such a gadget, we can construct a pure crossing gadget by adding external NOT gadgets (which are already known) to cancel the flips.

**Key insight:** If a gadget implements "crossing + flip", then crossing itself must be realizable.

## Design

### 1. Target Generation

Generate multiple target graphs by:

1. **Flip variants**: Apply logical flips to boundary pins of the CROSS graph
   - Flip = inverting a pin's logic value in the reduced alpha tensor
   - In tensor terms: `α̃_flipped[..., i_j, ...] = α̃_original[..., 1-i_j, ...]`
   - Consider symmetry to avoid equivalent configurations

2. **Extended variants**: Insert nodes on edges (copy gadget pattern)
   - Example: edge 1-3 becomes path 1-5-3
   - Allows controlling reduced alpha tensor offset

3. **Combination**: Generate flip variants for both base CROSS and extended versions

### 2. Symmetry Analysis

CROSS graph has symmetry group that makes some flip combinations equivalent:
- Rotation symmetry (90°, 180°, 270°)
- Reflection symmetry (horizontal, vertical, diagonals)

We only search non-equivalent flip patterns to avoid redundant computation.

### 3. Search Strategy

Use existing `search_multi_target_gadgets` infrastructure:
- Pre-compute reduced alpha tensors for all target variants
- Single pass over candidate graphs
- Each candidate's α̃ computed once, compared against all targets via inf-mask fingerprinting

### 4. Visualization

Generate Typst document with:
- Summary table of all found gadgets
- Graph visualizations showing structure and boundary pins
- Metadata: which target matched, vertex/edge counts, constant offset
- Output to desktop as PDF

## Implementation Plan

**Files to modify:**
1. New file: `src/utils/flip_variants.jl` - target generation functions
2. New file: `examples/crossing_with_flip.jl` - search script
3. New file: `compile_typst.jl` - visualization (already exists, may need updates)

**Key functions:**
- `flip_alpha_tensor(tensor, pins_to_flip)` - flip specified dimensions
- `generate_flip_variants(graph, boundary)` - generate non-equivalent flips
- `generate_extended_cross()` - create extended CROSS variants
- `generate_all_crossing_targets()` - combine all variants
- `visualize_results_typst(results, output_path)` - create Typst report

## Expected Outcomes

- Comprehensive search covering CROSS + all non-equivalent flip combinations
- Increased chance of finding crossing gadgets
- Clear visualization of results for analysis
- Reusable framework for future gadget searches with logical operations
