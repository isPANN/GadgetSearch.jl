# Flow journal — dynamic-cross-search

**GOAL:** Reproducibly find concrete unweighted CROSS embeddings on KSG and the triangular lattice without modifying the reduced-alpha verifier.
**Success test:** A recorded seed and bounded public search return lattice coordinates that pass `is_gadget_replacement` and G1-G4.
**Started:** 2026-08-08

## Levers & facts

- The verifier is fixed; outputs must be induced lattice patches; random frames and fixed-canvas subsets waste most evaluations.
- KSG has a compact solution; the triangular positive control has 37 sites and a non-radial track-swap frame.

## Trail

### Trial 1 — analyze — level 0
- **Action/outcome:** Random pins and rays were rejected post hoc; almost all tensor work went to unusable geometry. `{random frame + post-filter} ⇒ redundant dead states`.

### Trial 2 — simulate/what-if — level 1
- **Action/outcome:** One symmetric KSG frame cycled through 866 states; independent legal arm lengths plus frame-preserving edits found CROSS for 5/5 seeds. `{one symmetric frame} ⇒ basin`; `{independent legal arms} ⇒ escape`.

### Trial 3 — final-check — level 0
- **Action/outcome:** KSG seed 2 succeeds at evaluation 250 with 13 sites, 24 edges, offset 3; verifier and G1-G4 pass, as does the full test suite.

### Trial 4 — analyze — level 0
- **Action/outcome:** Recovered the paper's coordinate-drawn 37-site triangular witness; verifier returns `(true, 15.0)` and G1-G4 pass. Rays 180°, 60°, 60°, 0° prove `{radial-only frames} ⇒ excludes known solutions`.

### Trial 5 — simulate — level 1
- **Action/outcome:** Three 5,000-evaluation triangular runs with 45-site capacity retained only 16-18 sites. Removing early size pressure and seeding 25-40-site track-swap interiors reached one mask mismatch but never exact. `{compactness before correctness} ⇒ premature collapse`.

### Trial 6 — what-if — level 1
- **Action/outcome:** Raw-alpha dominance margins distinguished which state was one unit short; multi-site region rewrites moved the error to other states but did not remove it in 20,000 evaluations. `{local lattice rewrites} ⇒ one-bit plateau`.

### Trial 7 — analyze/backjump — level 1
- **Action/outcome:** Reducing the 37-site witness to an abstract core discarded lattice coordinates, pin geometry, and immediate embeddability. Subsequent work became specific to one topology. `{abstract core first} ⇒ wrong state representation`; reject this branch.

### Trial 8 — simulate/backjump — level 2
- **Action/outcome:** Abstract topology search could satisfy the tensor while failing to produce a concrete triangular or KSG embedding. Encoding those states as opaque graph strings made the trajectory unreadable. The experiments and search-facing encoding were deleted.

### Trial 9 — simulate — level 1
- **Action/outcome:** Lattice-native cluster rewrites and aligned parent crossover preserved concrete pin frames but still stopped one boundary state short. The next move must combine lattice regions while keeping every intermediate state embedded; it must not reintroduce an abstract-graph stage.

## Notes store

- Keep coordinates, ordered pins, and rays as the complete search state; construct legal geometry rather than filtering it; retain non-radial frames; treat repeated one-bit plateaus as a move-set conflict, not non-existence.

## Outcome

- **Status:** in progress
- **Result:** KSG is solved. The triangular positive control is independently verified. The abstract-core detour was rejected and removed; triangular search remains one boundary state short.
- **Next lever:** add lattice-native region recombination that never leaves concrete KSG or triangular coordinates.
