# Flow journal — dynamic-cross-search

**GOAL:** Build a dynamic lattice-native search that reproducibly finds a usable unweighted CROSS on KSG or the triangular lattice without modifying the existing reduced-alpha verifier.
**Success test:** With a recorded random seed and bounded evaluation budget, the public search returns a concrete lattice coordinate witness that passes `is_gadget_replacement` and all required four-port connection checks.
**Started:** 2026-08-08
**KB:** none

## Levers & facts (initial)

- Facts: the verifier is fixed; arbitrary-graph search is out of scope; fixed-grid subset enumeration is too large; random patches and rays waste most evaluations before topology matters.
- Initial distance estimate: high — correctness checks exist, but the proposal distribution does not construct useful frames and no valid CROSS has been rediscovered.

## Trail

### Trial 1 — analyze — level 0
- **Action:** Analyze the current random patch + random pin + random ray state.
- **Outcome:** Conflict: geometry variables dominate the state count while most combinations are immediately unusable; tensor work is spent on candidates that cannot become a connected crossing tile.
- **Distance:** high (was high) → no_progress = 1
- **Note (learned clause):** `{random pins, random rays, post-hoc frame filtering} ⇒ excessive redundant states and no reproducible CROSS`.

### Trial 2 — simulate — level 1
- **Action:** Hold a legal four-arm KSG frame fixed and dynamically add/remove only interior sites.
- **Outcome:** The regular frame starts with zero finite-offset spread, but the span-2 run exhausts 866 distinct legal connected patches at six mask mismatches and then cycles.
- **Distance:** high (was high) → no_progress = 2
- **Note (learned clause):** `{one fixed symmetric frame, single-site add/remove, elite-only retention} ⇒ a small closed basin; legal geometry alone does not provide enough topology or escape moves`.

### Trial 3 — what-if — level 1
- **Action:** Let the four legal KSG arms vary independently, then evolve only frame-preserving lattice edits.
- **Outcome:** Strong progress: all five tested seeds found a complete KSG CROSS within 1,000 evaluations; seed 2 first succeeds at evaluation 250 and returns a 13-site witness.
- **Distance:** low (was high) → no_progress = 0
- **Note (learned clause):** `{legal frame first, independent arm lengths, frame-preserving edits} ⇒ reproducible CROSS discovery without a fixed-grid subset scan`.

### Trial 4 — final-check — level 0
- **Action:** Run the public KSG search with seed 2 and a 400-evaluation budget, then reconstruct and independently check the returned graph.
- **Outcome:** Evaluation 250 yields a 13-site, 24-edge lattice patch; `is_gadget_replacement` returns `(true, 3.0)` and G1–G4 are all true. The full test suite passes.
- **Distance:** solved (was low)

## Notes store (learned clauses, deduplicated)

- `{random pins, random rays, post-hoc frame filtering} ⇒ excessive redundant states and no reproducible CROSS`.
- `{one fixed symmetric frame, single-site add/remove, elite-only retention} ⇒ a small closed basin; legal geometry alone does not provide enough topology or escape moves`.
- `{legal frame first, independent arm lengths, frame-preserving edits} ⇒ reproducible CROSS discovery without a fixed-grid subset scan`.

## Outcome

- **Status:** SOLVED
- **Result:** A bounded, reproducible KSG CROSS search now succeeds without changing the verifier.
- **Reasoning trail (clean):** Random port geometry wasted the budget; one symmetric frame trapped the search; independent legal arms plus frame-preserving edits exposed a short path to a compact verified crossing.
