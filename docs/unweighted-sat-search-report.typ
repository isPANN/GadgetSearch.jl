#set document(
  title: "Certified Search and Rule-Based Optimization for Unweighted MIS Gadgets",
  author: "GadgetSearch",
)
#set page(
  paper: "a4",
  margin: (left: 21mm, right: 21mm, top: 15mm, bottom: 17mm),
  numbering: "1",
  number-align: center + bottom,
)
#set text(font: "Libertinus Serif", size: 9.25pt)
#set math.equation(numbering: "(1)")
#set par(justify: true, leading: 0.58em)
#set heading(numbering: "1.1")
#set table(stroke: 0.45pt + rgb("#aeb7bf"), inset: (x: 5pt, y: 3.5pt))
#show heading.where(level: 1): it => block(
  above: 12pt, below: 5pt,
  text(font: "Libertinus Sans", size: 14pt, weight: "semibold", it),
)
#show heading.where(level: 2): it => block(
  above: 9pt, below: 3pt,
  text(font: "Libertinus Sans", size: 11pt, weight: "semibold", it),
)
#show figure.caption: set text(size: 8.2pt)
#show raw: set text(font: "Libertinus Mono", size: 8pt)

#let gray = rgb("#66727d")
#let fig(name, caption, width: 100%) = figure(
  image("figures/unweighted-sat/" + name + ".svg", width: width),
  caption: caption,
)

#align(center)[
  #text(font: "Libertinus Serif Display", size: 19pt, weight: "semibold")[
    Certified Search and Rule-Based Optimization for Unweighted MIS Gadgets
  ]
  #v(2pt)
  #text(size: 9pt)[GadgetSearch.jl]
]

#v(4pt)
#block(inset: (x: 11mm, y: 7pt), stroke: (
  top: .6pt + gray, bottom: .6pt + gray,
))[
  #text(size: 9pt)[#text(font: "Libertinus Sans", weight: "semibold")[Abstract.]
    We describe a certified synthesis framework for four-pin unweighted
    maximum-independent-set gadgets on the triangular and king's graph
    lattices. A joint SAT formulation discovers a first gadget by choosing its
    occupied sites and four ports in one finite window. A complementary
    frame-first formulation enumerates port geometry and solves a smaller
    occupancy problem for each frame. Explicit independent-set witnesses impose
    lower bounds on the reduced alpha tensor, frontier binary decision diagrams
    impose upper bounds, and direct layered reachability imposes connectivity.
    The second stage treats successful reductions as evidence from which
    to extract reusable, interface-preserving rewrite schemas. Their common
    invariant is equality, up to one constant, of the conditioned MIS signature
    seen across a small separator. Exact local replacements and
    interface-constrained SAT re-synthesis are instances of this rule. Every
    application is certified by the unchanged verifier. A reduction from a
    blank-window 28-site CROSS to 23 sites serves as a regression case, not as
    the definition of the optimization algorithm.]
]

= Problem formulation

Let $R$ be a target graph with four ordered boundary vertices
$partial R=(b_1,b_2,b_3,b_4)$. For a boundary state
$sigma=(sigma_1,...,sigma_4) in {0,1}^4$, the value $sigma_i=1$ requires
$b_i$ to belong to the independent set. Define
$
  alpha(R)_sigma = max{abs(I): I " independent in " R,
  I ∩ partial R = {b_i:sigma_i=1}}.
$
An infeasible boundary state has value $-infinity$. The *reduced alpha tensor*
$tilde(alpha)(R)$ removes entries dominated by a proper subconfiguration.
Denote the target tensor by $T=tilde(alpha)(R)$.

Let $Lambda$ be a lattice and let $H$ be the blockade graph of a finite lattice
window. A selected site set $S subset.eq V(H)$ induces the replacement graph
$G=H[S]$. Four ordered selected sites $P=(p_1,p_2,p_3,p_4)$ represent the
boundary vertices, and lattice directions $bold(r)=(r_1,r_2,r_3,r_4)$ describe
the exterior wires. The search seeks $S$, $P$, $bold(r)$, and an integer offset
$c$ satisfying
$
  tilde(alpha)(G) = T+c.
$

For each pin, let $s_i=p_i+r_i$ be the first interface site outside the gadget.
The frame must satisfy four geometric conditions:

- *(G1)* every $s_i$ is a strict vertex of the convex hull of
  $S union {s_1,s_2,s_3,s_4}$;
- *(G2)* the cyclic hull order alternates the interface pairs $(s_1,s_3)$ and
  $(s_2,s_4)$;
- *(G3)* every $r_i$ points strictly outward from the interface centroid;
- *(G4)* the exterior rays are unobstructed, touch only their own pins, and
  remain pairwise outside blockade range.

The selected graph must contain exactly $N$ sites, contain all four pins, and
be connected. A candidate is accepted only when the unchanged
`is_gadget_replacement` routine independently recomputes the tensor, agrees on
$c$, and the geometry checker confirms G1--G4.

= Two-stage algorithm

The recommended workflow has two stages: *search* first produces one or more
certified seed gadgets, and *optimization* then attempts to reduce each seed.
The public API exposes the two stages separately so that searches and optimizer
budgets can be scheduled independently; a complete synthesis run calls them in
this order.

#fig("pipeline", [The complete algorithm is Stage I search followed by Stage II
optimization. Joint and frame-first SAT are alternative engines inside the
search stage, not alternatives to optimization.])

Stage I has two search engines that share the tensor and connectivity encodings:

- *Joint occupancy-and-frame SAT* chooses $S$, $P$, and $bold(r)$ in one CNF.
  It is the primary blank-window discovery method.
- *Frame-first SAT* enumerates $(P,bold(r))$ geometrically and lets SAT choose
  only $S$. It is suited to systematic finite-window traversal and local
  re-synthesis around a known gadget.

Stage II consumes every verifier-valid `UnweightedGadget` returned by Stage I.
It never assumes that the starting gadget came from a particular solver run.
If no smaller certified descendant exists within the rule and SAT budgets, the
seed itself is the final result of the two-stage algorithm.

= Joint occupancy-and-frame SAT

Fix a lattice window $W$, atom count $N$, and offset $c$. For every site
$v in W$, introduce a selection variable $x_v$. For every geometrically
eligible site-direction pair and boundary label $i$, introduce a frame variable
$f_(i,v,d)$. The joint CNF imposes:

1. $sum_v x_v=N$;
2. exactly one frame choice for each label $i$;
3. each chosen pin is selected and the four pins are distinct;
4. local corridor, outward-direction, and alternating-interface constraints;
5. the reduced-alpha tensor constraints described below;
6. connectivity of the selected induced graph.

Only eligible port choices receive variables. This sparse allocation is much
smaller than allocating every $(i,v,d)$ combination and forbidding most of them.

One reference port is fixed by a canonical translation and ray direction. This
breaks translation and rotation symmetry, but it also defines the represented
finite instance. A known abstract gadget is absent if no symmetry-equivalent
embedding fits the anchored window. Canonical placement is therefore part of
the instance specification, not an innocuous implementation detail.

The formula is emitted as DIMACS and solved by Kissat. A model is materialized
and checked geometrically before the reduced tensor is recomputed. A rejected
model contributes a blocker and the solver is restarted; no rejected candidate
is reported as a gadget.

== Solver portfolios

Blank-window joint instances can be highly seed-sensitive even when they
represent the same finite search space. Independent Kissat seeds are therefore
a useful source of parallelism. A timed-out run is recorded as `UNKNOWN` and
may be rescheduled; it is never interpreted as an UNSAT certificate.

= Frame-first search

Frame-first search moves geometry outside the occupancy CNF. It constructs a
compatibility graph whose vertices are locally valid $(p,r)$ candidates. Four
mutually compatible candidates form a possible four-port frame. Valid label
orders are then checked against the complete G1--G4 conditions and canonicalized
under the lattice symmetries.

#fig("frame-enumeration", [Frame-first construction. Compatible port candidates
form a four-clique; a valid ordered frame determines the allowed host before
occupancy variables are introduced.])

For a surviving frame, each lattice site is tested against its pin corridors
and hull constraints. The allowed set $U(P,bold(r))$ defines the host
$H=G[U]$. The fixed-frame CNF contains site-selection, tensor, and connectivity
variables, but no frame-choice variables.

For the known 23-site triangular CROSS frame, Kissat selects the occupied sites
and returns a verifier-valid gadget at offset 7. The hard part of a blank
frame-first search is reaching this rare frame, not solving its occupancy
instance.

The production frame traversal has deterministic cursors and periodically
writes checkpoints containing the current frame, label order, offset, and
filter counters. Repeating the same bounded search with that checkpoint resumes
the traversal. This makes long finite traversals resumable, although their
total geometry space can still be very large.

= SAT encoding for a fixed host

Fix $(H,P,N,c)$. The formula is satisfiable exactly when $H$ contains a
connected $N$-site induced subgraph with reduced alpha tensor $T+c$.

== Atom selection

For every host vertex $v$, introduce $x_v$, with
$
  x_(p_i)=1 quad (i=1,...,4),
  quad sum_(v in V(H)) x_v=N.
$

== Direct layered connectivity

Let $z_v^t$ mean that selected vertex $v$ is reachable from $p_1$ by a path of
at most $t$ host edges. The initial layer contains only the root,
$
  z_v^0 arrow.l.r.double v=p_1.
$
For $t=0,...,N-2$, direct propagation is
$
  z_v^(t+1) arrow.l.r.double
  (z_v^t or (x_v and or.big_(u in N_H(v)) z_u^t)).
$
Finally,
$
  x_v => z_v^(N-1) quad (v in V(H)).
$
No per-edge arrival variables are used. For a host with $h=abs(V(H))$ vertices
and $m=abs(E(H))$ undirected edges, this encoding contributes $N h$ auxiliary
variables and
$
  2h+(N-1)(3h+2m)
$
clauses. The dependence on host edges is linear, matching the sparse TLSG and
KSG hosts.

#fig("connectivity-layers", [Direct layered reachability for a selected path
$p_1-u-v$ and an isolated host vertex $w$. Reachability persists and propagates
through selected vertices; $w$ cannot satisfy the final implication.], width: 78%)

The depth $N-1$ is exact: every vertex in a connected $N$-vertex graph has a
simple path of length at most $N-1$ from $p_1$, whereas a selected vertex in
another component can never become reachable.

== From the reduced tensor to bounds

For a selected set $S$ and boundary state $sigma$, let $A_S(sigma)$ be the
largest independent-set size with exactly the pins indicated by $sigma$.
Define the monotone completion
$
  C_T(sigma)=max{T_tau: tau subset.eq sigma, T_tau != -infinity}.
$
The equality $tilde(alpha)(G[S])=T+c$ is imposed by
$
  A_S(sigma) >= T_sigma+c quad "for finite " T_sigma,
$
and
$
  A_S(sigma) <= C_T(sigma)+c quad "for every " sigma.
$
Only lower and upper states not implied by other boundary states are emitted.

== Lower bounds: explicit witnesses

For each essential finite $T_sigma$, introduce witness variables
$y_v^sigma$. The clauses impose
$
  y_v^sigma => x_v,
  quad y_u^sigma+y_v^sigma <= 1 quad ({u,v} in E(H)),
$
$
  y_(p_i)^sigma=sigma_i,
  quad sum_(v in V(H)) y_v^sigma=T_sigma+c.
$
Thus the selected graph explicitly contains an independent set of the required
size and boundary state.

== Upper bounds: frontier BDDs

Upper bounds must exclude every oversized independent set. The implementation
constructs a frontier binary decision diagram under a fixed vertex order. A
path either skips or takes the current host vertex. A take edge is enabled only
when the site is selected and no taken frontier vertex conflicts with it.

After each prefix, BDD states with the same active frontier are merged. The
count is capped at $C_T(sigma)+c+1$, and every terminal state at the cap is
forbidden. The implementation tries three lattice projections and their
reversals and retains the smallest of the six BDDs.

#fig("bdd-example", [Frontier BDD for the path $a-b-c$ with upper bound one.
Skip and take branches represent every independent set; the two-vertex terminal
is forbidden.], width: 72%)

For fixed $(H,P,N,c)$, the complete occupancy formula is
$
  Phi = Phi_"select" and Phi_"connect"
  and ∧_(sigma in L) Phi_"witness"^sigma
  and ∧_(sigma in U) Phi_"BDD"^sigma,
$
where $L$ and $U$ are the essential lower and upper boundary states.

= Independent certification

A satisfying assignment is not returned directly. The selected coordinates
are materialized as the exact blockade graph, after which the production
verifier recomputes $tilde(alpha)(G)$ and checks the common offset. Connectivity
and G1--G4 are checked independently.

#fig("verification", [The SAT solver proposes; the unchanged verifier accepts.
The same certification boundary is used during discovery and after every
rewrite.])

This separation protects the result from an incomplete SAT-side optimization:
auxiliary encodings may change, but the semantic acceptance test does not.

= Reusable rewrite principle and implemented optimizer

A successful reduction should not be stored merely as a sequence of edits.
The useful object is the *semantic reason* why the changed region can be
replaced in any compatible context. For MIS gadgets, that reason is a boundary
signature across a small separator.

== Interface signatures

Let $Q$ be a connected subgraph of a certified gadget and let
$C=(c_1,...,c_k)$ contain every vertex through which $Q$ meets the unchanged
exterior. For an interface state $tau in {0,1}^k$, define
$
  A_Q^C(tau)=max{abs(I): I " independent in " Q,
  I ∩ C={c_i:tau_i=1}}.
$
The vector $"Sig"_C(Q)=(A_Q^C(tau))_tau$ is everything the exterior needs to
know about $Q$ when maximizing an independent set. A local replacement
$Q arrow.r Q'$ is context-independent when
$
  A_(Q')^C(tau)=A_Q^C(tau)+delta
  quad "for every feasible " tau,
$
and the infeasible-state pattern is unchanged. Gluing either patch to the same
exterior then changes every global conditioned optimum by the same $delta$.
Consequently, the reduced alpha tensor of the whole gadget changes only by the
allowed constant offset.

The signature condition is necessary but not sufficient for a lattice rewrite.
The replacement must also preserve the attachment coordinates or port roles,
introduce no unintended blockade edge to the exterior, maintain connectivity,
and satisfy G1--G4. The unchanged whole-gadget verifier remains the final
certificate.

== Generalizing one successful case

A successful before/after pair suggests a reusable rule only when its changed
region can be isolated behind a small separator and the two local interface
signatures differ by one constant. Absolute coordinates must then be replaced
by structural premises such as interface order, adjacency, pin role, channel
pairing, and lattice direction. This is the criterion for adding future rules:
the proof obligation is the interface contract, not reproduction of one edit
sequence. The current code does not automatically mine or minimize schemas.

== Operators implemented today

The current optimizer implements two operator families consistent with the
signature principle.

- *Exact local replacement.* Two closed-form rules contract an even boundary
  tail or a valid opposite pair of leaf pins. Each proposed whole gadget is
  certified; the implementation does not yet enumerate arbitrary equivalent
  local patches from their signatures.
- *Interface-constrained re-synthesis.* Hold the separator, external port roles,
  and required signature fixed, then ask SAT to synthesize a smaller interior.
  Moving a pin by one lattice step is one proposal mechanism for changing the
  interface geometry, not the semantic rewrite rule. The rule is to preserve
  the interface contract while re-solving the interior.

The optimizer searches the graph generated by these operators. It explores
all certified direct descendants to a closure and only then spends SAT budget
on re-synthesis. Exact coordinate, pin, and ray state keys merge identical
converging paths. This non-greedy search is essential: the smallest immediate
child may be a dead end, while a different equivalent child exposes a stronger
subsequent rewrite.

== What the CROSS case teaches

The reduction from the independently discovered 28-site CROSS to a 23-site
gadget exercises both implemented operator families: certified structural
contractions and one interface-constrained interior re-synthesis. Its general
lesson is to preserve an interface contract, explore alternative certified
descendants, and invoke exact re-synthesis when a structural step changes the
available interior. The particular atom counts and rule names are regression
evidence, not assumptions of a new search.

The same second stage is applied to CROSS+EDGE. Stage I finds a verified
9-site gadget; Stage II exhausts the direct contraction closure and a budget of
256 fixed-frame SAT evaluations without finding a smaller certified descendant.
The reported 9-site result is therefore the best result under the current rules
and budget, not a search-only result and not a proof of global minimality.

The final result is certified relative to the available rule schemas, host
radius, minimum atom count, and SAT budget. It is not a proof that no smaller
gadget exists.

= Results and regression coverage

#figure(
  table(
    columns: (1.35fr, 1.35fr, 1.35fr, 1.35fr, 2.3fr),
    align: (left, center, center, center, left),
    fill: (x, y) => if y == 0 { rgb("#edf4fa") } else { white },
    table.header(
      [*target*], [*Stage I search*], [*Stage II optimize*],
      [*offset*], [*final status*],
    ),
    [CROSS], [28 sites], [23 sites], [$10 arrow.r 7$],
      [reduced; verified],
    [CROSS+EDGE], [9 sites], [9 sites], [$1 arrow.r 1$],
      [unchanged after budget; verified],
  ),
  caption: [End-to-end results of the complete two-stage algorithm.],
)

Correctness is exercised at several levels:

- direct layered connectivity is compared with explicit connected-subset
  enumeration on small host graphs;
- fixed-frame and joint positive controls pass the unchanged verifier;
- Stage I search covers both CROSS and CROSS+EDGE, and Stage II optimization is
  run on both results;
- the historical 28-site CROSS satisfies the production CNF and verifier;
- a regression starts from that real 28-site result and checks every certified
  rewrite application down to 23 sites;
- the focused unweighted-search suite and the complete package test suite pass.

= Completeness and practical limits

The procedure is complete only within its explicitly scheduled finite space.

- A joint instance fixes one window, $N$, $c$, canonical anchor, and first ray.
- A frame-first run fixes a finite window, atom-count range, frame budget, and
  SAT-call budget.
- Kissat `UNKNOWN` records expired work, not a negative proof.
- The rewrite optimizer explores a bounded signature-preserving rule
  neighborhood and does not
  certify global minimality.

Frame-first enumeration can encounter millions of easy UNSAT frames before a
useful frame, with a small hard tail dominating wall time. Joint SAT avoids that
explicit traversal but exposes one larger, seed-sensitive CNF. The two modes are
therefore complementary rather than competing implementations.

= Recommended workflow

For a new four-pin unweighted target:

1. *Stage I -- search:* compute the target tensor, choose TLSG or KSG, and run
   joint or frame-first SAT until at least one seed gadget is independently
   verified;
2. *Stage II -- optimize:* run exact local replacements and
   interface-constrained re-synthesis on every distinct seed, retaining the
   smallest verified descendant or the unchanged seed at a fixed point;
