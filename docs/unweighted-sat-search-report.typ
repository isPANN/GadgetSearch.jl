#set document(
  title: "Unweighted MIS Gadget Search",
  author: "GadgetSearch",
)
#set page(
  paper: "a4",
  margin: (left: 20mm, right: 20mm, top: 16mm, bottom: 18mm),
  numbering: "1",
  number-align: center + bottom,
)
#set text(font: "Libertinus Serif", size: 10pt)
#set math.equation(numbering: "(1)")
#set par(justify: true, leading: 0.65em, spacing: 0.75em)
#set heading(numbering: "1.1")
#set table(stroke: 0.4pt + rgb("#b7c0c8"), inset: (x: 5.5pt, y: 3.6pt))
#show heading.where(level: 1): it => block(
  above: 14pt, below: 6pt,
  text(font: "Libertinus Sans", size: 13.5pt, weight: "semibold", it),
)
#show heading.where(level: 2): it => block(
  above: 10pt, below: 4pt,
  text(font: "Libertinus Sans", size: 11pt, weight: "semibold", it),
)
#show heading.where(level: 3): it => block(
  above: 8pt, below: 3pt,
  text(font: "Libertinus Sans", size: 10pt, weight: "semibold", it),
)
#show figure.caption: set text(size: 8.4pt)
#show raw: set text(font: "Libertinus Mono", size: 8.2pt)
#show math.equation: set text(size: 10.5pt)

#import "@preview/cetz:0.4.2"

#let ink = rgb("#1f2a33")
#let muted = rgb("#5b6770")
#let pin-fill = rgb("#d7e6f4")
#let pin-stroke = rgb("#2f5f86")
#let site-fill = rgb("#f7f4ee")
#let site-stroke = rgb("#4a5560")
#let is-fill = rgb("#f3d48a")
#let accent = rgb("#2f5f86")
#let pale = rgb("#f3f6f8")
#let rule = rgb("#d5dde3")

#let note(title, body) = block(
  width: 100%,
  inset: (x: 9pt, y: 8pt),
  fill: pale,
  stroke: 0.4pt + rule,
  radius: 3pt,
  {
    text(font: "Libertinus Sans", size: 8.6pt, weight: "semibold", fill: accent, title)
    v(3pt)
    set text(size: 9.2pt)
    body
  },
)

#let cetz-canvas = cetz.canvas.with(length: 1cm)
#let node-at(pos, label, fill: site-fill, paint: site-stroke, r: 0.22) = {
  import cetz.draw: circle, content
  circle(pos, radius: r, fill: fill, stroke: 0.85pt + paint)
  content(pos, text(size: 7.6pt, weight: "semibold", label))
}

#align(center)[
  #text(font: "Libertinus Serif Display", size: 18.5pt, weight: "semibold")[
    Finding Unweighted MIS Gadgets by SAT Search
  ]
  #v(4pt)
  #text(size: 9.4pt, fill: muted)[
    A walkthrough of SAT-based gadget search, with CROSS as the running example
  ]
]

#v(6pt)

This note explains the two-stage procedure that searches a finite lattice window
for an unweighted maximum-independent-set (MIS) gadget and then shrinks it. The
aim is to make the search problem, the SAT encoding, and the rewrite rules
readable, using one target graph throughout.

= The search problem

Let $R$ be a four-pin target with ordered boundary
$partial R = (b_1,b_2,b_3,b_4)$. Write $alpha(R)_sigma$ for the largest
independent-set size compatible with boundary state $sigma in {0,1}^4$, and
$tilde(alpha)(R)$ for its reduced tensor. The search target is
$T = tilde(alpha)(R)$. A graph $G$ with the same ordered pins is a valid
replacement when

$
  tilde(alpha)(G) = T + c
$

for a constant $c$.
The running example is CROSS (@fig-problem~(a)).

== Lattice embedding

A replacement $G$ is the unit-disk graph induced
by a finite set $S$ of lattice sites. Four occupied sites $P = (p_1,p_2,p_3,p_4)$
are the pins, and four lattice directions $bold(r) = (r_1,r_2,r_3,r_4)$ are the
outgoing wires. Write $s_i = p_i + r_i$ for the first site *outside* the gadget
along wire $i$. The pair $(P, bold(r))$ is a *crossing frame* when the four
wires can leave the gadget without interfering:

- *G1.* each $s_i$ is a strict vertex of the convex hull of $S union {s_1,s_2,s_3,s_4}$;
- *G2.* walking around that hull, the two logical channels alternate:
  odd labels ${s_1,s_3}$ and even labels ${s_2,s_4}$ interleave;
- *G3.* every $r_i$ points away from the centroid of the four $s_i$;
- *G4.* the four exterior rays miss every occupied site except their own pin,
  and they stay pairwise out of blockade range.

#figure(
  grid(
    columns: (1fr, 1fr),
    align: (center + horizon),
    column-gutter: 12pt,
    cetz-canvas({
      import cetz.draw: *
      let p1 = (0, 1.45)
      let p2 = (-1.45, 0)
      let p3 = (0, -1.45)
      let p4 = (1.45, 0)
      line(p1, p3, stroke: 1.15pt + ink)
      line(p2, p4, stroke: 1.15pt + ink)
      node-at(p1, [1], fill: pin-fill, paint: pin-stroke, r: 0.26)
      node-at(p2, [2], fill: pin-fill, paint: pin-stroke, r: 0.26)
      node-at(p3, [3], fill: pin-fill, paint: pin-stroke, r: 0.26)
      node-at(p4, [4], fill: pin-fill, paint: pin-stroke, r: 0.26)
    }),
    cetz-canvas({
      import cetz.draw: *
      let occupied = ((0, 0), (-1, 0), (0, 1), (1, 0), (0, -1))
      let pin-sites = ((-1, 0), (0, 1), (1, 0), (0, -1))
      let labs = ([1], [2], [3], [4])
      let interfaces = ((-2, 0), (0, 2), (2, 0), (0, -2))
      let sp = 0.52
      let pt(xy) = (xy.at(0) * sp, xy.at(1) * sp)
      let grid-paint = rgb("#8f9ba6")
      for x in range(-3, 4) {
        for y in range(-3, 4) {
          if x < 3 {
            line(pt((x, y)), pt((x + 1, y)), stroke: 0.6pt + grid-paint)
          }
          if y < 3 {
            line(pt((x, y)), pt((x, y + 1)), stroke: 0.6pt + grid-paint)
          }
        }
      }
      for x in range(-3, 4) {
        for y in range(-3, 4) {
          let p = (x, y)
          if occupied.contains(p) or interfaces.contains(p) { continue }
          circle(pt(p), radius: 0.075, fill: white, stroke: 0.45pt + rgb("#b7c0c8"))
        }
      }
      line(pt((-1, 0)), pt((1, 0)), stroke: 0.95pt + ink)
      line(pt((0, -1)), pt((0, 1)), stroke: 0.95pt + ink)
      for i in range(4) {
        line(
          pt(pin-sites.at(i)),
          pt((interfaces.at(i).at(0) * 1.35, interfaces.at(i).at(1) * 1.35)),
          stroke: (paint: pin-stroke, thickness: 0.85pt, dash: "dashed"),
        )
        circle(pt(interfaces.at(i)), radius: 0.08, fill: pin-stroke, stroke: none)
      }
      node-at(pt((0, 0)), [], r: 0.16)
      for i in range(4) {
        node-at(pt(pin-sites.at(i)), labs.at(i), fill: pin-fill, paint: pin-stroke, r: 0.18)
      }
      content(pt((-2.85, 0)), text(size: 7.4pt, fill: muted)[$s_1$])
      content(pt((0, 2.85)), text(size: 7.4pt, fill: muted)[$s_2$])
      content(pt((2.85, 0)), text(size: 7.4pt, fill: muted)[$s_3$])
      content(pt((0, -2.85)), text(size: 7.4pt, fill: muted)[$s_4$])
    }),
  ),
  caption: [(a)~CROSS, with edges $1$--$3$ and $2$--$4$; all four vertices are pins.
    (b)~A five-site plus on the square lattice, forming a valid crossing frame.
    The interfaces $s_i$ sit on the convex hull and alternate odd/even.],
) <fig-problem>

@fig-problem~(b) satisfies G1--G4. Adding a site on one
exterior ray blocks that corridor and violates G4.

The search problem is now concrete. Given a lattice (triangular TLSG, or square
KSG), a window $W$, an atom count $N$, and an offset $c$, find occupied sites
$S subset.eq W$, pins $P subset.eq S$, and rays $bold(r)$ such that

$
  abs(S) = N,
  quad S " is connected",
  quad (P,bold(r)) " satisfies G1"--"G4",
  quad tilde(alpha)(G[S]) = T + c.
$

= Two stages

The algorithm first finds a certified seed, then shrinks it.

#figure(
  {
    set text(font: "Libertinus Sans", size: 8.6pt)
    let boxy(title, body) = block(
      width: 100%,
      inset: (x: 8pt, y: 7pt),
      fill: white,
      stroke: 0.7pt + pin-stroke,
      radius: 3pt,
      {
        text(weight: "semibold", fill: accent, title)
        v(2pt)
        set text(font: "Libertinus Serif", size: 8.3pt)
        body
      },
    )
    grid(
      columns: (1fr, 18pt, 1fr, 18pt),
      align: center + horizon,
      boxy[Stage I --- search][Produce any certified seed: occupied sites, pins, and rays whose reduced tensor is $T+c$.],
      text(size: 14pt, fill: pin-stroke)[$arrow.r$],
      boxy[Stage II --- optimize][Apply certified local rewrites that preserve the tensor up to a new constant, seeking a smaller $N$.],
    )
  },
  caption: [The complete procedure.],
)

Stage I has two SAT engines that share the occupancy encoding of @occ-sat:

- *Joint SAT* chooses $S$, $P$, and $bold(r)$ in one formula. This is the
  blank-window discovery method.
- *Frame-first SAT* enumerates geometrically valid frames, then occupancy SAT
  chooses the remaining occupied sites. This is the systematic window traversal,
  and also the engine used later to re-synthesize an interior.

Stage II accepts any verifier-valid seed. When the rewrite neighbourhood and
SAT budget yield no smaller certified descendant, the seed itself is the answer.

= Occupancy SAT for a fixed frame <occ-sat>

Fix a finite lattice window and a geometrically valid frame: four pin sites
and their outgoing rays. The *candidate lattice graph* $H$ contains every
site in the window allowed by this frame. Two candidate sites are joined by
an edge exactly when their distance is within the blockade range.

SAT chooses a subset $S subset.eq V(H)$ of these sites. The resulting gadget
is the induced graph $H[S]$: it contains all selected sites and every edge
of $H$ between them. Edges are determined by geometry, so SAT chooses vertices,
not edges. Fix an atom count $N$ and an offset $c$. We seek a connected
$H[S]$ with $abs(S)=N$ and reduced tensor $T+c$.

== Variables and clauses

The candidate graph $H$ and its edges are fixed. SAT chooses which vertices remain
in the gadget. A clause is a disjunction of literals, such as
$(not y_u or not y_v)$; the complete CNF is the conjunction of all clauses.
A unit clause $(x)$ forces $x=1$, and $(not x)$ forces $x=0$.

The encoding uses two distinct kinds of choices:

- $x_v=1$ means that candidate site $v$ belongs to the gadget, $v in S$.
  These variables are shared by all boundary states.
- $y_v^sigma=1$ means that $v$ belongs to one independent set $I_sigma$
  used to certify the lower bound for state $sigma$. Each retained state
  gets its own fresh copy of these variables.

Adjacent sites may both belong to $S$: their edge is part of the gadget.
Thus the edge clauses below constrain $y^sigma$, not $x$.
For a candidate path `a-b-c`, the independent-set clauses are
$(not y_a^sigma or not y_b^sigma)$ and
$(not y_b^sigma or not y_c^sigma)$, while $x_a=x_b=x_c=1$ is allowed.

== Selecting atoms

Add four unit clauses $(x_(p_i))$, one per pin,
and impose

$
  sum_(v in V(H)) x_v = N.
$

The pins therefore belong to the selected graph, and the graph has exactly
$N$ vertices. The sum is shorthand for the counter clauses in @counter-cnf;
it is not a single CNF clause.

== Lower bounds: one independent-set witness per state

Define the *size* of a largest independent set with boundary $sigma$ by

$
  A_S (sigma) = max { abs(I) : I subset.eq S, I " independent in " H,
    I inter P = {p_i : sigma_i = 1} }.
$

An infeasible state has value $-infinity$. For a finite target entry,
we need $A_S (sigma) >= k_sigma$, where $k_sigma = T_sigma+c$.
To prove this lower bound, SAT only needs to exhibit one suitable independent
set. It does not need to prove that this witness is maximum.

For each retained state $sigma$, add the following four groups of constraints:

1. *Use only selected sites.* For every $v in V(H)$, append
   $
     (not y_v^sigma or x_v).
   $
   If the witness takes $v$, the gadget must contain $v$.
   A selected site need not occur in this particular witness.

2. *Fix the boundary.* For each pin $p_i$, append the unit clause
   $(y_(p_i)^sigma)$ if $sigma_i=1$, or $(not y_(p_i)^sigma)$ if $sigma_i=0$.
   A zero bit excludes a pin from this independent set; the pin still
   belongs to the gadget because $x_(p_i)=1$.

3. *Forbid adjacent takes.* For every candidate edge ${u,v} in E(H)$, append
   $
     (not y_u^sigma or not y_v^sigma).
   $
   These clauses make $I_sigma = {v : y_v^sigma=1}$ independent.

4. *Require the witness size.* Encode
   $
     sum_(v in V(H)) y_v^sigma = k_sigma.
   $
   This emits the counter clauses in @counter-cnf with a fresh set of
   auxiliary variables for this witness.

Together these groups imply $A_S (sigma) >= k_sigma$. Equality of the
*witness size* does not assert equality of the *maximum*. If a larger set
exists, non-pin vertices can be removed until its size is $k_sigma$, provided
$k_sigma >= abs(sigma)$, where $abs(sigma)$ counts the occupied pins.
If $k_sigma < abs(sigma)$ or $k_sigma > abs(V(H))$, the witness constraints
are unsatisfiable. Upper bounds require the separate encoding below.

=== Which lower states can be omitted?

The subset relation compares the sets of pins whose bits are $1$.
Given a witness for a strict superset $sigma subset tau$, remove its
$abs(tau)-abs(sigma)$ extra pins. The remaining set has boundary $sigma$,
so it proves

$
  A_S (sigma) >= T_tau + c - (abs(tau)-abs(sigma)).
$

Consequently a finite state $sigma$ can be omitted exactly when some finite strict superset $tau$ satisfies

$
  T_tau - (abs(tau)-abs(sigma)) >= T_sigma.
$

A superset alone is not sufficient: its target value must satisfy this
inequality. If a superset is itself omitted, the implication continues to a
larger superset and eventually reaches a retained state.

#note[Example: the four CROSS witnesses are in one SAT instance][
  Write bits in pin order $(p_1,p_2,p_3,p_4)$. The retained states are
  `1100`, `0110`, `0011`, and `1001`, each with $T_sigma=2$.
  At offset $c=7$, the `1100` copy has pin unit clauses
  $
    (y_(p_1)^(1100)), quad (y_(p_2)^(1100)), quad
    (not y_(p_3)^(1100)), quad (not y_(p_4)^(1100)),
  $
  plus the site and edge clauses above and a counter requiring exactly
  $9$ witness vertices: these two pins and seven internal vertices.
  The other three copies impose their own pin patterns and size-$9$ counts.
  All four copies share the same $x_v$, so one selected graph must support
  all four witnesses. They may reuse vertices across copies; the four sets
  are not required to be disjoint or simultaneously independent as a union.

  Removing $p_2$ from the `1100` witness gives an $8$-vertex witness for
  `1000`; removing both pins gives a $7$-vertex witness for `0000`.
  These are the required lower bounds $1+c$ and $c$, so separate copies for
  these states would be redundant.
]

=== How an exact count becomes CNF <counter-cnf>

Start with three candidate points: *choose exactly two of points 1, 2, and 3*.
Let $y_1,y_2,y_3$ say whether each point is chosen. The allowed answers are
`110`, `101`, and `011`.

Build two sequential counters: the first permits *at most two chosen points*;
the second permits *at most one unchosen point*. Together they require exactly
two chosen points. Each counter passes count thresholds along the ordered
points using auxiliary Boolean flags.

*First counter: at most two chosen points.*
Introduce three flags:

#table(
  columns: (auto, 1fr),
  table.header([*Flag*], [*When it must become true*]),
  [$q_(1,1)$], [After point 1: at least one point has been chosen.],
  [$q_(2,1)$], [After points 1 and 2: at least one point has been chosen.],
  [$q_(2,2)$], [After points 1 and 2: two points have been chosen.],
)

The at-most-two counter has these five clauses:

#table(
  columns: (1fr, 1.3fr),
  table.header([*Clause*], [*What it enforces*]),
  [$(not y_1 or q_(1,1))$], [Choosing point 1 sets the first flag.],
  [$(not y_2 or q_(2,1))$], [Choosing point 2 means the first two contain at least one chosen point.],
  [$(not q_(1,1) or q_(2,1))$], [A count already reached stays reached.],
  [$(not y_2 or not q_(1,1) or q_(2,2))$], [If point 1 was chosen and point 2 is chosen, the count reaches two.],
  [$(not y_3 or not q_(2,2))$], [If the first two were chosen, point 3 cannot be chosen.],
)

Try `111`: choosing point 1 forces $q_(1,1)=1$; choosing point 2 then
forces $q_(2,2)=1$. With point 3 also chosen, both literals in the last
clause are false. No assignment of the flags can rescue this choice.
For `110`, all three flags can be true and all five clauses hold.

This counter alone also permits `100` and `000`: it only imposes an upper
bound.
*Second counter: at most one unchosen point.*
Count $not y_1, not y_2, not y_3$ using two new flags: $r_1$ must be true
if point 1 is unchosen, and $r_2$ must be true if either of the first two
points is unchosen. Add these five clauses:

#table(
  columns: (1fr, 1.3fr),
  table.header([*Clause*], [*What it enforces*]),
  [$(y_1 or r_1)$], [Leaving point 1 unchosen sets the first flag.],
  [$(y_2 or r_2)$], [Leaving point 2 unchosen sets the second flag.],
  [$(not r_1 or r_2)$], [An unchosen point already counted stays counted.],
  [$(y_2 or not r_1)$], [If point 1 was unchosen, point 2 must be chosen.],
  [$(y_3 or not r_2)$], [If either of the first two was unchosen, point 3 must be chosen.],
)

Try `100`: point 2 is unchosen, so the second clause forces $r_2=1$.
Point 3 is also unchosen, making the last clause false. Thus choosing only
point 1 fails this counter. For `110`, set $r_1=r_2=0$: all five clauses
hold. Combined with the first counter, this accepts exactly two chosen points.

These ten clauses are the complete sequential-counter encoding for this
example. The flags are SAT variables constrained by the clauses; no separate
procedure enumerates choices or updates counts during formula construction.
The assignments above only illustrate how the clauses enforce the count.

For CROSS, the same construction limits each witness to *at most nine chosen
points* and *at most the number of candidate points minus nine unchosen points*.

*General counter clauses.* The formulas below describe the same construction
for any number of points and any required count.

Encode an exact count of $k$ among $m$ literals $(ell_1,...,ell_m)$
by imposing two at-most constraints, each with its own counter:

$
  sum_(i=1)^m ell_i <= k,
  quad sum_(i=1)^m (not ell_i) <= m-k.
$

Here a literal is counted as $1$ when true. The first constraint limits the number
of true literals; the second limits the number of false literals.
The same construction encodes both the atom count and each witness count.

For an at-most bound $b$ with $0<b<m$, introduce counter variables
$q_(i,j)$ for $1<=i<=m-1$ and $1<=j<=b$. A true count of at least $j$
among the first $i$ literals forces $q_(i,j)$ true. Add these clause families (an empty index range emits nothing):

$
  (not ell_i or q_(i,1))
  quad (1<=i<=m-1),
$
$
  (not q_(i-1,1) or q_(i,1))
  quad (2<=i<=m-1),
$
$
  (not ell_i or not q_(i-1,j-1) or q_(i,j))
  quad (2<=j<=b, j<=i<=m-1),
$
$
  (not q_(i-1,j) or q_(i,j))
  quad (2<=j<=b, j+1<=i<=m-1),
$
$
  (not ell_i or not q_(i-1,b))
  quad (b+1<=i<=m).
$

The first four families propagate count thresholds. The last forbids one
more true literal after the count has reached $b$. Only forward implications
are needed: any actual overflow forces a contradiction, and any assignment
within the bound can be extended by its actual prefix counts.
For $b=0$, add $(not ell_i)$ for every literal. For $b>=m$, no clauses
are needed. For $b<0$, add the empty clause, making the CNF UNSAT.

== Upper bounds: a frontier BDD

The lower witnesses establish existence. We must also rule out independent
sets that are too large for their boundary state. Define the monotone completion

$
  C_T(sigma) = max{ T_tau : tau subset.eq sigma, T_tau != -infinity }.
$

Then every boundary state is capped by

$
  A_S (sigma) <= C_T(sigma) + c.
$

The SAT encoding forbids every independent set larger than the cap. It does so
with a frontier binary decision diagram over a fixed vertex order. At vertex
$v$ the path either *skips* $v$ or *takes* $v$; a take edge exists when every
already-taken neighbour is compatible with $v$. After each prefix, states with the
same active frontier are merged. The count is truncated at
$C_T(sigma)+c+1$, and every terminal that reaches this overflow count is excluded. Among
six lattice-axis orders, the smallest BDD is kept.

For a retained upper state $sigma$, let $B=C_T(sigma)+c$.
Introduce a variable $r_(t,F,h)$ for
reachability after $t$ processed vertices, with active taken frontier $F$
and count $h$ truncated at $B+1$. Its clauses are:

- *Start:* the unit clause $(r_(0,emptyset,0))$.
- *Skip edge $a arrow.r b$:* $(not r_a or r_b)$.
- *Take edge $a arrow.r b$ at vertex $v$:* $(not r_a or not x_v or r_b)$.
- *Overflow terminal $b$ with count $B+1$:* the unit clause $(not r_b)$.

Here $a$ and $b$ abbreviate the full node indices. At a pin, the construction
omits edges whose take/skip decision disagrees with $sigma$; it also omits
all takes blocked by the active frontier. Skip edges need no $x_v$ guard:
an independent set may skip even a selected vertex.

These clauses force propagation along *every* available path, rather than
asking SAT to choose one path. If any independent set of the selected graph
exceeds $B$, its path forces an overflow terminal true, contradicting that
terminal's negative unit clause. If no such set exists, assigning precisely
the reachable nodes true satisfies this encoding.

An upper-bound state $sigma$ can be omitted if a strict subset $tau$
already implies its cap:
$C_T(tau)+abs(sigma)-abs(tau) <= C_T(sigma)$.
Removing the extra pins from any independent set proves this implication.

#figure(
  cetz-canvas({
    import cetz.draw: *
    // path a-b-c
    node-at((-4.6, 0.9), [a], r: 0.22)
    node-at((-3.4, 0.9), [b], r: 0.22)
    node-at((-2.2, 0.9), [c], r: 0.22)
    line((-4.38, 0.9), (-3.62, 0.9), stroke: 0.9pt + ink)
    line((-3.18, 0.9), (-2.42, 0.9), stroke: 0.9pt + ink)
    content((-3.4, 1.45), text(size: 7.6pt, fill: muted)[candidate path, cap $= 1$])

    let skip = rgb("#6a7d8a")
    let take = rgb("#b4532a")
    let n(pos, lab) = {
      circle(pos, radius: 0.28, fill: white, stroke: 0.8pt + ink)
      content(pos, text(size: 6.8pt, lab))
    }
    n((0.2, 2.0), [0])
    n((-1.3, 0.7), [0])
    n((1.7, 0.7), [1])
    n((-2.1, -0.7), [0])
    n((-0.5, -0.7), [1])
    n((0.9, -0.7), [1])
    n((-2.1, -2.0), [0])
    n((-0.5, -2.0), [1])
    n((0.9, -2.0), [1])
    n((2.6, -2.0), [$bot$])

    let arr(a, b, col) = {
      line(a, b, stroke: 0.75pt + col, mark: (end: "stealth", fill: col, scale: 0.45))
    }
    arr((0.05, 1.74), (-1.15, 0.96), skip)
    arr((0.38, 1.74), (1.55, 0.96), take)
    arr((-1.5, 0.46), (-2.05, -0.42), skip)
    arr((-1.1, 0.46), (-0.6, -0.42), take)
    arr((1.5, 0.46), (1.05, -0.42), skip)
    arr((-2.1, -0.98), (-2.1, -1.72), skip)
    arr((-1.9, -0.96), (-0.7, -1.74), take)
    arr((-0.5, -0.98), (-0.5, -1.72), skip)
    arr((0.9, -0.98), (0.9, -1.72), skip)
    arr((1.1, -0.96), (2.4, -1.74), take)

    content((-0.8, 1.55), text(size: 6.6pt, fill: skip)[skip $a$])
    content((1.35, 1.55), text(size: 6.6pt, fill: take)[take $a$])
    content((0.2, -2.55), text(size: 7.2pt, fill: muted)[labels are running IS size; $bot$ is the forbidden cap $2$])
  }),
  caption: [Frontier BDD for the path $a$--$b$--$c$ with upper bound $1$. Taking $a$ then $c$ would reach size $2$ and is cut off. Adjacent takes ($a$ then $b$, or $b$ then $c$) are omitted: those vertices blockade each other.],
) <fig-bdd>

== Connectivity

The selected sites must induce a connected graph. Let $z_v^t$ mean that
selected vertex $v$ is reachable from pin $p_1$ by a path of at most $t$ edges using only selected vertices. The encoding is direct layered reachability:

$
  z_v^0 arrow.l.r.double (v = p_1),
$
$
  z_v^(t+1) arrow.l.r.double
  (z_v^t or (x_v and or.big_(u in N_H(v)) z_u^t)),
$
$
  x_v => z_v^(N-1).
$

To express reachability in CNF, add the root unit clause
$(z_(p_1)^0)$ and $(not z_v^0)$ for every other vertex. For each
$t=0,...,N-2$ and vertex $v$, add

$
  (not z_v^t or z_v^(t+1)),
  quad (not z_v^(t+1) or x_v),
$
$
  (not z_v^(t+1) or z_v^t or or.big_(u in N_H(v)) z_u^t),
$

and, for every neighbour $u in N_H(v)$,

$
  (not x_v or not z_u^t or z_v^(t+1)).
$

Finally add $(not x_v or z_v^(N-1))$ for every vertex.
The forward clauses propagate reachability; the reverse clauses prevent SAT
from declaring an isolated selected vertex reachable without a path.

The depth $N-1$ is exact: in a connected $N$-vertex graph every vertex has a
simple path of length at most $N-1$ from $p_1$, and a selected vertex in
another component stays unreachable.

#figure(
  cetz-canvas({
    import cetz.draw: *
    let p = (0, 0)
    let u = (1.5, 0)
    let v = (3.0, 0)
    let w = (5.0, 0)
    line(p, u, stroke: 1pt + ink)
    line(u, v, stroke: 1pt + ink)
    node-at(p, [$p_1$], fill: pin-fill, paint: pin-stroke, r: 0.28)
    node-at(u, [$u$], fill: is-fill, r: 0.26)
    node-at(v, [$v$], fill: is-fill, r: 0.26)
    node-at(w, [$w$], fill: rgb("#e6e6e6"), paint: rgb("#8a8a8a"), r: 0.26)
    content((0, -0.7), text(size: 7.2pt)[$z^0$])
    content((1.5, -0.7), text(size: 7.2pt)[$z^1$])
    content((3.0, -0.7), text(size: 7.2pt)[$z^2$])
    content((5.0, -0.7), text(size: 7.2pt, fill: muted)[unreachable])
    content((1.5, 0.75), text(size: 7.6pt, fill: muted)[selected path])
    content((5.0, 0.75), text(size: 7.6pt, fill: muted)[isolated])
  }),
  caption: [Layered reachability. The selected path $p_1$--$u$--$v$ becomes reachable layer by layer. An isolated candidate vertex $w$ stays unselected: it remains unreachable.],
)

Putting the pieces together, the occupancy formula is

$
  Phi = Phi_"select" and Phi_"connect"
  and and.big_(sigma in L) Phi_"witness"^sigma
  and and.big_(sigma in U) Phi_"BDD"^sigma,
$

where $L$ and $U$ are the essential lower and upper states.

= Choosing the frame

The occupancy formula assumes the pins and rays are already known. Stage I
offers two ways to obtain them.

== Joint occupancy-and-frame SAT

One CNF chooses everything at once. Besides the occupancy variables $x_v$,
each geometrically eligible triple (boundary label, site, lattice direction)
gets a frame variable $f_(i,v,d)$. The extra clauses say:

1. exactly $N$ sites are occupied;
2. each label $i in {1,2,3,4}$ chooses exactly one pin-and-ray;
3. a chosen pin is occupied, and the four pins are distinct;
4. local corridor, outward-direction, and alternating-interface constraints
   (the SAT-side counterpart of G1--G4);
5. the tensor and connectivity constraints of @occ-sat, now relative to the
   chosen pins.

Only eligible port choices receive variables. One reference port is fixed by a
canonical translation and a first-ray direction. That choice specifies the
finite instance: a known abstract gadget appears in this CNF when a
symmetry-equivalent copy fits the anchored window.

A SAT solver searches for a satisfying assignment. A model is materialised as a concrete lattice
patch and checked geometrically *before* the reduced tensor is recomputed. A
rejected model is blocked and the solver is restarted. Timed-out runs are
recorded as unknown work.

== Frame-first SAT

The other engine moves geometry outside the CNF. Each locally valid pair
$(p,r)$ is a vertex of a compatibility graph; four mutually compatible ports
form a clique; surviving ordered frames are checked against the full G1--G4
conditions and canonicalised under lattice symmetries. For each remaining
frame the allowed candidate set

$
  U(P,bold(r)) = { v in W : "placing" v "preserves G1, G3, G4"}
$

is computed, and the occupancy SAT of @occ-sat runs on $H = G[U]$.

#note[Example: occupancy SAT is cheap once the CROSS frame is known][
  On the known 23-site triangular CROSS frame, occupancy SAT returns a
  verifier-valid gadget at offset $7$ immediately. A blank frame-first search
  spends most of its time enumerating frames that are easy UNSAT. Joint SAT
  avoids that explicit loop, at the cost of one larger, seed-sensitive CNF.
  The two engines are complementary.
]

#figure(
  {
    let sites = (
      (0, 3), (6, 6), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 3),
      (3, 4), (3, 5), (4, 1), (4, 2), (4, 4), (4, 5), (5, 1), (5, 2),
      (5, 3), (5, 4), (5, 5), (6, 2), (6, 3), (6, 5), (7, 4),
    )
    let pins = ((5, 1), (0, 3), (3, 5), (6, 6))
    let pin-labs = ([1], [2], [3], [4])
    let phys(pt) = {
      let x = pt.at(0) + (if calc.odd(pt.at(1)) { 0.5 } else { 0.0 })
      let y = pt.at(1) * calc.sqrt(3) / 2
      (x, y)
    }
    let scale = 0.42
    let origin = (1.55, 0.35)
    let P(pt) = {
      let q = phys(pt)
      (origin.at(0) + scale * q.at(0), origin.at(1) + scale * q.at(1))
    }
    cetz-canvas({
      import cetz.draw: *
      for a in sites {
        for b in sites {
          if a.at(0) < b.at(0) or (a.at(0) == b.at(0) and a.at(1) < b.at(1)) {
            let pa = phys(a)
            let pb = phys(b)
            let dx = pa.at(0) - pb.at(0)
            let dy = pa.at(1) - pb.at(1)
            if dx * dx + dy * dy < 1.21 {
              line(P(a), P(b), stroke: 0.7pt + rgb("#6e7c86"))
            }
          }
        }
      }
      for pt in sites {
        let is-pin = pins.contains(pt)
        circle(
          P(pt),
          radius: 0.13,
          fill: if is-pin { pin-fill } else { site-fill },
          stroke: 0.75pt + if is-pin { pin-stroke } else { site-stroke },
        )
      }
      for i in range(4) {
        let pos = P(pins.at(i))
        let dy = if i == 0 { -0.28 } else { 0.28 }
        content(
          (pos.at(0), pos.at(1) + dy),
          text(size: 7.4pt, weight: "semibold", fill: pin-stroke, pin-labs.at(i)),
        )
      }
    })
  },
  caption: [A certified 23-site triangular CROSS gadget (offset $7$). Pins are labelled; remaining disks are internal occupied sites.],
) <fig-cross23>

= Independent certification

Every satisfying assignment is materialised as the exact unit-disk graph. An
independent verification step then recomputes $tilde(alpha)(G)$, checks the common
offset, and re-runs G1--G4 and connectivity. The same boundary is used during
discovery and after every rewrite. Auxiliary SAT encodings may change; the
semantic test stays fixed.

= Shrinking a certified gadget

Stage II starts from a verifier-accepted gadget and looks for a smaller one
with the same reduced tensor up to a (possibly different) constant. A local
change is accepted when the *interface signature* of the edited region is
preserved.

== Interface signatures

Let $Q$ be a connected patch inside a certified gadget, and let
$C = (c_1,...,c_k)$ be the vertices through which $Q$ meets the unchanged
exterior. For an interface state $tau in {0,1}^k$,

$
  A_Q^C (tau)
  = max{ abs(I) : I " independent in " Q,
         I inter C = {c_i : tau_i = 1} }.
$

The vector of all such values is everything the exterior needs to know about
$Q$ when maximising an independent set. Replacing $Q$ by $Q'$ is valid in
every compatible context when the two signatures differ by one constant $delta$
and the infeasible-state pattern is unchanged. Gluing either patch to the same
exterior then shifts every global conditioned optimum by $delta$, so the
reduced tensor of the whole gadget shifts by an allowed offset.

On a lattice the replacement must also keep attachment coordinates or port
roles, leave the exterior blockade pattern unchanged, remain connected, and
satisfy G1--G4. The whole-gadget verifier is the final certificate.

== Two operator families

The optimizer uses two families consistent with that contract.

*Exact local replacement.* Two closed-form rules apply whenever their
structural premises hold, then the whole gadget is re-certified.

- *Even boundary-tail contraction.* If a pin is a leaf, and the unique path
  $p$--$m$--$e$ has $m$ and $e$ off the boundary, delete $p$ and $m$ and
  promote $e$ to the pin. The outgoing ray is updated so that it points
  through the deleted tail.
- *Opposite leaf-pin contraction.* If a pair of opposite pins $(p_1,p_3)$ or
  $(p_2,p_4)$ are both leaves, delete both leaves and promote their unique
  neighbours. This shortens both wires of one channel together, which is why
  the pins are treated as a pair.

#figure(
  cetz-canvas({
    import cetz.draw: *
    // before: tail p - m - e
    let e = (0, 0)
    let m = (1.15, 0)
    let p = (2.3, 0)
    line((-0.9, 0.55), e, stroke: 0.85pt + rgb("#9aa7b0"))
    line((-0.9, -0.55), e, stroke: 0.85pt + rgb("#9aa7b0"))
    line(e, m, stroke: 1pt + ink)
    line(m, p, stroke: 1pt + ink)
    line(p, (3.15, 0), stroke: (paint: pin-stroke, thickness: 0.85pt, dash: "dashed"))
    node-at((-0.9, 0.55), [], r: 0.16)
    node-at((-0.9, -0.55), [], r: 0.16)
    node-at(e, [e], r: 0.22)
    node-at(m, [m], r: 0.22)
    node-at(p, [p], fill: pin-fill, paint: pin-stroke, r: 0.22)
    content((1.15, -0.85), text(size: 8pt)[before: leaf tail])

    line((3.7, 0), (4.35, 0), stroke: 1.1pt + pin-stroke, mark: (end: "stealth", fill: pin-stroke, scale: 0.6))

    let e2 = (5.4, 0)
    line((4.5, 0.55), e2, stroke: 0.85pt + rgb("#9aa7b0"))
    line((4.5, -0.55), e2, stroke: 0.85pt + rgb("#9aa7b0"))
    line(e2, (6.35, 0), stroke: (paint: pin-stroke, thickness: 0.85pt, dash: "dashed"))
    node-at((4.5, 0.55), [], r: 0.16)
    node-at((4.5, -0.55), [], r: 0.16)
    node-at(e2, [e], fill: pin-fill, paint: pin-stroke, r: 0.22)
    content((5.6, -0.85), text(size: 8pt)[after: $e$ is the new pin])
  }),
  caption: [Even boundary-tail contraction. The two-edge tail is removed and the former internal endpoint becomes the pin. Opposite leaf-pin contraction does the same to both wires of one channel at once.],
) <fig-tail>

*Interface-constrained re-synthesis.* Hold a separator, the external port
roles, and the required signature fixed, then ask occupancy SAT for a smaller
interior. Moving one pin by a single lattice step proposes a new frame. The
semantic rule is to preserve the interface contract.

The optimizer first closes the graph of all certified direct descendants, then
spends SAT budget on re-synthesis from those states. Exact coordinate/pin/ray
keys merge identical converging paths. A smaller immediate child can be a dead
end; a different equivalent child can expose a later rewrite.

== Walking CROSS from 28 sites down to 23

The independently discovered 28-site triangular CROSS gadget is a complete
worked example of both operator families. Starting from that seed, the
optimizer produces three certified steps:

#figure(
  table(
    columns: (auto, auto, auto, 1fr),
    align: (center, center, left, left),
    fill: (x, y) => if y == 0 { rgb("#edf3f8") } else { white },
    table.header([*step*], [*$N$*], [*rule*], [*what changed*]),
    [seed], [$28 arrow.r 28$], [(Stage I)], [joint SAT, offset $10$, independently verified],
    [1], [$28 arrow.r 26$], [tail contraction], [an even boundary tail is removed; two sites leave],
    [2], [$26 arrow.r 25$], [re-synthesis], [one pin moves by a lattice step; occupancy SAT rebuilds a smaller interior],
    [3], [$25 arrow.r 23$], [opposite leaves], [a pair of opposite leaf pins is contracted],
  ),
  caption: [Certified rewrite trace for CROSS. The 23-site gadget in @fig-cross23 is the end of this path, at offset $7$.],
)

Each step is accepted after the unchanged verifier agrees on the new tensor
and geometry. The same Stage II applied to a 9-site CROSS+EDGE seed exhausts
the direct-contraction closure and a budget of 256 fixed-frame SAT calls; the
9-site gadget is the best result under the current rules and budget.

= Results

#figure(
  table(
    columns: (1.2fr, 1.15fr, 1.25fr, 0.9fr, 1.7fr),
    align: (left, center, center, center, left),
    fill: (x, y) => if y == 0 { rgb("#edf3f8") } else { white },
    table.header(
      [*target*], [*Stage I*], [*Stage II*], [*offset*], [*status*],
    ),
    [CROSS], [28 sites], [23 sites], [$10 arrow.r 7$],
      [reduced and verified],
    [CROSS+EDGE], [9 sites], [9 sites], [$1 arrow.r 1$],
      [unchanged after budget; verified],
  ),
  caption: [End-to-end outcome of search followed by optimization.],
)

#note[Scope of these numbers][
  Completeness holds inside the scheduled finite space: one window, one
  atom-count range, one canonical anchor, and a bounded rewrite neighbourhood.
  A SAT timeout is unfinished work. A Stage II fixed point is a local
  fixed point under the chosen rules and budget.
]

For a new four-pin target, first compute its reduced tensor and choose a
lattice window, atom count, and offset. Build the SAT formula either by
fixing a frame first or by including frame choices in the formula. Extract
the selected graph from a satisfying assignment and independently verify
its tensor, connectivity, and geometry. Then search for smaller replacements,
accepting each rewrite only after the same checks pass.
