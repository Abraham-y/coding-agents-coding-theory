# Failed ideas log (CLAUDE.md rule 7)

One line on what was tried and why it died. Nothing here gets retried without a
reason that addresses the stated cause of death.

---

**Extend VT_0(n) by one word to beat 316 at (s=1,n=12).**
Dead. `tests/test_verify.py` maximality probe: VT_0(n) is a MAXIMAL independent
set for every n = 6…20 — no word can be added to any of them. Any improvement at
(1,12) must be a structurally different code, not an extension of VT.

**Plain ILS (ARW-style 1,2-swaps) as the workhorse for s=1 cells.**
Dead as a primary method. Reaches 138 on (s=1,n=11) where the optimum is 172.
Cause: that optimum is *algebraic* — a global congruence condition (VT) — and
vertex-swapping local search never stumbles onto one. Local search is only
competitive on the s≥2 cells where the incumbents also came from MIS solvers.

**SAT on the full graph: clique-cover AtMost1 + AtLeast-k totalizer.**
Dead. On the *known-feasible* instance (s=2,n=11,k=24) the encoding is 2.2M
clauses and Cadical153 returned UNKNOWN after 462s. Cause: CDCL is poor at dense
cardinality-driven MIS; the totalizer over thousands of variables dominates.

**SAT on the complement-orbit graph (halves vars, halves the cardinality target).**
Dead too. (s=2,n=11), k=12 orbits (= size 24, known feasible): 724k clauses,
UNKNOWN after 389s. Halving the instance did not change the character of the
difficulty. Conclusion: SAT is the wrong rung for this problem, which is exactly
what the rewritten ladder in CLAUDE.md predicted. Do not escalate to SAT again
for MIS-shaped deletion-code cells; the reason is structural, not budgetary.

**Consequence — strategy change, recorded here because it is the important one.**
A better *lower* bound (a bigger code) is a construction and needs no
certificate; a better *upper* bound obtained from SAT would need a DRAT proof and
is out of scope (rule 11). But an upper bound obtained from LP DUALITY needs only
an exact rational dual vector, which is cheap to verify with `fractions.Fraction`
and fully satisfies rules 10 and 12. So the upper-bound side is the *more*
tractable side for this sprint, not the less. Pivoting effort there.

**Reversal-invariant codes at (s=2, n=12), motivated by a parity argument.**
Dead, and the motivating argument was sound but the conclusion was wrong. The
reasoning: complement orbits all have size 2, so complement-invariant codes have
even size, yet the record at (2,12) is **37, odd** — so the record code is
provably not complement-invariant. The reversal group has 64 weight-1 orbits
(palindromes) and can therefore produce odd sizes, making it the natural next
group to try. Result: **36**, no better than complement, though it did reach odd
sizes (35 on four of five seeds), confirming parity was never the obstacle.
Cause of death: reversal is a far harsher constraint than complement —
**924 of 2080 reversal orbits are internally conflicting** (a string and its
reverse usually share a length-10 common subsequence) versus only 37 of 2048 for
complement. Only 1156 usable orbits remain. Do not retry other small groups here
without first checking the internal-conflict count; it is a one-line diagnostic
and it predicts the outcome.

**First k=4 attempt for alpha_4: killed by memory exhaustion at ~7.6 GB.**
Three distinct inefficiencies, all fixed, recorded because each is easy to
reintroduce:
1. `leaf_stats` cached a 2 x 3^k integer array for every TOP-LEVEL DP state. Each
   such state is visited exactly once, so this was pure waste — ~32 GB at k=4 on
   its own. Fixed with `cache=False` at the top level only; sub-states, which are
   heavily shared, are still cached.
2. States were canonicalised only when used as a cache KEY, so every raw
   permutation of a configuration got its own entry in the counts and canonical
   caches. Raw states vastly outnumber canonical ones. Fixed by canonicalising
   immediately after `set_leaf`.
3. `canonical()` sorted children with `key=repr`, building a string per
   comparison. All three children of a node have equal height and so compare
   directly; the repr was removable overhead.
The module-level caches are now also cleared between fixed-point iterations.

**Congruence-family hill climbing from random weight vectors.**
Weak, not dead. Reached 196 on (s=1,n=12) where VT gives 316. The evaluator is
correct — fed the exact VT weights (1,…,12) mod 13 it returns 316, verified by
verify.py — so the failure is purely in the search over weight space. Needs
structured seeding (arithmetic progressions, Helberg-Ferreira recurrences) rather
than uniform random starts before it is worth more compute.

**Summary abstraction for k=5: memoise ot[] on the record triple (r0,r1,r2).**
Dead, and instructive about testing depth. The hypothesis was that a stable
configuration's optimum depends only on its three count fields, not on its child
structure — which would collapse the state space. It **holds exactly at k=2 and
k=3** (0 conflicts) and **fails at k=4**: 2086 summary collisions, 1750 with
conflicting ot values. Testing only up to k=3 would have produced a confident
wrong answer. Worse, it would not have helped anyway: at k=4 the summary maps
246,792 configurations to 244,706 distinct triples, i.e. 0.85% compression.

**Lazy evaluation of ot[] from the root.**
Dead. MNSX compute ot[a] eagerly for every configuration but return only the
root's, so on-demand evaluation looked like an easy win. Measured: **100% of
configurations are reached from the root** at k = 2, 3 and 4 alike. Every stable
configuration is reachable by some query sequence, so there is nothing to skip.
Lazy and eager agree on the root value, confirming the instrumentation was right.

**Admissible-heuristic branch and bound for k=5.**
Dead. The valid bound V(sigma) <= 2^(-k) * sum over unqueried v of
Pr[v sensitive | sigma] is non-negative termwise, so it fires only where every
unqueried leaf has sensitivity zero — exactly Rule A's trigger. Prunes nothing.

**Information-revealing relaxations for k=5.**
Dead, vacuously. Telling the tree which child of the root is the minority — the
mildest such relaxation — lets it query both majority children entirely, which
harvests all 2^k sensitive bits (they all live there) while never touching the
minority child: P_q = 2^k, P_m = 0, rho = 1 > 0 for every alpha, so alpha_hat is
infinite. Generalises: alpha_k measures the tension between harvesting sensitive
bits and avoiding the minority, so ANY hint about the minority's location
dissolves it. A sound relaxation must act on the configuration space, never on
the tree's information.

**Odd-cycle cutting planes to beat the published LP bounds at (2,13)+.**
Dead in one second of compute, and instructively so. At the clique-LP optimum of
G(10,2) the separation found ZERO violated odd-cycle inequalities. Cause: the LP
optimum spreads mass diffusely (typical x_v ~ 0.02), and odd-cycle cuts only cut
solutions with values concentrated near 1/2 along a cycle. The slack theta
exposes (20.47 -> 18.64 at (2,10)) lives in PSD structure invisible to cycle
cuts. Machinery kept (`search/cutplane.py`, with exact certification incl.
combinatorial re-validation of every cycle) since it may bite on OTHER MIS-type
problems, but for deletion graphs: certified UB improvements at n >= 13 need
SDP-scale tools, which do not fit this machine at 8192 vertices. Also checked:
the recursive bound A_s(n) <= 2*A_s(n-1) (valid, proof: append the shared last
bit to a common subsequence of the prefixes) gives 106 at (2,13) -- worse than
the LP's 93. No cheap certified route below the published numbers exists here.

**Beat 113 (quaternary single-edit, n=6) by raw ILS or symmetric-invariant search.**
Both fell short, with useful structure. Raw ILS: 6 workers x 25 min -> 111
(45s/seed gave 109; the last units do not fall to short local search against a
24h ReduMIS incumbent). All-subgroup sweep (30 subgroups of S_4, +/- reversal):
best invariant code **112, PROVEN optimal within the A_4-invariant class** —
one below the record. Every other proven class tops out far lower; the classes
that timed out (|H| = 2, 3, 4) reached only 86-106 as incumbents. Standing
insight (CORRECTED 2026-08-25: A_4-invariant sizes are 12b or 4+12b, so 113/114
are impossible in that class on arithmetic grounds and the next attainable value
above 112 is 120 -- "costs one codeword" was the wrong reading), and the record
itself is provably not invariant under any subgroup class we solved to
optimality. Follow-up in flight: symmetry-seeded ILS from the proven-112 code.

**Beat 352 (quaternary single-edit, n=7) by raw ILS.**
Fell short: 6 workers x 90 min -> 349 (seeds 343-349, all verified). Within
0.9% of the 24h-ReduMIS incumbent but no record. The softness argument (same
24h budget on a 4x graph) was directionally right — we got proportionally
closer than at n=6 — but not enough without a stronger engine or much more time.

**Beat the quaternary records with the record-setting engine itself (KaMIS).**
Record attempts closed. Built ReduMIS/OnlineMIS locally (clang + conda OpenMP,
four portability patches), exported our validated graphs to METIS, ran 14+
seeds per cell. Outcome: n=6 record 113 MATCHED three times, n=7 record 352
MATCHED once (seed spread 348-352), n=8 reached 1129 vs published 1130. Never
exceeded anywhere. Key discovery: ReduMIS CONVERGES in ~12 minutes on these
graphs — the published "24h" is a timeout label, not depth — so the records are
best-of-3 draws from a quick-convergence distribution, and 14 draws without an
exceedance suggests the ceiling is real. OnlineMIS is strictly weaker here
(342-345 at n=7). Do not retry with more seeds expecting a different outcome;
any future attempt needs a genuinely different method, not more draws.
CORRECTION to the entry above, same day: "do not retry with more seeds" is
earned at n=6 and n=7 (14 draws, records matched, never exceeded). At n=8 only
TWO seeds were drawn (1129, 1123) against a record of 1130, and the n=7 seed
spread was 4 wide — so a final bounded batch at n=8 is a different bet, not a
re-roll. Running 4 more seeds as the closing action; if none exceeds 1130, the
campaign conclusion stands unmodified.
CLOSED (recorded 2026-08-13, from results/kamis_n8_s*.log): the closing batch
DID exceed. Nine n=8 seeds ran, spread 1123-1131, and seeds 4 and 9 both
reached **1131** — one above W-H's 1130. So the "never exceeded" conclusion is
earned at n=6 and n=7 only; at n=8 the published 1130 was a best-of-3 sampling
artifact, not a ceiling, and the hedge above was the right call. Still not a
record: the literature's 1132 (Ashlock-Houghten-Brown-Orth 2012) stands.
Artifact best/edit4_n8_size1131.txt re-verified in a fresh process
(sha256 483b9385aed78e8695dda5e0b721d4620ce4f9dcb07ea844ad8f4e5242ec824c):
verify_edit4 True/1131, and zero violating pairs under an independent pairwise
Levenshtein >= 3 check. writeup.md section 6 was stale on this row until now.

**Computing alpha_5 for the 3-MAJ lower bound — ABANDONED, closed by mathematics.**
Final status after adversarial review (2026-08-07): not just infeasible (~79
years at 10^15 configurations) but pointless even if free. Goos-Jayram (CCC
2016) prove R(3-MAJ_h) >= Omega(2.5962^h) via conical-junta degree, and every
defensible extrapolation of the alpha_k ceiling (2.575-2.580; certified only
<= 2.649) sits below 2.5962. The alpha_k program cannot reach a bound the field
has had since 2016. Our own earlier claim "2.57143 unchanged since 2013" was
wrong — the first status check missed Goos-Jayram; caught by the workflow's
deeper literature sweep. Lesson reinforced: one ad-hoc search is not a status
check.

**Beat Blaeser's 19 for <3,3,3> over F_3 via the Wang framework (generalized by us).**
Fell short, with the full negative documented in tensor-rank-f3/README.md. The
port itself passed a six-cell regression gate (byte-identical F_2 certificate;
Table 9 orbit censuses reproduced at F_3 including 736 for 3x3; de Groote's 7
at F_3 <2,2,2>), so the machinery is right; the RESULT at the paper's own
search parameters after 25.5 h on 10 cores is R >= 18 (verified), one below
the 2003 classical bound. Cause: the framework's two strongest techniques
degrade with q (forced-product enumerations q^{s(t-s)}; substitution branching
(q^d-1)/(q-1) per level under a fixed step cap). Escalating to the author's
100M two-pass regime projects to WEEKS on this hardware with no assurance of
clearing even 19. Do not retry at laptop scale; a genuinely new technique or
serious compute would be needed. The durable outputs: the field-generic port
(answers the paper's "applies verbatim" question: code yes, bound-power no),
four verified F_3 certificates, and the q-degradation calibration table.

**The all-subgroup sweep at quaternary n=6 never covered 15 of its 60 cells.**
Not a failed idea but a hole in a recorded conclusion, found 2026-08-13. The
sweep in search/edit4.py guards its exact MILP with --exact-cap (default 900
orbits) and SKIPS any cell above it. At n=6 that silently drops 15 of the 60
(subgroup x reversal) cells, and they are exactly the weakest-symmetry ones:
|H|=1 (both), all 9 order-2-without-reversal, all 4 order-3-without-reversal.
So "best invariant 112, and the record is provably not invariant" was never a
statement about a single alphabet involution or a single 3-cycle -- those were
never evaluated at all. The "|H| = 2, 3, 4 timed out at 86-106" line recorded
earlier refers to the REVERSAL-PAIRED versions (groups of order 4 and 6), which
are strictly stronger constraints. Separately, where the MILP was pushed past
the cap it degrades visibly: a 1021-orbit cell in
results/edit4_sweep_n6_20260806T025537.json returns an incumbent of 8.

**Filling those 15 cells with a hand-rolled (omega,1) local search.**
Failed, and the calibration cell is why it was caught. Best over all 15 cells
was 100 (vs 113 W-H, 114 literature) -- but the same solver scores only 99 on
the TRIVIAL cell, which IS the unrestricted problem where 111 (repo ILS) and
113 (KaMIS) are known attainable. A solver at 88% of a known-attainable value
measures itself, not the mathematics, so none of the 15 numbers meant anything.
Cause, and the durable lesson: the "insert v when it outweighs every selected
neighbour" rule cannot express "remove one, add two". (1,2)-swaps and plateau
moves are what MIS local search runs on, and search/ils.py:Searcher has both --
which is precisely why it reaches 111 where the rule-only version stalls at 92.
Always put an instance of KNOWN optimum in the batch; the trivial subgroup is
free and it invalidated the whole run in one line.
Log: results/edit4_orbit_ls_n6_20260813T190813.json.

**Certified theta upper bound for 1zc.1024 (Sloane "Most wanted", [112,117] since 1998).**
Dead in 94 seconds of compute. Plain Lovasz theta = 128.67 float — ELEVEN above
the 1998 combinatorial UB of 117. Calibration at n=7 (theta 20.67 vs alpha 18,
gap 2.7) misled: the gap grows with n because the classical asymmetric-code
bounds (Kloeve/Borden weight-class partitioning) exploit structure invisible to
the vanilla SDP. theta-plus/Schrijver refinements typically gain 1-3, nowhere
near 11.7. The graphs and the machinery stay (independently rebuilt DIMACS
instances, edge-for-edge; ReduMIS reproduces 18/36/62/112 exactly), but the UB
route for this family needs weight-structured SDPs, not theta. Redirecting to
EXACT closure: at 16,640 edges the 1024-vertex instance is 10x sparser than the
deletion graphs where our MILP stalled — worth one serious attempt.

**Exact closure of 1zc cells via generic MILP.**
Dead in one probe: edge-formulation HiGHS could not close even 1zc.256 (256
vertices, 2816 edges, known alpha=36) within 600 s, so 1zc.1024 is far out of
reach this way. Consistent with the deletion-graph experience: edge
formulations of MIS are weak, and the 2005-era closures of comparable
challenge cells (Borchers, 1dc.1024, 298 CPU-hours) used dedicated
branch-and-bound with clique-cover bounds, i.e. a specialized exact solver we
do not have in the toolchain. The 1zc thread closes with three documented
negatives (theta +11.7 off, MILP, and the UB structure argument); the durable
outputs are the independently-rebuilt instances, METIS exports, and exact
ReduMIS reproduction of all four known lower bounds (18/36/62/112).

**Beat 113 at quaternary n=6 via the weak-symmetry classes the sweep skipped.**
Fell short at 112, but the cell is NOT closed and the reason is worth keeping.
Discovery first: `search/edit4.py --sweep` guards each orbit MILP with
`--exact-cap` (default 900 orbits), and at n=6 that cap silently SKIPS 15 of
the 60 (subgroup x reversal) cells -- exactly the weakest-symmetry ones (all 9
order-2-without-reversal, all 4 order-3-without-reversal, and |H|=1 with/without
reversal). So the recorded "best invariant = 112, proven optimal in the A_4
class" never covered them; the "86-106 incumbents" in the entry above are the
REVERSAL-PAIRED order-2/3 cells, i.e. groups of order 4 and 6. A single alphabet
involution had never been tested.
Conjugacy collapses those 15 cells to 4 distinct problems (if H' = sHs^-1 then
s acting symbol-wise is a graph automorphism carrying H-invariant codes onto
H'-invariant ones). Corroborated by the censuses: all 6 transposition cells
build identical 2080/1359/96182 orbit graphs.
First attempt was solver-limited and its own internal check caught it: values
inside a conjugacy class must be EQUAL and came back with spread 5 (88,87,86,
84,88,83). Cause: an "(omega,1) insertion" move set cannot express remove-one-
add-two, which is what MIS local search runs on. Porting the (1,2)-swap +
plateau moves from the repo's own Searcher fixed it (that class went 98 -> 108
in 8 seconds).
Calibrated result (`search/edit4_orbit_ls.py --batch`, 4 seeds x 15 min x 4
workers, results/edit4_orbit_batch_n6_20260814T013427.json). Gate: the trivial
cell (= unrestricted problem) reached 110 at 60 core-min vs the repo ILS's 111
at 150 core-min, so the solver is comparable to a proven engine, and lands
~1-3 short of attainable.

| cell | usable orbits | best |
|---|---|---|
| trivial (calibration) | 4096 | 110  (bar 111 / KaMIS 113) |
| reversal only | 1498 | 107 |
| transposition (0 1) | 1359 | 95 |
| double transposition (0 1)(2 3) | 2046 | **112** |
| 3-cycle (0 1 2) | 1295 | 109 |

STRUCTURAL FINDING, cheap and checkable: two of these four classes cannot
contain a code of size 113 AT ALL. Under (0 1)(2 3) the involution is
fixed-point-free on words -- all 2048 orbits have size exactly 2 -- so every
invariant code has EVEN size. Under the 3-cycle, orbits have size 3 except the
single fixed word 333333, so sizes are 0 or 1 mod 3, and 113 = 2 mod 3. In both
classes the next attainable value above 112 is exactly 114, the literature
record. So the (0 1)(2 3) cell at 112 is not "one short of 113"; 113 is
impossible there and the question is 112 vs 114.
STILL OPEN, and the sharpest remaining target in this repo: that cell reduces to
ONE unweighted MIS instance -- 2046 vertices, need alpha >= 57 (=> 114) against
the 56 found. Half the size of the 4096-vertex graph where ReduMIS gets 113.
Not attempted: exact closure with clique-cover bounds. A proven 56 closes the
class; a proven 57 ties the 2012 record with a machine-verifiable artifact
(their codeword lists were never published).

**Transfer the diagonal-symmetry method from n=6 to n=7 (beat E4(7,3)=356).**
Fell one short. All 30 of the 33 conjugacy classes of S_4 x C_2 swept (the enumeration bug was found later; see sota.md) with the calibrated
orbit local search, 4 seeds x up to 600 s each, 3.2 h wall
(results/edit4_n7_diag.json, search/edit4_n7.py). Best **355** in an |H|=6 class
(2218 usable orbits), verified. Not an arithmetic wall: that class can attain
357, and 10 of the 30 classes admit a +1, so this is a search shortfall. Context
that keeps it honest: 355 still beats W-H Table 3's 352 by three and beats every
prior attempt here (raw ILS 349, local KaMIS 352), so the method is the
strongest we have at n=7 -- just not enough. The n=6 win (116 vs 114) does NOT
generalise for free. Untested: whether a 10x compute increase closes the 1-gap;
the free order-4 classes need +4 (360) and are the wrong place to spend it, the
mixed-orbit-size classes needing +1 are the right one.

**Reporting 112 as the ceiling of the free-acting V_4 class at n=6.**
Wrong, and it cost four codewords. The sweep gave each class 4 seeds; V_4 returned
112 on all four and was recorded as topping out there, while the two classes that
returned 116 were later given 10 more seeds each. Nobody pushed V_4. An
adversarial re-run with a stronger operator (ARW 2-improvements + plateau) found
116 on every seed and 120 on two in six; our own operator at 10 seeds finds 116 on
five of ten and never 120 (results/repro_v4_120.json). Lesson: a four-seed result
is a statement about the search, and recording it as a property of the cell is the
same failure this repo already documented once at the parameterisation level.

**Close the free order-4 classes at n=6 (120 vs 124) with Lovasz theta on the orbit graph.**
Dead in one run (2026-09-07): theta(V_4 orbit graph, 1021 vertices) = 45.1 float
(SCS, 39 min; results/theta_edit4_n6_orbit_class20.json) against the 31 orbits
that 124 needs. Same failure shape as the 1zc.1024 attempt: vanilla theta is
loose by a large factor on these graphs, so certifying it would prove nothing.
The class stays open at 120 with no upper-bound evidence either way.

**Long single seeds of the C ARW engine (search/mis_ils.c).**
Wrong shape for this operator, measured twice. n=6 free classes: 15 seeds x 600 s
all 116, while 30 seeds x 60 s gave 120 four times (results/edit4_cmis_n6_push124_*,
results/edit4_cmis_n6_shortbatch.log). n=7 order-6 class: 2 seeds x 600 s gave
358, 358; 3 seeds x 900 s gave 355, 355, 355. The incumbent-return rule kept the
walk in one basin; a fresh-restart branch was added (--restart), but the
operational rule stands: many short seeds beat few long ones here.

**Close the remaining E4(6,4) classes at 30 by MaxSAT (RC2).**
Plateau, recorded 2026-09-07. RC2 (python-sat, glucose4 backend, adapt/exhaust/minz)
proves optimum 30 on the 518-orbit order-4 diagonal class in 21 s
(results/editq_n6_d4_exact_class23.json), but on the 1498-orbit class 30
(478k hard clauses) it had not finished after 2 h 8 min and was killed; the
2016-orbit class 29 was never started. Core-guided MaxSAT scales badly when the
optimum is a small fraction of the soft weight (30 of 1952). Whether 31 exists
in those two classes, or in the unrestricted 4096-vertex graph (table UB 32),
stays open; a clique-cover branch-and-bound or a SAT "at least 31" encoding
with symmetry breaking is the next rung, not more RC2 time.

**RC2 MaxSAT closure of the (7,4) order-8 diagonal class (1277 orbits, 523k edges).**
Killed after 75 minutes without converging (2026-09-08), the same stall as the
1498-orbit (6,4) class; the 265-, 518- and 708-orbit classes closed in
seconds. Core-guided MaxSAT on these instances scales to a few hundred orbits
only. Whether 72 is optimal in that class stays unknown.

**E4(9,6) in the V_4 class.**
Not runnable with the current tools (2026-09-08): the orbit graph has 58,435
usable orbits and 1.51 billion edges (a 16 GB edge list, density 88%), too
dense for the exact clique route (~200M compatible pairs) and too large for the
C engine's adjacency arrays. The order-8 diagonal class (5813 orbits) is
exactly closed at 32; a V_4 attempt would need an implicit-adjacency search
(conflict test by Levenshtein on the fly) rather than an explicit edge list.
The 16 GB file was deleted; the orbit list is kept.

**Exact closure of the (7,5) reversal class by max clique of its compatibility graph.**
Not feasible with networkx (2026-09-08): 2248 orbits with 957,612 compatible
pairs (43% dense complement); max_weight_clique ran 2 h 11 min without
finishing and was killed. The clique route works when the compatibility graph
is sparse (a few hundred thousand pairs on ≤ 12k vertices: (7,5) and (8,6)
free classes, (9,6)/(9,7) order-8 class), not at this density.

**NuMVC (vertex-cover local search with edge weighting and configuration checking) as a second operator.**
Calibration on the n=6 V_4 orbit graph (1021 vertices, density 19%), 60 s
per seed: 108, 108 after the max-score fix (96, 92 with sampled removal),
against 116-120 for the ARW (1,2)-swap engine and 116 for ReduMIS. NuMVC's
strength is very large sparse graphs; these orbit graphs are small and dense.
Kept (search/mis_numvc.c) because the n=9 orbit graphs are far sparser (0.7%
density) and it may behave differently there; that run is separate.
NuMVC on the n=9 C_4-product orbit graph (65,533 vertices, 0.7% density),
2 seeds x 900 s: 3524, 3516 (881/879 orbits) against ARW 3688-3724 and
ReduMIS 3712-3720 on the same graph. Sparser did not help; NuMVC is not the
right operator for these instances at any size tried.
