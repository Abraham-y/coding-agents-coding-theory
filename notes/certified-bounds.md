# Certified upper bounds on A_s(n)

Each row is the floor of an **exact rational** LP dual objective, validated
against all 2^n dual constraints in `fractions.Fraction` arithmetic. The float
LP only proposes the dual vector; it is never trusted (CLAUDE.md rules 10, 12).
Certificates are in `results/cert_n<N>_s<S>.json` and are re-checkable with
`python search/lp_bound.py --verify-cert <file>`.

| s | n | certified UB | exact dual objective | published LP UB | best known code | agree |
|---|---|---|---|---|---|---|
| 1 | 6 | 10 | 205003/20000 | 10 | 10 | yes |
| 1 | 7 | 17 | 1715023/100000 | 17 | 16 | yes |
| 1 | 8 | 30 | 303237/10000 | 30 | 30 | yes |
| 1 | 9 | 53 | 5356523/100000 | 53 | 52 | yes |
| 1 | 10 | 96 | 4826261/50000 | 96 | 94 | yes |
| 1 | 11 | 175 | 17519677/100000 | 175 | 172 | yes |
| 1 | 12 | 321 | 2008033/6250 | 321 | 316 | yes |
| 1 | 13 | 593 | 14838071/25000 | 593 | 586 | yes |
| 2 | 6 | 4 | 4/1 | 4 | 4 | yes |
| 2 | 7 | 5 | 568903/100000 | 5 | 5 | yes |
| 2 | 8 | 8 | 103151/12500 | 8 | 8 | yes |
| 2 | 9 | 12 | 64533/5000 | 12 | 11 | yes |
| 2 | 10 | 20 | 1023461/50000 | 20 | 16 | yes |
| 2 | 11 | 33 | 1670329/50000 | 33 | 24 | yes |
| 2 | 12 | 55 | 2773999/50000 | 55 | 37 | yes |
| 2 | 13 | 93 | 233707/2500 | 93 | 56 | yes |
| 3 | 6 | 2 | 2/1 | 2 | 2 | yes |
| 3 | 7 | 2 | 2/1 | 2 | 2 | yes |
| 3 | 8 | 4 | 4/1 | 4 | 4 | yes |
| 3 | 9 | 5 | 137501/25000 | 5 | 5 | yes |
| 3 | 10 | 7 | 188261/25000 | 7 | 6 | yes |
| 3 | 11 | 10 | 34063/3125 | 10 | 8 | yes |
| 3 | 12 | 16 | 1627641/100000 | 16 | 12 | yes |
| 3 | 13 | 24 | 246037/10000 | 24 | 16 | yes |

All certificates re-checked: 25/25 PASS (24 plain cells + the (2,12) clique-cuts certificate; recount 2026-08-09).

---

## Semidefinite (Lovasz theta) bounds — strictly stronger than the LP

Added 2026-08-05, prompted by reading OpenAI's Astra formalisation (2026-08-01).
Their binary-codes result chains `binaryRate <= combinedVariationalRate < mrrwRate`,
i.e. it beats MRRW — which IS the Delsarte LP bound — by moving to Gram-matrix /
positive-semidefinite arguments. Their setting is the Hamming metric in the
asymptotic rate regime and does not touch these deletion cells, but the
structural lesson transfers: **where an LP bound is loose, the semidefinite
strengthening is the next rung.**

`search/theta.py` computes the Lovasz theta bound, with the SDP used ONLY to
propose a dual; the reported number comes from rounding that dual to exact
rationals, inflating until the matrix is exactly positive semidefinite (verified
by rational LDL^T), and taking the floor.

| cell | exact value | certified LP bound | **certified theta bound** |
|---|---|---|---|
| (s=2, n=6) | 4 | 4 | **4** |
| (s=2, n=7) | 5 | 5 | **5** |
| (s=2, n=8) | 7 | 8 | **7**  <- strictly better than LP, and tight |

### A floating-point trap this caught, worth recording

The float SDP returns theta(G(6,2)) = **3.999997**, where the true value is 4.
Flooring that float gives 3 and would assert A_2(6) <= 3, which is **false** —
the true value is 4. Any pipeline that floors a solver's float here silently
emits a wrong bound. This is precisely the failure CLAUDE.md rule 10 forbids, and
it is not hypothetical: it happens on the very first cell tried. Only the
exactly-certified dual is ever floored in this repo.

### Scaling

Plain theta is expensive: 64 vertices took 2.8s, 128 took 80s, 256 took 467s,
and 512 did not finish inside 10 minutes. The formulation carries one variable
per edge (17,183 already at n=8; 1.46M at n=12), so the open cell (2,12) is far
out of reach this way. Symmetry reduction is the obvious next step, but the
automorphism group here is only the Klein four-group, worth at most a 4x block
reduction — not enough on its own for 4096 vertices.

### Theta additions (2026-08-06, after the n=9 certifications finished)

| cell | exact | certified LP | certified theta | note |
|---|---|---|---|---|
| (s=2, n=9) | 11 | 12 | **11** | tight; third cell where theta strictly beats the LP |
| (s=1, n=9) | 52 | 53 | 53 | theta ~53.03 fails to beat even the LP here |

Float-only theta (informative, NOT certified): (2,10) 18.64 [exact 16, LP 20],
(1,10) 95.98 [exact 94, LP 96]. Pattern: theta is tight only on small cells;
Hypothesis T remains refuted. Certification cost is the binding constraint:
the (1,9) run took 9.5 h, nearly all in the exact rational LDL^T.
