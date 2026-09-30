# SAT certificate: A_4-invariant quaternary single-edit code at n=6 is at most 112

Produced 2026-08-06 by an adversarial-review agent, independently of the
original HiGHS MILP pipeline (whose [OPT] statuses proved hard to reproduce and
should be treated as fragile evidence).

Reduction: an A_4-invariant code is a union of complete A_4 orbits; on the
usable-orbit graph (270 orbits: 269 of weight 12, 1 of weight 4; 17,992 edges;
every conflict edge joins two weight-12 orbits) the maximum equals 4 + 12*m
where m = max independent set among weight-12 orbits.

- `no10orbits.cnf`  — DIMACS encoding of "m >= 10" (sha256 a6587790...)
- `no10orbits.drup` — CaDiCaL refutation, VERIFIED by drat-trim (sha256 815f38a9...)

m >= 10 is UNSAT (certified) and m = 9 is achieved by a verified witness, so the
maximum invariant code has exactly 4 + 12*9 = 112 codewords. Confirmed
independently by RC2 exact MaxSAT (optimum 112, witness re-verified).

Artifact: ../edit4_n6_A4invariant_112.txt
