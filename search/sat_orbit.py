"""SAT on the ORBIT graph: is there an H-invariant code of a given size?

Rationale for the retarget: encoding the full 2^n-vertex problem with an
AtLeast-k totalizer produced 2.2M clauses and could not even re-find a known
feasible code at (s=2, n=11, k=24) in 460s. Prescribing complement invariance
halves the variables and, more importantly, halves the cardinality target: a
complement-invariant code of size 38 is 19 orbits, and "independent set of size
19 in a 2011-vertex graph" is a far easier formula than "38 in 4096".

    variables : y_O for each usable orbit O
    clauses   : (-y_O v -y_O')  for every conflicting orbit pair
                AtLeast k { y_O }
    code size : sum of |O| over chosen orbits  (= 2k for the complement group)

SAT here yields a genuine construction, checked by verify.py and verify2.py.
UNSAT here is doubly weak: it is a nonexistence claim (rule 11, needs DRAT) AND
it is restricted to H-invariant codes, so it says nothing about the general case.
Both caveats are printed, and neither is recorded as a result.

Run:
    python search/sat_orbit.py --n 12 --s 2 --group complement --k 19 --timeout 900
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402
from pysat.solvers import Cadical153  # noqa: E402

from search.orbit import OrbitGraph  # noqa: E402
from verify import canonical_form, checksum, verify  # noqa: E402
from verify2 import verify2  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")
BEST_DIR = os.path.join(ROOT, "best")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--group", choices=["complement", "reversal", "both"], default="complement")
    ap.add_argument("--k", type=int, required=True, help="number of ORBITS, not codewords")
    ap.add_argument("--timeout", type=float, default=900)
    ap.add_argument("--tag", type=str, default="satorbit")
    args = ap.parse_args()

    t0 = time.time()
    og = OrbitGraph(args.n, args.s, args.group)
    nv = og.num_vertices
    weights = sorted(set(og.weight))
    print(f"orbit graph: {nv} usable orbits, weights {weights}, {og.dropped} dropped")

    pool = IDPool(start_from=nv + 1)
    ys = [i + 1 for i in range(nv)]
    clauses = []
    edges = 0
    for i in range(nv):
        m = og.adj[i]
        while m:
            low = m & -m
            j = low.bit_length() - 1
            m ^= low
            if j > i:
                clauses.append([-(i + 1), -(j + 1)])
                edges += 1
    enc = CardEnc.atleast(lits=ys, bound=args.k, vpool=pool, encoding=EncType.totalizer)
    clauses.extend(enc.clauses)
    build_time = time.time() - t0
    print(f"encoded: {edges} edges, {len(clauses)} clauses, k={args.k} orbits, build {build_time:.1f}s")

    solver = Cadical153(bootstrap_with=clauses)
    t1 = time.time()
    solver.conf_budget(int(args.timeout) * 20000)
    res = solver.solve_limited(expect_interrupt=False)
    solve_time = time.time() - t1
    verdict = {True: "SAT", False: "UNSAT", None: "UNKNOWN (budget exhausted)"}[res]
    print(f"result: {verdict} in {solve_time:.1f}s")

    record = {
        "method": f"SAT on orbit graph, group={args.group}",
        "solver": "Cadical153",
        "n": args.n,
        "s": args.s,
        "group": args.group,
        "orbit_target_k": args.k,
        "usable_orbits": nv,
        "orbit_edges": edges,
        "clauses": len(clauses),
        "build_seconds": round(build_time, 2),
        "solve_seconds": round(solve_time, 2),
        "verdict": verdict,
    }

    if res is True:
        model = set(l for l in solver.get_model() if l > 0)
        chosen = [i for i in range(nv) if (i + 1) in model]
        words = sorted(og.words_of(chosen))
        ok1, sc1, r1 = verify(words, args.n, args.s)
        ok2, sc2, r2 = verify2(words, args.n, args.s)
        print(f"  code size : {len(words)} from {len(chosen)} orbits")
        print(f"  verify.py : valid={ok1} score={sc1} {r1}")
        print(f"  verify2.py: valid={ok2} score={sc2} {r2}")
        record.update({"size": len(words), "verify": ok1, "verify2": ok2, "words": words})
        if ok1 and ok2:
            os.makedirs(BEST_DIR, exist_ok=True)
            path = os.path.join(BEST_DIR, f"n{args.n}_s{args.s}_size{sc1}_{args.group}.txt")
            with open(path, "w") as fh:
                fh.write(
                    canonical_form(
                        words,
                        args.n,
                        args.s,
                        f"SAT on {args.group}-orbit graph, {len(chosen)} orbits, Cadical153",
                    )
                )
            print(f"  artifact  : {path}")
            print(f"  checksum  : {checksum(words, args.n, args.s)}")
            record["artifact"] = path
    elif res is False:
        print("  NOTE: UNSAT here is (a) a nonexistence claim needing a drat-trim-checked")
        print(f"        certificate per rule 11, and (b) restricted to {args.group}-invariant")
        print("        codes only. It bounds nothing about general codes. Not recorded.")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    logpath = os.path.join(
        RESULTS_DIR, f"{args.tag}_{args.group}_n{args.n}_s{args.s}_k{args.k}_{stamp}.json"
    )
    with open(logpath, "w") as fh:
        json.dump(record, fh, indent=2)
    print(f"  log       : {logpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
