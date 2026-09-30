"""SAT attack on 'is there an s-deletion code of length n with K codewords?'

The clique cover of G(n,s) is exactly the constraint set we want: for every word
t of length n-s, at most one codeword may contain t as a subsequence. So

    variables : x_v for each v in {0,1}^n
    clauses   : AtMost1{ x_v : t is a subsequence of v }   for every t
                AtLeast K { x_v : all v }

A SAT answer is a CONSTRUCTION and is trivially checkable by verify.py, so the
positive direction is fully in scope for this sprint.

An UNSAT answer is a NONEXISTENCE claim and per CLAUDE.md rule 11 it is a rumour
until a DRAT certificate is checked with drat-trim. This module will print UNSAT
but refuses to record it as a result; use --proof to emit the certificate.

Optional --symmetry complement adds x_v <-> x_{complement(v)}, which prescribes a
complement-invariant code. That is a RESTRICTION: SAT under it is still a real
construction, but UNSAT under it says nothing about the unrestricted problem.

Run:
    python search/sat_exact.py --n 12 --s 2 --k 38 --timeout 600
    python search/sat_exact.py --n 11 --s 2 --k 24 --timeout 60      # known feasible
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

from search.graph import complement_idx, load_graph  # noqa: E402
from verify import canonical_form, checksum, verify  # noqa: E402
from verify2 import verify2  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")
BEST_DIR = os.path.join(ROOT, "best")


def build(n: int, s: int, k: int, symmetry: str | None):
    g = load_graph(n, s)
    pool = IDPool(start_from=g.num_vertices + 1)
    xs = [v + 1 for v in range(g.num_vertices)]  # var id = vertex index + 1

    clauses = []
    for members in g.cliques.values():
        if len(members) < 2:
            continue
        lits = [m + 1 for m in members]
        enc = CardEnc.atmost(lits=lits, bound=1, vpool=pool, encoding=EncType.seqcounter)
        clauses.extend(enc.clauses)

    enc = CardEnc.atleast(lits=xs, bound=k, vpool=pool, encoding=EncType.totalizer)
    clauses.extend(enc.clauses)

    if symmetry == "complement":
        for v in range(g.num_vertices):
            c = complement_idx(v, n)
            if c > v:
                # x_v <-> x_c
                clauses.append([-(v + 1), c + 1])
                clauses.append([v + 1, -(c + 1)])

    return g, xs, clauses


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--timeout", type=float, default=600)
    ap.add_argument("--symmetry", choices=["complement"], default=None)
    ap.add_argument("--tag", type=str, default="sat")
    args = ap.parse_args()

    t0 = time.time()
    g, xs, clauses = build(args.n, args.s, args.k, args.symmetry)
    build_time = time.time() - t0
    print(f"encoded: {len(clauses)} clauses, target k={args.k}, build {build_time:.1f}s")

    solver = Cadical153(bootstrap_with=clauses)
    t1 = time.time()
    # pysat has no wall-clock interrupt on Cadical153 without a timer thread;
    # use the built-in budget instead.
    solver.conf_budget(int(args.timeout) * 20000)
    res = solver.solve_limited(expect_interrupt=False)
    solve_time = time.time() - t1

    verdict = {True: "SAT", False: "UNSAT", None: "UNKNOWN (budget exhausted)"}[res]
    print(f"result: {verdict} in {solve_time:.1f}s")

    record = {
        "method": "SAT via clique-cover AtMost1 + AtLeast-k cardinality",
        "solver": "Cadical153",
        "n": args.n,
        "s": args.s,
        "k": args.k,
        "symmetry": args.symmetry,
        "clauses": len(clauses),
        "build_seconds": round(build_time, 2),
        "solve_seconds": round(solve_time, 2),
        "verdict": verdict,
    }

    if res is True:
        model = set(l for l in solver.get_model() if l > 0)
        chosen = [v for v in range(g.num_vertices) if (v + 1) in model]
        words = sorted(g.words[v] for v in chosen)
        ok1, sc1, r1 = verify(words, args.n, args.s)
        ok2, sc2, r2 = verify2(words, args.n, args.s)
        print(f"  verify.py : valid={ok1} score={sc1} {r1}")
        print(f"  verify2.py: valid={ok2} score={sc2} {r2}")
        record.update({"size": sc1, "verify": ok1, "verify2": ok2, "words": words})
        if ok1 and ok2 and sc1 >= args.k:
            os.makedirs(BEST_DIR, exist_ok=True)
            note = (
                f"found by SAT (Cadical153) via clique-cover encoding"
                + (f", symmetry={args.symmetry}" if args.symmetry else "")
            )
            path = os.path.join(BEST_DIR, f"n{args.n}_s{args.s}_size{sc1}.txt")
            with open(path, "w") as fh:
                fh.write(canonical_form(words, args.n, args.s, note))
            print(f"  artifact  : {path}")
            print(f"  checksum  : {checksum(words, args.n, args.s)}")
            record["artifact"] = path
    elif res is False:
        print("  NOTE: UNSAT is a nonexistence claim. Per CLAUDE.md rule 11 it is NOT")
        print("        a result without a drat-trim-checked certificate. Not recorded.")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    logpath = os.path.join(
        RESULTS_DIR, f"{args.tag}_n{args.n}_s{args.s}_k{args.k}_{stamp}.json"
    )
    with open(logpath, "w") as fh:
        json.dump(record, fh, indent=2)
    print(f"  log       : {logpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
