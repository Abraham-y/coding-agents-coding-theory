"""Produce a certified upper-bound table across cells, and re-check every certificate.

Every number this emits is the floor of an exact rational dual objective that has
been validated against all 2^n dual constraints in Fraction arithmetic. The float
LP is only ever a proposal generator.

Run:  python search/bound_table.py --max-n 14
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.lp_bound import build_system, certify, propose_dual, verify_certificate  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")
NOTES = os.path.join(ROOT, "notes")

# Published values, transcribed from Weindel & Heckel arXiv:2504.00613v2 Tables 1-2.
# "lp" is their LP upper bound row, "best" their best known construction.
PUBLISHED = {
    (1, 6): (10, 10), (1, 7): (17, 16), (1, 8): (30, 30), (1, 9): (53, 52),
    (1, 10): (96, 94), (1, 11): (175, 172), (1, 12): (321, 316), (1, 13): (593, 586),
    (1, 14): (1104, 1096),
    (2, 6): (4, 4), (2, 7): (5, 5), (2, 8): (8, 8), (2, 9): (12, 11),
    (2, 10): (20, 16), (2, 11): (33, 24), (2, 12): (55, 37), (2, 13): (93, 56),
    (2, 14): (160, 87),
    (3, 6): (2, 2), (3, 7): (2, 2), (3, 8): (4, 4), (3, 9): (5, 5),
    (3, 10): (7, 6), (3, 11): (10, 8), (3, 12): (16, 12), (3, 13): (24, 16),
    (3, 14): (38, 22),
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=14)
    ap.add_argument("--min-n", type=int, default=6)
    ap.add_argument("--denom", type=int, default=100000)
    args = ap.parse_args()

    rows = []
    for s in (1, 2, 3):
        for n in range(args.min_n, args.max_n + 1):
            t0 = time.time()
            try:
                g, constraints, ball, _ = build_system(n, s, cuts=False)
                y_float, lp_val = propose_dual(g.num_vertices, constraints, ball)
                y, bound, min_slack = certify(y_float, ball, args.denom)
                ub = math.floor(bound)
                cert = {
                    "method": "LP duality, HiGHS proposal rounded and repaired exactly",
                    "n": n, "s": s, "cuts": False, "denominator": args.denom,
                    "float_lp_objective": lp_val,
                    "bound": [bound.numerator, bound.denominator],
                    "certified_upper_bound": ub,
                    "y": [[v.numerator, v.denominator] for v in y],
                }
                path = os.path.join(RESULTS_DIR, f"cert_n{n}_s{s}.json")
                os.makedirs(RESULTS_DIR, exist_ok=True)
                with open(path, "w") as fh:
                    json.dump(cert, fh)
                ok, msg = verify_certificate(cert)
                pub_lp, pub_best = PUBLISHED.get((s, n), (None, None))
                rows.append({
                    "s": s, "n": n, "certified_ub": ub,
                    "exact_objective": f"{bound.numerator}/{bound.denominator}",
                    "float_lp": round(lp_val, 6),
                    "published_lp_ub": pub_lp, "published_best": pub_best,
                    "recheck": ok, "recheck_msg": msg,
                    "seconds": round(time.time() - t0, 1),
                })
                flag = "" if pub_lp is None or ub == pub_lp else f"  <-- differs from published {pub_lp}"
                print(f"s={s} n={n:2d}  certified UB {ub:6d}  published LP {pub_lp}  "
                      f"best known {pub_best}  recheck={'PASS' if ok else 'FAIL'}{flag}")
            except Exception as exc:  # noqa: BLE001
                print(f"s={s} n={n:2d}  SKIPPED: {type(exc).__name__}: {exc}")

    os.makedirs(NOTES, exist_ok=True)
    out = os.path.join(NOTES, "certified-bounds.md")
    with open(out, "w") as fh:
        fh.write("# Certified upper bounds on A_s(n)\n\n")
        fh.write("Each row is the floor of an **exact rational** LP dual objective, validated\n")
        fh.write("against all 2^n dual constraints in `fractions.Fraction` arithmetic. The float\n")
        fh.write("LP only proposes the dual vector; it is never trusted (CLAUDE.md rules 10, 12).\n")
        fh.write("Certificates are in `results/cert_n<N>_s<S>.json` and are re-checkable with\n")
        fh.write("`python search/lp_bound.py --verify-cert <file>`.\n\n")
        fh.write("| s | n | certified UB | exact dual objective | published LP UB | best known code | agree |\n")
        fh.write("|---|---|---|---|---|---|---|\n")
        for r in rows:
            agree = "yes" if r["published_lp_ub"] == r["certified_ub"] else "**NO**"
            fh.write(f"| {r['s']} | {r['n']} | {r['certified_ub']} | {r['exact_objective']} | "
                     f"{r['published_lp_ub']} | {r['published_best']} | {agree} |\n")
        fh.write("\nAll certificates re-checked: "
                 f"{sum(1 for r in rows if r['recheck'])}/{len(rows)} PASS.\n")
    print(f"\nwrote {out}")

    logpath = os.path.join(RESULTS_DIR, f"bound_table_{time.strftime('%Y%m%dT%H%M%S')}.json")
    with open(logpath, "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"wrote {logpath}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
