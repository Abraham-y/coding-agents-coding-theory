"""Symmetry-seeded search: extract the proven-optimal A_4-invariant code (112),
then break the symmetry with unconstrained local search. Target: 113+ at n=6.

Rationale. The all-subgroup sweep proved the maximum A_4-invariant code has
exactly 112 codewords -- one below the unconstrained record 113 -- while random
ILS restarts plateau at 109-111. Symmetric optima sit in a different, provably
good basin; perturbing FROM one is the classic route to records that plain
restarts miss.

Run:  python search/edit4_hybrid.py --n 6 --seconds 900 --workers 6
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.edit4 import (  # noqa: E402
    PUBLISHED,
    Edit4Graph,
    all_subgroups_s4,
    orbit_problem,
    orbits_of,
    solve_orbit_exact,
    verify_edit4,
    word_transforms,
)
from search.ils import Searcher  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def a4_optimal_code(n: int, time_limit: int = 600):
    """Solve the A_4 (order-12) orbit MILP exactly; return (words, proven)."""
    g = Edit4Graph(n)
    a4 = next(s for s in all_subgroups_s4() if len(s) == 12)
    transforms = word_transforms(a4, n, with_reversal=False)
    orbits = orbits_of(g, transforms)
    usable, edges = orbit_problem(g, orbits)
    chosen, proven = solve_orbit_exact(usable, edges, time_limit)
    words = sorted(g.words[v] for i in chosen for v in usable[i])
    return g, words, proven


def seeded_run(args):
    n, seed, seconds, seed_words = args
    g = Edit4Graph(n)
    idx = {w: i for i, w in enumerate(g.words)}
    s = Searcher(g, seed)
    for w in seed_words:
        s.add(idx[w])
    s.local_search()
    best = s.current()
    best_size = len(best)
    t0 = time.time()
    since = 0
    while time.time() - t0 < seconds:
        s.perturb(strength=1 + s.rng.randrange(4))
        s.local_search()
        s.plateau(moves=s.rng.randrange(10))
        s.local_search()
        if s.size > best_size:
            best, best_size, since = s.current(), s.size, 0
        else:
            since += 1
            if since > 300:  # restart from the symmetric seed, not from scratch
                s.clear()
                for w in seed_words:
                    s.add(idx[w])
                s.local_search()
                since = 0
    return seed, sorted(g.words[v] for v in best)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--seconds", type=float, default=900)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=40)
    args = ap.parse_args()

    print("extracting proven-optimal A_4-invariant code...")
    g, sym_words, proven = a4_optimal_code(args.n)
    ok, sz, reason = verify_edit4(sym_words, args.n)
    print(f"A_4 code: size={sz}, proven_opt_in_class={proven}, verify={ok} {reason}")
    if not ok:
        return 1

    # keep the symmetric object itself; proven-optimal-in-class is worth having
    os.makedirs(os.path.join(ROOT, "best"), exist_ok=True)
    sym_path = os.path.join(ROOT, "best", f"edit4_n{args.n}_A4invariant_{sz}.txt")
    with open(sym_path, "w") as fh:
        fh.write(f"# quaternary single-edit code, n={args.n}, size {sz}\n")
        fh.write(f"# A_4-invariant, PROVEN optimal among A_4-invariant codes\n")
        for w in sym_words:
            fh.write(w + "\n")
    print(f"symmetric artifact: {sym_path}")

    jobs = [(args.n, args.seed + i, args.seconds, sym_words) for i in range(args.workers)]
    with Pool(args.workers) as pool:
        results = pool.map(seeded_run, jobs)

    pub = PUBLISHED.get(args.n, {}).get("best")
    best_sz, best_words = 0, None
    for seed, words in sorted(results, key=lambda r: -len(r[1])):
        ok, k, reason = verify_edit4(words, args.n)
        print(f"seed {seed}: size={len(words)} verify={ok}" + ("" if ok else f" {reason}"))
        if ok and k > best_sz:
            best_sz, best_words = k, words
    print(f"BEST HYBRID: {best_sz}  (published best known: {pub})")

    stamp = time.strftime("%Y%m%dT%H%M%S")
    with open(os.path.join(ROOT, "results", f"edit4_hybrid_n{args.n}_{stamp}.json"), "w") as fh:
        json.dump({"n": args.n, "seed_code_size": sz, "workers": args.workers,
                   "seconds": args.seconds,
                   "sizes": {str(s): len(w) for s, w in results},
                   "best": best_sz, "published_best": pub}, fh, indent=2)

    if best_words and pub and best_sz > pub:
        path = os.path.join(ROOT, "best", f"edit4_n{args.n}_size{best_sz}.txt")
        with open(path, "w") as fh:
            fh.write(f"# quaternary single-edit code, n={args.n}, size {best_sz}\n")
            fh.write("# found by symmetry-seeded ILS from the optimal A_4-invariant code\n")
            for w in best_words:
                fh.write(w + "\n")
        print(f"CANDIDATE RECORD saved: {path} -- novelty check required before any claim")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
