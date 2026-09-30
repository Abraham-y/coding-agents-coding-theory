"""Heavy parallel ILS on the raw quaternary single-edit graph. Target: beat 113 at n=6.

A proper script (not python -c) because multiprocessing under macOS spawn mode
must be able to import the worker function.

Run:  python search/edit4_ils.py --n 6 --seconds 1200 --workers 6 --seed 10
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.edit4 import PUBLISHED, Edit4Graph, verify_edit4  # noqa: E402
from search.ils import Searcher  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def work(args):
    n, seed, seconds = args
    g = Edit4Graph(n)
    s = Searcher(g, seed)
    best, size, iters, el = s.run(seconds=seconds)
    return seed, sorted(g.words[v] for v in best), iters


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--seconds", type=float, default=1200)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--seed", type=int, default=10)
    args = ap.parse_args()

    t0 = time.time()
    jobs = [(args.n, args.seed + i, args.seconds) for i in range(args.workers)]
    with Pool(args.workers) as pool:
        results = pool.map(work, jobs)

    best_sz, best_words, best_seed = 0, None, None
    for seed, words, iters in sorted(results, key=lambda r: -len(r[1])):
        ok, sz, reason = verify_edit4(words, args.n)
        print(f"seed {seed}: size={len(words)} iters={iters} verify={ok}"
              + ("" if ok else f"  REJECTED {reason}"))
        if ok and sz > best_sz:
            best_sz, best_words, best_seed = sz, words, seed

    pub = PUBLISHED.get(args.n, {}).get("best")
    print(f"BEST RAW: {best_sz}  (published best known: {pub})  wall={time.time()-t0:.0f}s")

    stamp = time.strftime("%Y%m%dT%H%M%S")
    with open(os.path.join(ROOT, "results", f"edit4_ils_n{args.n}_{stamp}.json"), "w") as fh:
        json.dump({"n": args.n, "seconds": args.seconds, "workers": args.workers,
                   "base_seed": args.seed,
                   "sizes": {str(s): len(w) for s, w, _ in results},
                   "best": best_sz, "published_best": pub}, fh, indent=2)

    if best_words and pub and best_sz > pub:
        path = os.path.join(ROOT, "best", f"edit4_n{args.n}_size{best_sz}.txt")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(f"# quaternary single-edit code, n={args.n}, size {best_sz}, "
                     f"raw ILS seed {best_seed}\n")
            for w in best_words:
                fh.write(w + "\n")
        print(f"CANDIDATE RECORD saved: {path} -- novelty check required before any claim")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
