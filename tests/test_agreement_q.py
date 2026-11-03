"""verify.py, verify_editq.py and verify_editq2.py must agree on q = 2 and q = 3
edit-metric artifacts, in both directions (CLAUDE.md rules 2 and 3).

verify.py gained q = 2, 3 support on 2026-10-30 (ball marking restricted to the
q-ary alphabet).  The two editq verifiers are pairwise and were q-generic from
the start, so they are the independent check on the new path.

For every q < 4 artifact in best/ and results/fullgraph_q/ this runs:
  accept : the artifact as written
  reject : d raised by one (each artifact's true minimum distance is checked to
           equal its declared d first, so this must fail)
  reject : one codeword moved to distance d-1 of another
  reject : one symbol replaced by q (outside the alphabet)
  reject : a duplicated codeword
plus 300 random small sets with random q, n, d, where the verdict is whatever
the three agree on, and a brute-force distance check decides who is right.

    .venv/bin/python tests/test_agreement_q.py
"""
from __future__ import annotations

import glob
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable


def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def write(path, q, n, d, words):
    with open(path, "w") as fh:
        fh.write(f"# edit-metric test artifact\n# n = {n}\n# d = {d}\n# q = {q}\n")
        fh.write("".join(w + "\n" for w in words))


def verdicts(path):
    out = []
    for script in ("verify.py", "verify_editq.py", "verify_editq2.py"):
        p = subprocess.run([PY, os.path.join(ROOT, script), path], capture_output=True, text=True)
        out.append(p.returncode == 0)
    return out


def read(path):
    n = d = q = None
    words = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            b = line[1:].replace(" ", "")
            for k in ("n", "d", "q"):
                if b.startswith(k + "=") and b[2:].isdigit():
                    if k == "n": n = int(b[2:])
                    if k == "d": d = int(b[2:])
                    if k == "q": q = int(b[2:])
            continue
        words.append(line)
    return q, n, d, words


def close_word(w, q, d, rng):
    """A q-ary word of the same length at edit distance exactly d-1 from w (d >= 2),
    by d-1 substitutions at distinct positions, rechecked by brute force."""
    for _ in range(1000):
        pos = rng.sample(range(len(w)), d - 1)
        x = list(w)
        for p in pos:
            x[p] = rng.choice([c for c in "0123"[:q] if c != x[p]])
        x = "".join(x)
        if lev(w, x) == d - 1:
            return x
    return None


def main():
    rng = random.Random(20261030)
    fails = 0
    n_acc = n_rej = 0
    arts = [p for p in glob.glob(os.path.join(ROOT, "best", "*.txt")) +
            glob.glob(os.path.join(ROOT, "results", "fullgraph_q", "edit*_size*.txt"))]
    with tempfile.TemporaryDirectory() as tmp:
        for path in sorted(arts):
            q, n, d, words = read(path)
            if q not in (2, 3) or d is None:
                continue
            true_d = min(lev(a, b) for i, a in enumerate(words) for b in words[i + 1:]) if len(words) > 1 else 99
            name = os.path.basename(path)
            if true_d != d:
                print(f"NOTE {name}: true minimum distance {true_d} != declared {d}; "
                      f"skipping the d+1 case")
            v = verdicts(path)
            n_acc += 1
            if v != [True] * 3:
                print(f"FAIL accept {name}: {v}"); fails += 1
            cases = []
            if true_d == d:
                cases.append(("d+1", q, n, d + 1, words))
            if d >= 2 and len(words) > 1:
                x = close_word(words[0], q, d, rng)
                if x is not None and x not in words:
                    cases.append(("close pair", q, n, d, [x] + words[1:] + [words[0]]))
            bad = list(words)
            bad[0] = bad[0][:-1] + str(q)
            cases.append(("symbol q", q, n, d, bad))
            cases.append(("duplicate", q, n, d, words + [words[-1]]))
            for label, qq, nn, dd, ws in cases:
                f = os.path.join(tmp, "case.txt")
                write(f, qq, nn, dd, ws)
                v = verdicts(f)
                n_rej += 1
                if v != [False] * 3:
                    print(f"FAIL reject {name} [{label}]: {v}"); fails += 1
        # random sets, judged by brute force
        n_rand = 0
        for t in range(300):
            q = rng.choice((2, 3))
            n = rng.randint(3, 8)
            d = rng.randint(1, n)
            k = rng.randint(2, 12)
            ws = sorted({"".join(rng.choice("0123"[:q]) for _ in range(n)) for _ in range(k)})
            truth = all(lev(a, b) >= d for i, a in enumerate(ws) for b in ws[i + 1:])
            f = os.path.join(tmp, "rand.txt")
            write(f, q, n, d, ws)
            v = verdicts(f)
            n_rand += 1
            if v != [truth] * 3:
                print(f"FAIL random q={q} n={n} d={d} {ws}: verifiers {v}, truth {truth}"); fails += 1
    print(f"{n_acc} accepts, {n_rej} rejects, {n_rand} random sets, {fails} disagreements")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
