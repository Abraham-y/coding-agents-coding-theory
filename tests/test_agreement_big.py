"""Agreement test for the edit-metric verifiers (CLAUDE.md rule 3).

verify.py (ground truth; ball-ownership map, general d since 2026-10-10),
verify_big.py (ball-collision sort, numpy) and verify_big2 (pairwise Myers
bit-parallel edit distance, C; written in a separate session from the definition
only) must return the same verdict on:
  ACCEPT  every edit-metric artifact in best/ (headed copies for the two legacy
          files that state n only in prose)
  REJECT  corrupted copies of the small ones: near / dup / len / sym / dplus1
  RANDOM  small random sets, half valid by construction (fixed seed)
Exit codes must agree exactly (0 valid, 1 invalid).

    .venv/bin/python tests/test_agreement_big.py     (build verify_big2 first)
"""
import glob, os, random, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
G = [PY, os.path.join(ROOT, "verify.py")]
A = [PY, os.path.join(ROOT, "verify_big.py")]
B = [os.path.join(ROOT, "verify_big2")]
fails, n_cases = [], {"accept": 0, "reject": 0, "random": 0}


def rc(cmd, path):
    return subprocess.run(cmd + [path], capture_output=True, text=True).returncode


def ref_lev(a, b):
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def both(path, want, label, words=None):
    ra, rb = rc(A, path), rc(B, path)
    # verify.py tells binary from quaternary by alphabet, so a set using only 0/1
    # is (correctly) refused as an undeclared binary artifact; skip it there.
    binary_only = words is not None and set("".join(words)) <= set("01")
    rg = want if binary_only else (0 if rc(G, path) == 0 else 1)
    if not (ra == rb == want and rg == want):
        fails.append(f"{label}: verify.py={rg} verify_big={ra} verify_big2={rb} expected={want}")


def write(tmp, name, n, d, words):
    p = os.path.join(tmp, name)
    open(p, "w").write(f"# quaternary edit-metric code\n# n = {n}\n# d = {d}\n" + "\n".join(words) + "\n")
    return p


with tempfile.TemporaryDirectory() as tmp:
    for f in sorted(glob.glob(os.path.join(ROOT, "best", "edit4_*.txt")) +
                    glob.glob(os.path.join(ROOT, "best", "editq_*.txt"))):
        lines = open(f).read().splitlines()
        head = "\n".join(l for l in lines if l.startswith("#"))
        words = [l.strip() for l in lines if l.strip() and not l.startswith("#")]
        n = len(words[0])
        m = re.search(r"^#\s*d\s*=\s*(\d+)", head, re.M)
        d = int(m.group(1)) if m else 3
        base = os.path.basename(f)[:-4]
        both(write(tmp, base + ".txt", n, d, words), 0, f"ACCEPT {base}", words)
        n_cases["accept"] += 1
        if len(words) > 600:
            continue
        a = words[0]
        near = "".join(str((int(c) + 1) % 4) if i < d - 1 else c for i, c in enumerate(a))
        variants = {"dup": words + [a], "len": [a[:-1]] + words[1:],
                    "sym": ["4" + a[1:]] + words[1:]}
        if near not in words:
            variants["near"] = [a, near] + words[2:]
        for k, ws in variants.items():
            both(write(tmp, f"{base}_{k}.txt", n, d, ws), 1, f"REJECT {k} {base}", ws)
            n_cases["reject"] += 1
        if len(words) <= 130:
            dmin = min(ref_lev(x, y) for i, x in enumerate(words) for y in words[i + 1:])
            if dmin == d:
                both(write(tmp, f"{base}_dplus1.txt", n, d + 1, words), 1, f"REJECT dplus1 {base}", words)
                n_cases["reject"] += 1

    rng = random.Random(20261010)
    for t in range(200):
        n, d = rng.choice([4, 5, 6, 7]), rng.choice([2, 3, 4, 5])
        pool = list(dict.fromkeys("".join(rng.choice("0123") for _ in range(n)) for _ in range(300)))
        if t % 2:
            code = []
            for w in pool:
                if all(ref_lev(w, c) >= d for c in code):
                    code.append(w)
                if len(code) >= 30:
                    break
        else:
            code = pool[: rng.randint(2, 30)]
        truth = 0 if all(ref_lev(x, y) >= d for i, x in enumerate(code) for y in code[i + 1:]) else 1
        both(write(tmp, f"rand{t}.txt", n, d, code), truth, f"RANDOM {t} n={n} d={d}", code)
        n_cases["random"] += 1

print(f"accept {n_cases['accept']}, reject {n_cases['reject']}, random {n_cases['random']}; "
      f"disagreements or wrong verdicts: {len(fails)}")
for f in fails[:20]:
    print("  FAIL", f)
sys.exit(1 if fails else 0)
