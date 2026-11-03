"""Two-directional tests for verify_big.py (CLAUDE.md rule 2).

ACCEPT: every edit-metric artifact in best/ (they pass verify.py / verify_editq.py).
REJECT: corrupted copies of accepted codes, one defect each:
  near    one codeword replaced by a word within distance d-1 of another
  dup     a codeword repeated
  len     a codeword truncated by one symbol
  sym     a symbol outside 0-3
  dplus1  a code whose true minimum distance is exactly d, checked at d+1
  nod     '# d' removed and no --d given (must exit 2, not guess)
  sha     header sha256 that does not match the codewords
RANDOM: verify_big.check() against a plain Wagner-Fischer reference written here,
on random small sets (fixed seed), both valid and invalid.

    .venv/bin/python tests/test_verify_big.py      exit 0 iff every expectation holds
"""
import glob, os, random, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import verify_big  # noqa: E402

PY = sys.executable
fails, counts = [], {"accept": 0, "reject": 0, "random": 0}


def ref_lev(a, b):
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def ref_min_dist(words):
    return min(ref_lev(a, b) for i, a in enumerate(words) for b in words[i + 1:])


def run(path, extra=()):
    return subprocess.run([PY, os.path.join(ROOT, "verify_big.py"), path, *extra],
                          capture_output=True, text=True).returncode


def read(path):
    lines = open(path).read().splitlines()
    head = [l for l in lines if l.startswith("#")]
    words = [l.strip() for l in lines if l.strip() and not l.startswith("#")]
    m = re.search(r"#\s*n\s*=\s*(\d+)", "\n".join(head))
    if m:
        n = int(m.group(1))
    else:  # two old artifacts state n only in prose ("n=6"); flagged by the caller
        n = -int(re.search(r"\bn=(\d+)", "\n".join(head)).group(1))
    m = re.search(r"#\s*d\s*=\s*(\d+)", "\n".join(head))
    return head, words, n, (int(m.group(1)) if m else None)


def write(head, words, tmp, name):
    p = os.path.join(tmp, name)
    open(p, "w").write("\n".join(head + words) + "\n")
    return p


def strip_hdr(head):
    """Drop size and sha lines so a deliberately altered body is judged on the codes."""
    return [h for h in head if not re.search(r"size\s*=|sha256", h)]


def expect(cond, what):
    if not cond:
        fails.append(what)


files = sorted(glob.glob(os.path.join(ROOT, "best", "edit4_*.txt")) +
               glob.glob(os.path.join(ROOT, "best", "editq_*.txt")))
with tempfile.TemporaryDirectory() as tmp:
    for f in files:
        head, words, n, d = read(f)
        if n < 0:  # no machine-readable '# n': must be refused, then accepted once headed
            rc = run(f, ("--d", "3")); counts["reject"] += 1
            expect(rc == 2, f"REJECT missing-n failed (rc={rc}): {os.path.basename(f)}")
            n = -n
            head = head + [f"# n = {n}"]
            f = write(head, words, tempfile.gettempdir(), "headed_" + os.path.basename(f))
        extra = () if d is not None else ("--d", "3")
        d = d if d is not None else 3
        rc = run(f, extra)
        counts["accept"] += 1
        expect(rc == 0, f"ACCEPT failed (rc={rc}): {os.path.basename(f)}")
        if len(words) > 600:
            continue
        h = [f"# n = {n}", f"# d = {d}"]   # minimal explicit header for corrupted copies
        base = os.path.basename(f)[:-4]
        # near: replace words[1] by words[0] with d-1 substitutions
        a = words[0]
        b = "".join(str((int(c) + 1) % 4) if i < d - 1 else c for i, c in enumerate(a))
        if b not in words:
            bad = [a, b] + words[2:]
            expect(ref_lev(a, b) <= d - 1, "test bug: near word not near")
            rc = run(write(h, bad, tmp, base + "_near.txt")); counts["reject"] += 1
            expect(rc == 1, f"REJECT near failed (rc={rc}): {base}")
        rc = run(write(h, words + [words[0]], tmp, base + "_dup.txt")); counts["reject"] += 1
        expect(rc == 1, f"REJECT dup failed (rc={rc}): {base}")
        rc = run(write(h, [words[0][:-1]] + words[1:], tmp, base + "_len.txt")); counts["reject"] += 1
        expect(rc == 1, f"REJECT len failed (rc={rc}): {base}")
        rc = run(write(h, ["4" + words[0][1:]] + words[1:], tmp, base + "_sym.txt")); counts["reject"] += 1
        expect(rc == 1, f"REJECT sym failed (rc={rc}): {base}")
        h_nod = [x for x in h if not re.search(r"#\s*d\s*=", x)]
        rc = run(write(h_nod, words, tmp, base + "_nod.txt")); counts["reject"] += 1
        expect(rc == 2, f"REJECT nod failed (rc={rc}): {base}")
        rc = run(write(h + ["# sha256(codewords) = " + "0" * 64], words, tmp, base + "_sha.txt"))
        counts["reject"] += 1
        expect(rc == 1, f"REJECT sha failed (rc={rc}): {base}")
        if len(words) <= 130:
            true_d = ref_min_dist(words)
            expect(true_d >= d, f"test bug: {base} true min distance {true_d} < {d}")
            if true_d == d:
                h_up = [x for x in h if not re.search(r"#\s*d\s*=", x)] + [f"# d = {d + 1}"]
                rc = run(write(h_up, words, tmp, base + "_dplus1.txt")); counts["reject"] += 1
                expect(rc == 1, f"REJECT dplus1 failed (rc={rc}): {base}")

rng = random.Random(20261009)
for trial in range(400):
    n = rng.choice([4, 5, 6])
    d = rng.choice([1, 2, 3, 4, 5])
    pool = ["".join(rng.choice("0123") for _ in range(n)) for _ in range(200)]
    pool = list(dict.fromkeys(pool))
    if trial % 2:            # build a valid code greedily with the reference DP
        code = []
        for w in pool:
            if all(ref_lev(w, c) >= d for c in code):
                code.append(w)
            if len(code) >= 25:
                break
    else:                    # random set, usually invalid
        code = pool[: rng.randint(2, 25)]
    truth = len(code) < 2 or ref_min_dist(code) >= d
    got, _ = verify_big.check(code, n, d)
    counts["random"] += 1
    expect(got == truth, f"RANDOM disagree n={n} d={d} size={len(code)} truth={truth} got={got}")

print(f"accept checked: {counts['accept']}, reject checked: {counts['reject']}, "
      f"random agreement trials: {counts['random']}, failures: {len(fails)}")
for f in fails[:20]:
    print("  FAIL", f)
sys.exit(1 if fails else 0)
