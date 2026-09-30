"""Accept/reject/agreement suite for the GROUND-TRUTH edit-metric verifiers.

Adversarial review found that every existing edit-metric test imported
`search.edit4.verify_edit4` -- the verifier that lives inside the search module,
which CLAUDE.md rules 1 and 3 exclude from ground-truth status. So the paths that
actually decide `python verify.py best/<artifact>` had only their inline selftests
and no cross-agreement test. This file covers verify.verify_edit and
verify2.verify2_edit directly, in both directions, and requires them to agree.

Run:  python tests/test_groundtruth_edit.py
"""
from __future__ import annotations
import os, random, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from verify import verify_edit, edit_ball, checksum          # noqa: E402
from verify2 import verify2_edit, levenshtein                # noqa: E402

FAILURES: list[str] = []
CHECKS = 0
SEED = 20260905

def check(label, got, want):
    global CHECKS
    CHECKS += 1
    if got != want:
        FAILURES.append(f"{label}: got {got!r}, want {want!r}")

# --- accept set: the real record artifacts, through GROUND TRUTH ---------------
RECORDS = ["edit4_n6_C4_product_size116.txt",
           "edit4_n6_C4_diagonal_size116.txt",
           "edit4_n6_C4_diagonal_b_size116.txt"]
loaded = []
for fn in RECORDS:
    p = os.path.join(ROOT, "best", fn)
    if not os.path.exists(p):
        continue
    W = [l.strip() for l in open(p) if l.strip() and not l.startswith("#")]
    loaded.append((fn, W))
    check(f"verify.verify_edit accepts {fn}", verify_edit(W, 6)[:2], (True, 116))
    check(f"verify2.verify2_edit accepts {fn}", verify2_edit(W, 6)[:2], (True, 116))
    # the header's checksum must match what ground truth computes
    hdr = "".join(l for l in open(p) if l.startswith("#"))
    if "sha256(codewords)" in hdr:
        claimed = hdr.split("sha256(codewords) = ")[1].split()[0]
        check(f"checksum matches header for {fn}", checksum(W), claimed)
check("all three records present", len(loaded), 3)

# --- reject set: deliberate corruptions, both verifiers must reject ------------
rng = random.Random(SEED)
if loaded:
    _, base = loaded[0]
    corruptions = {
        "duplicate word": base[:1] + base[:1] + base[2:],
        "wrong length": [base[0][:5]] + base[1:],
        "symbol outside alphabet": ["01234x"] + base[1:],
        # replace a word with one at distance <= 2 from another codeword
        "colliding word": [base[1][:5] + ("0" if base[1][5] != "0" else "1")] + base[1:],
    }
    for label, W in corruptions.items():
        check(f"verify.verify_edit rejects {label}", verify_edit(W, 6)[0], False)
        check(f"verify2.verify2_edit rejects {label}", verify2_edit(W, 6)[0], False)

# --- differential agreement on random inputs ----------------------------------
A = "0123"
for trial in range(300):
    k = rng.randint(2, 5)
    W = ["".join(rng.choice(A) for _ in range(6)) for _ in range(k)]
    a = verify_edit(W, 6)[0]
    b = verify2_edit(W, 6)[0]
    check(f"agreement trial {trial}", a, b)
    # and both must agree with the definition, computed a third way
    ok = len(set(W)) == len(W) and all(
        levenshtein(W[i], W[j]) >= 3 for i in range(len(W)) for j in range(i + 1, len(W)))
    check(f"matches definition trial {trial}", a, ok)

# --- edit_ball is exactly the radius-1 Levenshtein ball ------------------------
for w in ["000", "012", "333"]:
    ball = edit_ball(w)
    for x in ball:
        check(f"ball member {w}->{x}", levenshtein(w, x) <= 1, True)

print(f"checks run : {CHECKS}")
print(f"failures   : {len(FAILURES)}")
for f in FAILURES[:10]:
    print("  ", f)
raise SystemExit(1 if FAILURES else 0)
