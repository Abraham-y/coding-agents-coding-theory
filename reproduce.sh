#!/usr/bin/env bash
# Reproduces the paper's claims: the thirteen improved lower bounds of Table 2
# (E4(6,3) >= 120, E4(7,3) >= 364, E4(6,4) >= 30, E4(7,4) >= 72, E4(8,3) >= 1148,
# E4(9,3) >= 3760, E4(8,4) >= 200, E4(9,4) >= 520, E4(7,5) >= 20, E4(8,5) >= 44,
# E4(9,5) >= 100, E4(8,6) >= 16, E4(9,6) >= 32), the encoding check, and the
# 33-class classification.
#
# Default mode is fast and fully deterministic: it re-verifies the artifacts
# from their definition in a fresh process. The search that FOUND them is
# stochastic; pass --search to re-run it (~10 min/seed, finds 116 in roughly
# 1 seed out of 10).
#
# Usage:  ./reproduce.sh [--search]
set -euo pipefail
cd "$(dirname "$0")"
PY=.venv/bin/python
[ -x "$PY" ] || PY=python3

echo "=== 1. verifier test suites backing this result (accept + reject sets) ==="
# Only the suites this paper's claims depend on, so the check runs with just
# numpy+scipy. ./run_tests.sh covers the whole repo but pulls in optional deps
# (cvxpy, python-sat) for tracks unrelated to these claims.
$PY tests/test_edit4_verify.py 2>&1 | tail -3
$PY tests/test_verify.py       2>&1 | tail -3

echo
echo "=== 2. the record code and the three 116-word codes, from the definition ==="
$PY - <<'PYEOF'
import hashlib
A='0123'
def ball(w):
    B={w}
    for i in range(len(w)):
        B.add(w[:i]+w[i+1:])
        for c in A:
            if c!=w[i]: B.add(w[:i]+c+w[i+1:])
    for i in range(len(w)+1):
        for c in A: B.add(w[:i]+c+w[i:])
    return B
def lev(a,b):
    prev=list(range(len(b)+1))
    for i,ca in enumerate(a,1):
        cur=[i]
        for j,cb in enumerate(b,1):
            cur.append(min(prev[j]+1,cur[j-1]+1,prev[j-1]+(ca!=cb)))
        prev=cur
    return prev[-1]
ok=True; best=0
for f in ['best/edit4_n6_V4_size120.txt','best/edit4_n6_C4_product_size116.txt','best/edit4_n6_C4_diagonal_size116.txt','best/edit4_n6_C4_diagonal_b_size116.txt']:
    W=[l.strip() for l in open(f) if l.strip() and not l.startswith('#')]
    body="".join(w+"\n" for w in sorted(W))
    coll=sum(1 for i in range(len(W)) for j in range(i+1,len(W))
             if ball(W[i])&ball(W[j]))
    mind=min(lev(W[i],W[j]) for i in range(len(W)) for j in range(i+1,len(W)))
    # size comes from the filename, so this cannot silently accept a short code
    want = int(f.rsplit('size',1)[1].split('.')[0])
    good = len(W)==want and len(set(W))==want and coll==0 and mind==3
    ok &= good
    best = max(best, len(W)) if good else best
    print(f"  {f.split('/')[-1]:<36} n={len(W)} collisions={coll} "
          f"minEdit={mind} sha256={hashlib.sha256(body.encode()).hexdigest()[:16]} "
          f"{'PASS' if good else 'FAIL'}")
print(f"  => E4(6,3) >= {best}" if ok else "  => VERIFICATION FAILED")
raise SystemExit(0 if ok else 1)
PYEOF

echo
echo "=== 2b. the 2026-09-07 records at n=7 and d=4, from the definition ==="
$PY tests/test_editq_verify.py 2>&1 | tail -2
$PY - <<'PYEOF2'
import hashlib
def lev(a,b):
    prev=list(range(len(b)+1))
    for i,ca in enumerate(a,1):
        cur=[i]
        for j,cb in enumerate(b,1):
            cur.append(min(prev[j]+1,cur[j-1]+1,prev[j-1]+(ca!=cb)))
        prev=cur
    return prev[-1]
CELLS={'best/edit4_n7_reversal_size359.txt':(7,3,356),
       'best/edit4_n7_C3xC2_size358.txt':(7,3,356),
       'best/edit4_n6_C4_product_size120.txt':(6,3,114),
       'best/edit4_n6_C4_diagonal_size120.txt':(6,3,114),
       'best/editq_n6_d4_size30_class30.txt':(6,4,28),
       'best/editq_n7_d4_size72_class12.txt':(7,4,65),
       'best/editq_n7_d5_size20_class20.txt':(7,5,18),
       'best/edit4_n7_C3xC2_size364.txt':(7,3,356),
       'best/edit4_n8_V4_size1148.txt':(8,3,1132),
       'best/editq_n8_d4_size200_class12.txt':(8,4,176),
       'best/editq_n8_d6_size16_class20.txt':(8,6,14),
       'best/editq_n8_d5_size44_class20.txt':(8,5,38),
       'best/edit4_n9_C4_diagonal_size3760.txt':(9,3,3451),
       'best/editq_n9_d4_size520_class12.txt':(9,4,495),
       'best/editq_n9_d6_size32_class12.txt':(9,6,27),
       'best/editq_n9_d5_size100_class12.txt':(9,5,97)}
ok=True
for f,(n,d,published) in CELLS.items():
    W=[l.strip() for l in open(f) if l.strip() and not l.startswith('#')]
    import re as _re
    want=int(_re.search(r'size(\d+)', f).group(1))
    body="".join(w+"\n" for w in sorted(W))
    mind=min(lev(W[i],W[j]) for i in range(len(W)) for j in range(i+1,len(W)))
    good=(len(W)==want and len(set(W))==want and all(len(w)==n and set(w)<=set('0123') for w in W)
          and mind>=d and want>published)
    ok&=good
    print(f"  {f.split('/')[-1]:<34} E4({n},{d}) size={len(W)} minEdit={mind} "
          f"published={published} sha256={hashlib.sha256(body.encode()).hexdigest()[:16]} "
          f"{'PASS' if good else 'FAIL'}")
print("  => E4(7,3) >= 364, E4(8,3) >= 1148, E4(6,4) >= 30, E4(7,4) >= 72, E4(7,5) >= 20, E4(8,4) >= 200, E4(8,5) >= 44, E4(8,6) >= 16, E4(9,3) >= 3760, E4(9,4) >= 520, E4(9,5) >= 100, E4(9,6) >= 32" if ok else "  => VERIFICATION FAILED")
raise SystemExit(0 if ok else 1)
PYEOF2

echo
echo "=== 3. encoding check: graph vs definition, ALL pairs + exact small cells ==="
$PY search/encoding_check.py

echo "=== 4. arithmetic classification: all 33 classes of S_4 x C_2 ==="
FUNMATH_RESULTS_DIR="${TMPDIR:-/tmp}" $PY search/edit4_diag.py --n 6 --classes 2>&1 | tail -4

if [ "${1:-}" = "--search" ]; then
  echo
  echo "=== 5. re-running the search (stochastic; ~1 seed in 10 reaches 116) ==="
  FUNMATH_RESULTS_DIR="${TMPDIR:-/tmp}" $PY search/edit4_diag.py --n 6 --sweep --seeds 4 --workers 4 --max-seconds 300
fi
echo
echo "reproduce.sh: OK"
