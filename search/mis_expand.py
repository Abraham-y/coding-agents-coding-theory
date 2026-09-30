"""Expand a ReduMIS .mis output (one 0/1 per line, vertex order = orbit-graph
vertex order) on an orbit graph into a code, and verify it with the general-d
verifiers.  Saves to best/ only when it beats the given published value.

    python search/mis_expand.py RESULT.mis ORBITS.json n d PUBLISHED [--out FILE]
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from verify_editq import check as verify_editq  # noqa: E402
from verify_editq2 import validate as verify_editq2  # noqa: E402


def main() -> int:
    mis, orbits_path, n, d, published = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
    sel = [i for i, l in enumerate(open(mis)) if l.strip() == "1"]
    orbits = json.load(open(orbits_path))
    words = sorted("".join("0123"[(idx // 4 ** i) % 4] for i in reversed(range(n))) for o in sel for idx in orbits[o])
    ok1, size, why = verify_editq(words, n, d)
    ok2, why2 = verify_editq2(words, n, d, 4)
    print(f"{mis}: {len(sel)} orbits -> {size} words; verify_editq={ok1} {why} verify_editq2={ok2} {why2}; published {published}")
    if out and ok1 and ok2 and size > published:
        body = "".join(x + "\n" for x in words)
        h = hashlib.sha256(body.encode()).hexdigest()
        with open(out, "w") as fh:
            fh.write(f"# quaternary edit-metric code, minimum Levenshtein distance {d}\n# n = {n}\n# d = {d}\n# q = 4\n# size = {size}\n")
            fh.write(f"# found by ReduMIS on the orbit graph ({os.path.basename(mis)}); expanded by search/mis_expand.py\n")
            fh.write(f"# sha256(codewords) = {h}\n{body}")
        print("saved", out, h)
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    raise SystemExit(main())
