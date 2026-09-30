"""Varshamov-Tenengolts code construction. TEST FIXTURE ONLY.

verify.py must never import this. This exists so the test suite can build
known-good objects from the literature independently of any search code.

    VT_a(n) = { v in {0,1}^n : sum_{i=1}^{n} i * v_i  ==  a  (mod n+1) }

Levenshtein proved every VT_a(n) is single-deletion-correcting.
"""

from __future__ import annotations


def vt_code(n: int, a: int = 0) -> list[str]:
    """All length-n binary words whose position-weighted sum is a mod (n+1)."""
    modulus = n + 1
    out = []
    for x in range(1 << n):
        # bit i-1 of x is the symbol at 1-indexed position i
        total = 0
        for i in range(1, n + 1):
            if (x >> (i - 1)) & 1:
                total += i
        if total % modulus == a % modulus:
            out.append(format(x, f"0{n}b")[::-1])
    return out


# Published |VT_0(n)|, transcribed from Table 1 ("VT_0(n) code" row) of
# Weindel & Heckel, "LLM-Guided Search for Deletion-Correcting Codes",
# arXiv:2504.00613v2 (5 Jun 2026). Checked 2026-08-03.
PUBLISHED_VT0_SIZES = {
    6: 10,
    7: 16,
    8: 30,
    9: 52,
    10: 94,
    11: 172,
    12: 316,
    13: 586,
    14: 1096,
    15: 2048,
    16: 3856,
    17: 7286,
    18: 13798,
    19: 26216,
    20: 49940,
}

# Exact optimal sizes A_s(n) that the literature reports as SETTLED, i.e. where a
# published upper bound meets a published construction. Same source, Tables 1-2.
#   s = 1: settled for n <= 11 (VT optimality proven only to n = 11)
#   s = 2: settled for n <= 11 (ILP upper bound meets best known construction)
#
# CORRECTION, 2026-08-03. A_2(8) was transcribed here as 8, read off the
# "ILP lower bound" row of Table 2. That is wrong: **A_2(8) = 7**, established
# here by two independent exact methods that share no code --
#   (a) set-packing MILP solved to proven optimality by HiGHS, and
#   (b) exact maximum clique on the complement graph via networkx --
# and corroborated by the paper's own repository
# (github.com/MLI-lab/DistributedFunSearch, analysis/utils/best_known.py lists
# (8,2): 7) and by both MIS-solver rows of Table 2, which also give 7.
# A *lower* bound of 8 is impossible when the maximum is 7, so either the table
# entry is a typo or it was misread here; either way 7 is the verified value.
SETTLED_OPTIMA = {
    1: {6: 10, 7: 16, 8: 30, 9: 52, 10: 94, 11: 172},
    2: {6: 4, 7: 5, 8: 7, 9: 11, 10: 16, 11: 24},
    3: {6: 2, 7: 2, 8: 4, 9: 5, 10: 6, 11: 8},
}

# Cells this repo has closed EXACTLY, by two independent exact methods each.
# These are proofs of optimality, not best-known constructions.
EXACTLY_SETTLED_HERE = {
    (2, 6): 4,
    (2, 7): 5,
    (2, 8): 7,
    (2, 9): 11,
}

# Open cells: (s, n) -> (best known construction, best published upper bound)
OPEN_CELLS = {
    (1, 12): (316, 320),
    (1, 13): (586, 592),
    (1, 14): (1096, 1104),
    (1, 15): (2048, 2251),
    (1, 16): (3856, 4202),
    (2, 12): (37, 53),
    (2, 13): (56, 93),
    (2, 14): (87, 160),
    (2, 15): (136, 277),
    (2, 16): (215, 487),
    (3, 12): (12, 16),
    (3, 13): (16, 24),
}
