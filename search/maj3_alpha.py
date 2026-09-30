"""Exact computation of the MNSX constant alpha_k for recursive 3-majority.

Background
----------
Magniez, Nayak, Santha, Sherman, Tardos and Xiao (arXiv:1309.7565) prove

    R_delta(3-MAJ_h)  >=  (1 - 2*delta) * (alpha_k / 2^k) * (2 + alpha_k^(-1/k))^h

so the base of the exponential is 2 + alpha_k^(-1/k), and a SMALLER alpha_k^(1/k)
gives a BETTER lower bound. They computed

    alpha_1 = 2,  alpha_2 = 24/7,  alpha_3 = 12231/2203,  alpha_4 = 2027349/216164

giving bases 2.5, 2.54006, 2.56474 and 2.57143. The last is the best known lower
bound on the randomized query complexity of recursive majority, and it has stood
since 2013. They stopped at k=4 because their DP is over "stable configurations"
whose count obeys N_k = C(N_{k-1}+1,2) + C(N_{k-1}+2,3), giving
N_4 = 246,792 but N_5 = 2,505,258,478,767,772, and they wrote that "a much more
efficient algorithm would be needed to make the calculation for k = 5 feasible."

Definitions implemented here
---------------------------
H_k      hard inputs of height k: recursively, at every internal node the three
         children values are not all equal, so exactly one child is a minority.
m(x)     absolute minority: follow the minority child from the root to a leaf.
sensitive(x)  the 2^k leaves whose flip flips 3-MAJ_k(x).

For a deterministic decision tree C on 3^k variables querying at least one bit,

    rho_alpha(C) = 2^(-k) * E[# sensitive bits queried] - alpha * Pr[m(x) queried]

with both taken over x uniform on the 0-hard inputs H_k^0, and

    alpha_k = max_C alpha_C   where   alpha_C  is the root of rho_alpha(C) = 0.

Equivalently alpha_k > alpha iff max_C rho_alpha(C) > 0.

Method
------
Incremental DP over configurations. Querying variable v at configuration gamma
contributes 2^(-k)*Pr[v sensitive | gamma] - alpha*Pr[v = m(x) | gamma] and then
splits into the two child configurations, so

    V(gamma) = max( 0, max_v [ contribution(v) + sum_b Pr[x_v=b|gamma] V(gamma+{v=b}) ] )

with the stop option unavailable at the root. Everything is fractions.Fraction:
no float is compared anywhere (CLAUDE.md rule 10), and the reported constant is
the exact rational, not a decimal (rule 12).

alpha_k is found by the same fixed-point iteration the authors used: run the DP
at a trial alpha, read off the maximiser's (P_q, P_m), set alpha <- 2^(-k) P_q/P_m,
repeat until stable.

This file is SEARCH code, not ground truth. Its outputs are checked against the
published alpha_1..alpha_4 in tests/test_maj3.py.

Run:
    python search/maj3_alpha.py --k 1
    python search/maj3_alpha.py --k 2
    python search/maj3_alpha.py --k 3 --verbose
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from functools import lru_cache

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(ROOT, "results")

PUBLISHED_ALPHA = {
    1: Fraction(2),
    2: Fraction(24, 7),
    3: Fraction(12231, 2203),
    4: Fraction(2027349, 216164),
}

# N_k = number of equivalence classes of stable configurations (MNSX section 3.2)
PUBLISHED_N = {0: 1, 1: 2, 2: 7, 3: 112, 4: 246792, 5: 2505258478767772}


def stable_config_counts(kmax: int) -> dict[int, int]:
    """N_0 = 1, N_k = C(N_{k-1}+1, 2) + C(N_{k-1}+2, 3)."""
    out = {0: 1}
    n = 1
    for k in range(1, kmax + 1):
        n = (n + 1) * n // 2 + (n + 2) * (n + 1) * n // 6
        out[k] = n
    return out


# ---------------------------------------------------------------------------
# hard inputs
# ---------------------------------------------------------------------------


def maj3(a: int, b: int, c: int) -> int:
    return 1 if (a + b + c) >= 2 else 0


def build_hard_inputs(k: int):
    """Return (inputs, values, minority_pos, sensitive_sets) for height k.

    inputs[i] is a tuple of 3^k bits. Only HARD inputs are generated: at every
    internal node the three child values are not all equal.
    """
    if k == 0:
        inputs = [(0,), (1,)]
        values = [0, 1]
        minority = [0, 0]          # the single leaf is trivially the minority path end
        sensitive = [frozenset({0}), frozenset({0})]
        return inputs, values, minority, sensitive

    sub_inputs, sub_values, sub_minority, sub_sensitive = build_hard_inputs(k - 1)
    width = 3 ** (k - 1)

    inputs, values, minority, sensitive = [], [], [], []
    for i, xi in enumerate(sub_inputs):
        for j, xj in enumerate(sub_inputs):
            for l, xl in enumerate(sub_inputs):
                vi, vj, vl = sub_values[i], sub_values[j], sub_values[l]
                if vi == vj == vl:
                    continue  # not a hard input
                root_val = maj3(vi, vj, vl)
                # exactly one child is the minority
                child_vals = (vi, vj, vl)
                minority_child = child_vals.index(1 - root_val)

                word = xi + xj + xl
                inputs.append(word)
                values.append(root_val)

                # absolute minority: descend into the minority child
                sub_idx = (i, j, l)[minority_child]
                minority.append(minority_child * width + sub_minority[sub_idx])

                # a leaf is sensitive iff every node on its root path has the
                # same value as the root, i.e. it lies in a MAJORITY child and is
                # sensitive within that child.
                sens = set()
                for pos, (child_val, child_idx) in enumerate(
                    zip(child_vals, (i, j, l))
                ):
                    if child_val == root_val:
                        for leaf in sub_sensitive[child_idx]:
                            sens.add(pos * width + leaf)
                sensitive.append(frozenset(sens))

    return inputs, values, minority, sensitive


# ---------------------------------------------------------------------------
# DP
# ---------------------------------------------------------------------------


class AlphaSolver:
    def __init__(self, k: int):
        self.k = k
        self.nvars = 3 ** k
        inputs, values, minority, sensitive = build_hard_inputs(k)
        # restrict to the 0-hard inputs
        keep = [i for i, v in enumerate(values) if v == 0]
        self.inputs = [inputs[i] for i in keep]
        self.minority = [minority[i] for i in keep]
        self.sensitive = [sensitive[i] for i in keep]
        self.num_hard = len(self.inputs)
        self.weight = Fraction(1, 2) ** k  # 2^{-k}
        self.states_seen = 0

    # -- the recursion -----------------------------------------------------
    def solve(self, alpha: Fraction):
        """Return (value, Pq, Pm) for the optimal tree at this alpha."""
        self.states_seen = 0
        all_inputs = frozenset(range(self.num_hard))
        return self._value(all_inputs, frozenset(), alpha, root=True)

    def _value(self, S: frozenset, Q: frozenset, alpha: Fraction, root: bool = False):
        return self._memo_value(S, Q, alpha, root)

    def _memo_value(self, S, Q, alpha, root):
        key = (S, Q, root)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        self.states_seen += 1

        best = None if root else (Fraction(0), Fraction(0), Fraction(0))
        size = len(S)
        for v in range(self.nvars):
            if v in Q:
                continue
            zeros = [i for i in S if self.inputs[i][v] == 0]
            ones = [i for i in S if self.inputs[i][v] == 1]
            n_sens = sum(1 for i in S if v in self.sensitive[i])
            n_min = sum(1 for i in S if self.minority[i] == v)

            contrib = self.weight * Fraction(n_sens, size) - alpha * Fraction(n_min, size)
            pq = Fraction(n_sens, size)
            pm = Fraction(n_min, size)

            Q2 = Q | {v}
            total = contrib
            for branch in (zeros, ones):
                if not branch:
                    continue
                p = Fraction(len(branch), size)
                sub_val, sub_pq, sub_pm = self._memo_value(
                    frozenset(branch), Q2, alpha, False
                )
                total += p * sub_val
                pq += p * sub_pq
                pm += p * sub_pm

            if best is None or total > best[0]:
                best = (total, pq, pm)

        if best is None:
            best = (Fraction(0), Fraction(0), Fraction(0))
        self._cache[key] = best
        return best

    # -- fixed point on alpha ---------------------------------------------
    def find_alpha(self, alpha0: Fraction = Fraction(0), max_iter: int = 40, verbose=False):
        alpha = alpha0
        history = []
        for it in range(max_iter):
            self._cache = {}
            value, pq, pm = self.solve(alpha)
            if pm == 0:
                raise RuntimeError("optimal tree queries the minority with probability 0")
            new_alpha = self.weight * pq / pm
            history.append(
                {"iter": it, "alpha_in": str(alpha), "rho": str(value),
                 "Pq": str(pq), "Pm": str(pm), "alpha_out": str(new_alpha),
                 "states": self.states_seen}
            )
            if verbose:
                print(f"  iter {it}: alpha={alpha} -> rho={value} "
                      f"Pq={pq} Pm={pm} -> alpha={new_alpha} ({self.states_seen} states)")
            if new_alpha == alpha:
                return alpha, history
            alpha = new_alpha
        raise RuntimeError("alpha iteration did not converge")


# ---------------------------------------------------------------------------


def base_from_alpha(alpha: Fraction) -> float:
    """The lower-bound base 2 + alpha^(-1/k) is irrational; float is DISPLAY only."""
    return float(alpha)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print(f"stable-configuration counts N_k: {stable_config_counts(6)}")

    t0 = time.time()
    solver = AlphaSolver(args.k)
    print(f"k={args.k}: {solver.nvars} variables, {solver.num_hard} zero-hard inputs")

    alpha, history = solver.find_alpha(verbose=args.verbose)
    elapsed = time.time() - t0

    published = PUBLISHED_ALPHA.get(args.k)
    match = (published is not None and alpha == published)
    print(f"alpha_{args.k} = {alpha}  (= {float(alpha):.10f})")
    if published is not None:
        print(f"published  = {published}  -> {'MATCH' if match else 'MISMATCH'}")
    root = float(alpha) ** (-1.0 / args.k)
    print(f"lower-bound base 2 + alpha^(-1/k) = {2 + root:.8f}  [float, display only]")
    print(f"elapsed {elapsed:.1f}s, {solver.states_seen} DP states at final iteration")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    path = os.path.join(RESULTS_DIR, f"maj3_alpha_k{args.k}.json")
    with open(path, "w") as fh:
        json.dump(
            {
                "k": args.k,
                "alpha_numerator": alpha.numerator,
                "alpha_denominator": alpha.denominator,
                "alpha_str": str(alpha),
                "published": str(published) if published else None,
                "matches_published": match,
                "base_float_display_only": 2 + root,
                "num_hard_inputs": solver.num_hard,
                "elapsed_seconds": round(elapsed, 2),
                "iterations": history,
            },
            fh,
            indent=2,
        )
    print(f"log: {path}")
    return 0 if (published is None or match) else 1


if __name__ == "__main__":
    raise SystemExit(main())
