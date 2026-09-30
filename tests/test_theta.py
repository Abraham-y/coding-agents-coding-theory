"""Test the Lovasz theta implementation against graphs with KNOWN theta.

This closes a real gap: theta bounds were being reported before the solver had
ever been checked against a value from the literature. CLAUDE.md rule 2 says test
in both directions; for a numeric tool that means known answers, not just
self-consistency.

Known values used:
    K_n (complete)        theta = 1                    (alpha = 1)
    empty graph on n      theta = n                    (alpha = n)
    C_5 (5-cycle)         theta = sqrt(5) ~ 2.2360680  (Lovasz 1979; alpha = 2)
    C_7 (7-cycle)         theta = 7*cos(pi/7)/(1+cos(pi/7)) ~ 3.3177    (alpha = 3)
    Petersen graph        theta = 4                    (alpha = 4)

Also checked structurally on the real deletion graphs:
    alpha(G) <= theta(G) <= (this repo's LP clique-cover bound)

Run:  python tests/test_theta.py
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from search.theta import certify, edges_of, solve_theta  # noqa: E402

FAIL: list[str] = []
CHECKS = 0
TOL = 1e-4


class SimpleGraph:
    """Minimal stand-in for ConflictGraph: bitset adjacency over 0..n-1."""

    def __init__(self, n: int, edges):
        self.num_vertices = n
        self.adj = [0] * n
        for u, v in edges:
            self.adj[u] |= 1 << v
            self.adj[v] |= 1 << u
        self.words = [str(i) for i in range(n)]


def check_close(label, got, want, tol=TOL):
    global CHECKS
    CHECKS += 1
    if abs(got - want) > tol:
        FAIL.append(f"{label}: got {got:.8f}, want {want:.8f}")


def check(label, cond):
    global CHECKS
    CHECKS += 1
    if not cond:
        FAIL.append(label)


def complete(n):
    return SimpleGraph(n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def empty(n):
    return SimpleGraph(n, [])


def cycle(n):
    return SimpleGraph(n, [(i, (i + 1) % n) for i in range(n)])


def petersen():
    outer = [(i, (i + 1) % 5) for i in range(5)]
    spokes = [(i, i + 5) for i in range(5)]
    inner = [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    return SimpleGraph(10, outer + spokes + inner)


def test_known_values():
    for n in (3, 5):
        g = complete(n)
        t, _, _ = solve_theta(g)
        check_close(f"theta(K_{n}) == 1", t, 1.0)

    for n in (4, 6):
        g = empty(n)
        t, _, _ = solve_theta(g)
        check_close(f"theta(empty_{n}) == {n}", t, float(n))

    g = cycle(5)
    t, _, _ = solve_theta(g)
    check_close("theta(C_5) == sqrt(5)", t, math.sqrt(5))

    g = cycle(7)
    t, _, _ = solve_theta(g)
    c = math.cos(math.pi / 7)
    check_close("theta(C_7) == 7cos(pi/7)/(1+cos(pi/7))", t, 7 * c / (1 + c))

    g = petersen()
    t, _, _ = solve_theta(g)
    check_close("theta(Petersen) == 4", t, 4.0)


def test_certification_is_an_upper_bound():
    """The exact certified value must be >= the true theta, never below."""
    cases = [
        ("C_5", cycle(5), math.sqrt(5), 2),
        ("Petersen", petersen(), 4.0, 4),
        ("K_4", complete(4), 1.0, 1),
        ("empty_5", empty(5), 5.0, 5),
    ]
    for name, g, true_theta, alpha in cases:
        t_float, _, (E, a) = solve_theta(g)
        t_exact = certify(g, t_float, E, a)
        check(f"certified theta({name}) >= true theta", float(t_exact) >= true_theta - 1e-9)
        check(f"floor(certified theta({name})) >= alpha = {alpha}",
              math.floor(t_exact) >= alpha)


def test_sandwich_on_deletion_graphs():
    """alpha <= theta <= LP objective, on the real graphs, at cells we know exactly.

    Two subtleties, both of which produced false failures the first time:

    1. The upper side of the sandwich is the LP's EXACT DUAL OBJECTIVE, not the
       integer bound reported after flooring it. At (s=2,n=7) the objective is
       568903/100000 = 5.689 while the reported bound is 5; theta is 5.242, which
       respects the real sandwich and violates the floored one. Comparing against
       the floor compares incompatible quantities.
    2. The lower side must be checked against the CERTIFIED theta, not the float.
       The float solver returns 3.999997 at (s=2,n=6) where theta is exactly 4, so
       a tight tolerance on the float fails for accuracy reasons alone. The
       certified value is the one required to dominate alpha, and it does.
    """
    import json
    from fractions import Fraction

    from search.graph import load_graph

    known = {(6, 2): 4, (7, 2): 5, (7, 1): 16}
    for (n, s), alpha in known.items():
        g = load_graph(n, s)
        t_float, _, (E, a) = solve_theta(g)

        with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                               "results", f"cert_n{n}_s{s}.json")) as fh:
            lp_obj = float(Fraction(*json.load(fh)["bound"]))

        # float side: only a loose sanity check, at solver accuracy
        check(f"theta float >= alpha - tol at (s={s},n={n})", t_float >= alpha - 1e-3)
        check(f"theta <= LP exact objective at (s={s},n={n})", t_float <= lp_obj + 1e-6)

        # the check that actually matters: the CERTIFIED value is a valid bound
        t_exact = certify(g, t_float, E, a)
        check(f"certified theta >= alpha at (s={s},n={n})", float(t_exact) >= alpha)
        check(f"floor(certified theta) >= alpha at (s={s},n={n})",
              math.floor(t_exact) >= alpha)


def main() -> int:
    test_known_values()
    test_certification_is_an_upper_bound()
    test_sandwich_on_deletion_graphs()
    print(f"checks run : {CHECKS}")
    print(f"failures   : {len(FAIL)}")
    for f in FAIL:
        print("  FAIL " + f)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
