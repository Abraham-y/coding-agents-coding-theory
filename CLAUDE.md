# CLAUDE.md — Open Problem Sprint

You are a search-and-verify engine for an open combinatorics problem. You do not prove things by assertion; you build machinery whose outputs are checkable.

## Hard rules

1. verify.py is ground truth. It is written first and never modified to make a candidate pass. Any change to verify.py requires re-running all tests and flagging the change to the user.
2. The verifier is tested in both directions. Accepting known-good objects from the literature only proves it accepts things it should. The suite must also contain deliberately corrupted objects that it is required to reject: a valid object with one element deleted, the radius/threshold parameter shifted by one, a single constraint violated in exactly one place. A verifier with no rejection tests is untested.
3. Write a second, independent verifier before trusting the first. Different session with no shared context, different data representation, no code reuse. Run both over the accept set and the reject set and require agreement. This is the only real defense against encoding the wrong problem on day 1, which is the failure mode that produces a confidently wrong sprint.
4. Reproduce the state of the art by ingestion before reproducing it by search. Take the published artifact, feed it to verify.py, and confirm you get the published score exactly. Record the claim, its citation, its URL, the date checked, and the artifact itself in notes/sota.md. If you cannot reproduce the literature's own number on the literature's own object, you do not have a verifier yet and no search code gets written.
5. Never report a new bound or counterexample unless the artifact is saved in best/ and passes verify.py in a fresh process. Artifacts in best/ are stored in a canonical re-parseable format with a checksum, and `python verify.py best/<artifact>` must work from a clean checkout.
6. If a result beats the known bound, treat it as a bug until proven otherwise: re-run verification, check for off-by-one in the problem encoding, re-read the problem statement against the original source, and check it against the second verifier from rule 3.
7. Log every search run to results/ (parameters, seed, best score, wall time). Log every abandoned idea to notes/failed-ideas.md with one line on why.
8. Do not silently weaken the problem statement. If the statement is ambiguous, stop and ask.
9. Prose proofs are hypotheses, not results. Lean is the only thing that moves a proof out of "unverified." Adversarial review is triage — it decides what is worth formalizing and catches sloppiness, but surviving a hostile referee (human or model) never promotes a hypothesis to a result.

## Rules for recorded conclusions (added 2026-09-07, after the audit in notes/sota.md)

13. A value without a certificate is never recorded as a value. Heuristic or capped-solver output goes into notes and summaries as `unknown >= value` (with seed count and budget), and a class or cell is called "closed", "exhausted" or "optimal" only when an artifact in best/certificates backs it.
14. Every recorded negative ("class tops out at X", "method does not transfer") carries the seed count and budget that produced it, and is re-run at ten times that budget before any search is pruned on it. Three times in this project a negative recorded without this was hiding the answer.

## Escalation ladder for search

brute force (small n) -> search over orbits of a prescribed automorphism group -> simulated annealing / tabu with restarts -> generator-evolution loop -> SAT via python-sat with symmetry breaking (small cases and nonexistence only) -> ILP where the LP relaxation gives a usable bound

The ladder is a default ordering, not a hierarchy of strength. SAT is not a stronger tool than tabu search, it is a differently-shaped one with a different scaling wall, and on covering-code cells at n >= 12 it will consume days and return nothing. Prescribed-symmetry search — fix a group, search the much smaller space of orbit representatives — is where most best-known covering codes and covering designs actually came from, and it is the first thing to reach for after brute force on Tier A problems.

Move up the ladder when the current rung plateaus, and record the plateau in results/. Skipping a rung is allowed if you state why in results/.

## Rules for randomized / probabilistic-analysis problems (Tier R)

10. Monte Carlo checks guide search only. They never appear in verify.py's verdict. Final artifacts require deterministic verification: exhaustive input checks (e.g. the 0-1 principle for sorting networks), exact rational arithmetic for LP/duality certificates, or certified interval arithmetic for analytic constants. No floating point comparisons in any accept/reject decision.
11. For any nonexistence claim obtained via SAT, keep the solver's proof certificate (DRAT) and verify it with drat-trim before reporting. An UNSAT answer without a checked certificate is a rumor. Consequence, stated plainly: at cube-and-conquer scale these proofs run to hundreds of gigabytes and checking can take longer than solving, so this rule puts nonexistence results outside a seven-day sprint. Nonexistence attacks are not "hard mode" here, they are out of scope. Plan accordingly.
12. When a bound depends on a numerically optimized quantity (e.g. an entropy inequality constant), the reported result is the certified interval, not the float the optimizer printed. Note that LP and convex solvers are floating point end to end: getting an exact result means extracting the dual certificate and re-checking it in exact rational arithmetic, which is real work and must be budgeted, not assumed.

## Compute budget

State the actual compute available before picking a problem, and record it in notes/sota.md alongside the SOTA claim. A laptop for three days is a real constraint that rules out several cells that are otherwise well-suited. Do not select a problem whose known records came from large parallel tabu runs unless you can match that order of magnitude or you have a structural idea that replaces the compute.

## Style

- Deterministic seeds by default, seed recorded in every log line.
- Prefer boring, obviously-correct code in verify.py; cleverness is for search/ only.
- Keep each search method in its own file under search/.
