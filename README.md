# Supplementary material for "Coding Agents for Coding Theory"

Companion code and data for the NeurIPS 2026 MATH-AI workshop paper, [arXiv:2609.39081](https://arxiv.org/abs/2609.39081). MIT licensed (see LICENSE).

Everything needed to check the paper's claims. No network access required.

    pip install numpy scipy      # the paper's claims need only these
    ./reproduce.sh               # a few minutes

`reproduce.sh` runs, in order: the verifier test suites (accept AND reject
sets); every improved lower bound of Table 2 re-verified from the definition in
a fresh process, with no imports from this code base; the encoding check
(confusability-graph adjacency against the Levenshtein definition on all
8,386,560 vertex pairs at n=6, plus the exact small cells); and the arithmetic
classification of all 33 conjugacy classes of subgroups of S_4 x C_2.

## The claims

The thirteen codes of Table 2, listed in full in Appendix C of the paper, are in
best/ (for example best/edit4_n6_V4_size120.txt for E4(6,3) >= 120). The sha256
of each file's sorted codewords, one per line, matches the value printed in the
paper. The four solver-reported class optima are best/*_exact.txt; none has a
checked certificate.

## Layout

    verify.py, verify2.py   ground-truth verifiers (binary deletion codes; verify.py
                            also checks the edit-metric artifacts)
    verify_editq.py, verify_editq2.py
                            independent verifiers for q-ary edit-metric codes, any d
    tests/                  two-directional suites: objects they must ACCEPT,
                            deliberately corrupted objects they must REJECT
    search/edit4.py         confusability graph and the original sweep (contains the
                            --exact-cap defect described in Section 5)
    search/edit4_diag.py    the 33-class S_4 x C_2 classification, incl. diagonal classes
    search/edit4_orbit_ls.py  our Python orbit-graph local search
    search/mis_ils.c        C reimplementation of the ARW local search
    search/orbit_graph.c    C orbit-graph builder used for n >= 8
    search/encoding_check.py  the all-pairs encoding check
    best/                   every code artifact, plus certificate hashes
    results/                run records behind the paper's numbers, including
                            results/redumis/ (the ReduMIS controls) and
                            edit4_sweep_n6_20260806T025537.json -- the run whose
                            MILP reported 100, 92 and 8 for three isomorphic problems
    notes/sota.md           dated project log (the source of Figure 1)
    notes/failed-ideas.md   log of abandoned approaches with causes of death
    CLAUDE.md               the protocol document
    lean/                   the orbit-reduction lemma in Lean 4 + mathlib

## Full test suite

`./run_tests.sh` covers the whole repository, including tracks unrelated to this
paper (recursive-majority constants, Lovasz theta bounds). Those extra suites
need the optional dependencies in requirements.txt; without them they fail on
import, which is expected and does not affect any claim in the paper.

## After the paper (updated 2026-11-03)

The current bounds for all five of Houghten's edit-metric tables, with a code
file for every improved entry, are kept at
[abrahamyeung.com/edit-codes](https://abrahamyeung.com/edit-codes/).

`verify.py` has changed three times since the paper. Each change was followed by
a full re-run of the test suite.

- 2026-10-10: it now enforces the `# d` declared in an edit-metric artifact.
  The paper-era version checked d = 3 whatever the header said, so the d >= 4
  codes of Table 2 were checked only by `verify_editq.py`, `verify_editq2.py`
  and `reproduce.sh`. All of them pass the new check.
- 2026-10-12: a bitset path, so a code of 500,882 words verifies in about a minute.
- 2026-10-30: binary and ternary codes, declared with `# q = 2` or `# q = 3`.

`verify_big.py` and `verify_big2.c` are two further verifiers, written
independently, for codes too large for the pairwise checkers.
`tests/test_agreement_big.py` (quaternary) and `tests/test_agreement_q.py`
(binary and ternary) run all of the verifiers over each code in best/, over
deliberately corrupted copies of them, and over random sets, and require them
to agree.

New codes in best/, each accepted by `verify.py` in a fresh process and by the
independent verifiers. The previous value is the lower bound in Houghten's table.

| cell | new | previous | file |
|---|---:|---:|---|
| E2(11,3) | 106 | 104 | best/edit2_n11_d3_size106.txt |
| E2(11,4) | 33 | 32 | best/edit2_n11_d4_size33.txt |
| E2(12,3) | 180 | 179 | best/edit2_n12_d3_size180.txt |
| E2(12,4) | 55 | 54 | best/edit2_n12_d4_size55.txt |
| E3(7,3) | 83 | 71 | best/edit3_n7_d3_size83.txt |
| E3(7,4) | 22 | 21 | best/edit3_n7_d4_size22.txt |
| E3(8,3) | 201 | 179 | best/edit3_n8_d3_size201.txt |
| E3(8,4) | 45 | 41 | best/edit3_n8_d4_size45.txt |
| E3(9,3) | 503 | 438 | best/edit3_n9_d3_size503.txt |
| E3(9,4) | 97 | 88 | best/edit3_n9_d4_size97.txt |
| E3(9,5) | 26 | 25 | best/edit3_n9_d5_size26.txt |
| E3(10,3) | 1261 | 1133 | best/edit3_n10_d3_size1261.txt |
| E3(10,4) | 214 | 202 | best/edit3_n10_d4_size214.txt |
| E3(10,5) | 54 | 49 | best/edit3_n10_d5_size54.txt |
| E3(10,6) | 19 | 18 | best/edit3_n10_d6_size19.txt |
| E3(11,3) | 3205 | 2972 | best/edit3_n11_d3_size3205.txt |
| E3(11,4) | 525 | 471 | best/edit3_n11_d4_size525.txt |
| E3(11,5) | 108 | 98 | best/edit3_n11_d5_size108.txt |
| E3(11,6) | 36 | 29 | best/edit3_n11_d6_size36.txt |
| E4(12,3) | 141943 | 40604 | best/edit4_n12_lexicode_size141943.txt |
| E4(12,4) | 14522 | 4574 | best/editq_n12_d4_lexicode_size14522.txt |
| E4(13,3) | 500882 | 40604 | best/edit4_n13_lexicode_size500882.txt |
| E4(13,4) | 46908 | 7744 | best/editq_n13_d4_lexicode_size46908.txt |
| E4(14,4) | 153767 | 15206 | best/editq_n14_d4_lexicode_size153767.txt |

The binary and ternary codes came from local search on the whole conflict graph
(`search/fullgraph_q.py`, built with `search/trie_graph.c`), except the n = 10 (d >= 5)
and n = 11 ternary codes, which are unions of orbits of a group of symbol permutations and
word reversal (`search/orbit_q.py`, built with `search/trie_orbit.c`; each file
names its group). The quaternary codes at
n >= 12 are lexicodes from `search/big_greedy.c`.

## Citation

```bibtex
@misc{yeung2026codingagents,
  title         = {Coding Agents for Coding Theory},
  author        = {Yeung, Abraham},
  year          = {2026},
  eprint        = {2609.39081},
  archivePrefix = {arXiv},
  primaryClass  = {cs.IT},
  note          = {6th Workshop on Mathematical Reasoning and AI (MATH-AI), NeurIPS 2026}
}
```
