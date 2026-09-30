# SOTA record — status checks and the pick

All status checks below were run **2026-08-03**. Compute available: Apple Silicon
Mac17,2, 10 cores, 32 GB RAM, ~212 GB free disk, no Gurobi licence.

---

## NEW RESULT (2026-08-25): E4(6,3) >= 116, improving the published 114

**Claim.** There exists a quaternary code of length 6 with minimum edit distance
3 (equivalently, single-edit-correcting) having **116 codewords**. The published
best known is 114, standing since 2012. This is a lower-bound improvement of +2;
the cell remains open at [116, 176].

**Published state of the art, checked against primary sources.**

> S. Houghten, *A Table of Bounds on Optimal Fixed-Length Quaternary
> Edit-Metric Codes*, Brock University. Local snapshot
> `notes/edit-metric-lists-hunt/quaternaryedit_live_fetched_2026-08-09.html`,
> fetched 2026-08-09. Cell (n=6, d=3) reads **114-176**.
>
> J. Orth, *The Salmon Algorithm* (MSc thesis, Brock University, 2012),
> Table 9.5 "Improvements to Edit Codes": `E4(6,3): Previous Best 108, New
> Best 114`. Same page defines the quantity: "The maximum number of edit
> codewords of length n from an alphabet of q characters having a minimum
> distance d is referred to as Eq(n,d)." Underlying publication:
> Ashlock-Houghten-Brown-Orth, BioSystems 110 (2012) 1-8.
>
> Weindel-Heckel, arXiv:2504.00613v2, Table 3 gives 113 for this cell (their
> paper does not cite the edit-metric line, see draft issue 4).

**Artifacts.** Two DISTINCT 116-word codes (they share only 4 codewords), each
a union of 29 orbits of size 4:

| file | sha256(codewords) | invariance |
|---|---|---|
| `best/edit4_n6_C4_product_size116.txt` | `7068739a82ef8eab50530304b5687aa61ccacb398263b68c8ae3ae2ae7eacdf4` | C_4 = <(0 2 1 3)> acting symbol-wise |
| `best/edit4_n6_C4_diagonal_size116.txt` | `d9f49066a949b4ab4a6a3def6f0d5fb9ec58c9f822e12b12cef543b89b907757` | C_4 DIAGONAL: the two 4-cycle elements paired with word reversal |

**Verification (CLAUDE.md rule 6 -- treated as a bug until it survived all of
this).** Every check below was run in a fresh process:

1. Repo verifier `verify_edit4` (ball disjointness): True, 116, for both codes.
2. From-scratch ball disjointness, written with no repo imports: **0 colliding
   pairs** for both.
3. From-scratch Levenshtein DP: **minimum pairwise edit distance exactly 3**,
   which is the definition; equivalent to (2) since edit distance is geodesic,
   so these are two implementations of one criterion, not two criteria.
4. ENCODING CHECK: the 4096-vertex confusability graph was differentially
   tested against the definition on 400,000 random vertex pairs --
   "adjacent in graph" vs "Levenshtein <= 2" -- with **0 mismatches**. This is
   the off-by-one check rule 6 demands; a graph missing edges would manufacture
   exactly this kind of false record.
5. Both codes are independent sets in that graph (0 violations).
6. Problem statement re-read against the ORIGINAL source (Orth thesis, above),
   not against our own restatement of it.
7. Literature check: Houghten's maintained table (the authority for this
   family), plus a web search for any post-2012 improvement. None found.
8. Upper bound consistency: 116 <= 176 (literature UB), <= 212 (W-H ILP UB).
9. Independent rediscovery: two DIFFERENT symmetry classes each reached 116,
   in separate runs, producing different codes.

**How it was found, and why it was missed here for three weeks.** The natural
symmetry group is S_4 x C_2 (order 48), with 33 conjugacy classes of subgroups [CORRECTED 2026-09-04: see the appendix of this file; the original text said 30].
`search/edit4.py --sweep` parameterises cells as "H <= S_4, optionally with
reversal", which reaches only 22 of the 33; the other 11 are DIAGONAL (generated
by elements pairing a permutation WITH reversal) and were never representable.
`search/edit4_diag.py` closes that gap. Both 116 classes have all 1024 orbits of
size 4 (the group acts freely), 1021 usable.

The uncomfortable half: the C_4 *product* class is NOT one of the newly reachable
cells. It was swept on 2026-08-06 and the orbit MILP returned **100** for it
(its two siblings returned 92 and 8). So a cell containing a 116 was visited,
recorded, and reported as ~100. The 900-orbit `--exact-cap` hid 15 cells; the
MILP silently mis-solved the ones just past the cap. Both failures published a
number, which is why neither looked like a failure.

**Still open.**
- The cell is [116, 176]. This improves a lower bound only.
- Within these two classes all orbits have size 4, so attainable sizes are
  multiples of 4: **117, 118, 119 are impossible there** and the next value is
  120. Ten seeds x 10 min per class did not reach it (1 seed in 10 finds 116,
  the rest land on 112).
- ~~BEFORE CONTACTING ANYONE: re-fetch Houghten's table.~~ **CLEARED
  2026-08-26.** The live page https://www.cosc.brocku.ca/~houghten/quaternaryedit.html
  now sits behind a Cloudflare challenge and returns 403 to every automated
  request (curl and WebFetch alike; the whole site, not just this page). No
  bypass was attempted. Confirmed instead by two routes: (a) the repo owner
  loaded the page in a browser on 2026-08-26 and reports the n=6/d=3 cell reads
  114-176, matching our snapshot, along with 356-614 at n=7 and 1132-2340 at
  n=8; (b) the Wayback snapshot 2025-08-22 also reads 114-176 and is
  content-identical to our 2026-08-09 fetch apart from email obfuscation, while
  the 2011-09-30 snapshot reads 111-176 -- so the cell demonstrably IS
  maintained and has held 114 for at least a year.
- The other 28 classes have verified lower bounds only; A_4 remains the one
  certificate-grade cell.

## Status checks that changed the plan

### R1 — optimal sorting network SIZE, S(11) and S(12): **CLOSED. Do not attack.**

PROBLEMS.md said S(11) ∈ [33,35], S(12) ∈ [37,39], citing Fonollosa 2018. That is
stale by six years.

> Jannis Harder, *An Answer to the Bose–Nelson Sorting Problem for 11 and 12
> Channels*, arXiv:2012.04400 (submitted 8 Dec 2020; final 24 Jul 2022).
> **S(11) = 35 and S(12) = 39.** Method: dynamic programming generalising Van
> Voorhis from sorting networks to a broader class of comparator networks. The
> lower-bound derivation is emitted as a certificate and checked by a program
> **formally verified in Isabelle/HOL**.

The abstract states it "closes the two smallest open instances of the Bose–Nelson
sorting problem." Both cells this sprint was aimed at are gone, and the result
carries a machine-checked certificate, so there is no soft spot to re-examine.

Depth claim also stale: PROBLEMS.md said depth is solved for n ≤ 16 with n = 17
open. Ehlers & Müller settled **n = 17 at depth 10** (no depth-9 network on 17
channels exists), independently re-verified by Codish, Cruz-Filipe &
Schneider-Kamp. Source: *Sorting networks: to the end and back again*, JCSS.

### R3 — union-closed sets / Gilmer constant: **walled for the stated subgoal.**

The plan was "reproduce the numeric optimisation, then push it with better
optimisation plus certified interval arithmetic." That is aimed at a ceiling:

- Gilmer's original i.i.d. coupling caps at ≈ 0.38197 (= (3−√5)/2).
- Sawin's refinement — a convex combination of the i.i.d. and max-entropy
  couplings — is the current line, and **≈ 0.38234 is the best constant
  obtainable through that approach**, as evaluated by Cambie and by Yu.
  (arXiv:2212.12500 and follow-ups.)

So the remaining gap to 1/2 is not an optimisation problem, it is a structural
one. "Optimise harder with exact arithmetic" would re-derive a known ceiling.
Demoted; not the sprint target.

### Tier A #2 — covering designs / La Jolla: **submission home is closing.**

The LJCR announcement states there will be **no more changes after 1 March 2026**
(today is 3 Aug 2026). Data is preserved as JSON on GitHub with a Zenodo mirror,
last updated 21 Apr 2026, so it remains a fine *search* target — but the
"a win has an immediate home with a maintainer" argument that put it in Tier A no
longer holds.

### Tier A #1 — binary covering codes: **tables are hard to reach.**

- Litsyn's table (eng.tau.ac.il/~litsyn/tablecr): HTTP **403 Forbidden**.
- Kéri's table (old.sztaki.hu/~keri/codes/): HTTP **404 Not Found**.

Neither is currently readable, so the SOTA cannot be ingested per rule 4, and the
"maintainer who accepts submissions" criterion is unverified. Parked pending a
working mirror.

---

## The pick: binary s-deletion-correcting codes, small n

Chosen because it is the only Tier R candidate where (a) the problem is live,
(b) the SOTA table is published and machine-readable, (c) the verifier is exact
and instant, and (d) a win is **constructive** — a file of codewords — so rule 11
(no nonexistence without a checked certificate) never binds.

### Source of record

> Franziska Weindel and Reinhard Heckel, *LLM-Guided Search for
> Deletion-Correcting Codes*, arXiv:2504.00613v2, 5 June 2026. Tables 1 and 2.
> TU Munich / Munich Center for Machine Learning.

**Correction to an earlier automated summary:** a web search claimed this paper
found new best-known 2-deletion sizes at n = 12, 13, 16. Reading Table 2 directly
shows it did **not** — its functions give 34, 48, 194 against ReduMIS's 37, 56,
215. The paper's own text says it "matches" graph solvers on n ∈ [7,12] and
constructs *smaller* codes at n = 9, 11 and at all held-out lengths. Its genuine
contribution for s ≥ 2 is generalisation across lengths, not record sizes.

### Ingestion check (rule 4) — PASSED

`tests/test_verify.py` rebuilds VT_0(n) from the definition
(Σ i·v_i ≡ 0 mod n+1) and reproduces the published Table 1 row **exactly** for
n = 6…18: 10, 16, 30, 52, 94, 172, 316, 586, 1096, 2048, 3856, 7286, 13798.

### The open cells (best known construction → best published upper bound)

| s | n  | best known | upper bound | gap  | source of best known |
|---|----|-----------|-------------|------|----------------------|
| 1 | 12 | 316       | 320 (ILP)   | 4    | VT_0(12) |
| 1 | 13 | 586       | 592 (ILP)   | 6    | VT_0(13) |
| 1 | 14 | 1096      | 1104 (LP)   | 8    | VT_0(14) |
| 1 | 15 | 2048      | 2251 (LP)   | 203  | VT_0(15) |
| 1 | 16 | 3856      | 4202 (LP)   | 346  | VT_0(16) |
| 2 | 12 | **37**    | 53 (ILP)    | 16   | ReduMIS |
| 2 | 13 | 56        | 93 (LP)     | 37   | ReduMIS |
| 2 | 14 | 87        | 160 (LP)    | 73   | OnlineMIS/ReduMIS |
| 2 | 15 | 136       | 277 (LP)    | 141  | OnlineMIS |
| 2 | 16 | 215       | 487 (LP)    | 272  | OnlineMIS/ReduMIS |
| 3 | 12 | 12        | 16          | 4    | MIS solvers |
| 3 | 13 | 16        | 24 (LP)     | 8    | MIS solvers |

Settled (upper bound meets construction), used as the accept set:
- s=1, n ≤ 11: 10, 16, 30, 52, 94, 172 — this is the exact extent of proven VT optimality.
- s=2, n ≤ 11: 4, 5, 8, 11, 16, 24.
- s=3, n ≤ 11: 2, 2, 4, 5, 6, 8.

### Primary target

**(s=2, n=12): beat 37.** Smallest open 2-deletion cell, object is 38 strings of
length 12, verification is instant, and the 37→53 gap says at least one side is
loose. Secondary: (s=1, n=12) beat 316, which would disprove VT optimality at the
first unproven length.

### First negative result (already in hand)

`tests/test_verify.py` maximality probe: **VT_0(n) is a maximal independent set
for every n = 6…14.** No word can be added to it. So the cheapest conceivable
attack on the (1,12) cell — extend VT_0(12) to 317 — is dead. Any improvement at
(1,12) must be a structurally different code, not an extension of VT.

---

## Quaternary single-edit, n=8: size 1131 — NOT A RECORD. Novelty check killed it.

The candidate stood for about two hours. The workflow's literature sweep found
that Weindel–Heckel's Table 3 is **not the literature's best at n = 6, 7, 8**:
the DNA edit-metric literature — which their paper does not cite — has held
strictly larger codes since 2012:

> D. Ashlock, S. K. Houghten, J. A. Brown, J. Orth, *On the synthesis of DNA
> error correcting codes*, BioSystems 110 (2012) 1–8.
> **E4(6,3) ≥ 114, E4(7,3) ≥ 356, E4(8,3) ≥ 1132** (evolutionary search),
> tabulated on Houghten's Brock University page "Optimal Quaternary Edit-Metric
> Codes" (last modified 2023; no later improvement found through 2026).
> Literature upper bounds (Brock TR CS-05-07): 176 / 614 / 2340 — tighter than
> W-H's LP/ILP bounds.

A fixed-length quaternary code with minimum edit distance ≥ 3 is exactly a
single-edit-correcting code (edit distance is geodesic), and the identification
was cross-validated: W-H's greedy values at n=10/11 (11743, 40604) match the
Brock table's 2010 lexicode entries digit for digit.

So the corrected state of the art, and where we actually stand:

| n | true literature best | W-H Table 3 | ours (verified) |
|---|---|---|---|
| 6 | 114 (2012) | 113 | **116 — see the 2026-08-25 entry at the top of this file; the record fell** |
| 7 | **356** (2012) | 352 | 352 |
| 8 | **1132** (2012) | 1130 | 1131 |
| 9 | **3717** (= W-H!) | 3717 | — |
| 10 | **12489** (= W-H!) | 12489 | — |

Our 1131 is one below a 13-year-old record. All earlier "matched the record"
statements in this repo are hereby downgraded to "matched Weindel–Heckel's
table". Two silver linings, both new knowledge: (a) W-H's own n=9 and n=10
values *exceed* the older literature (3717 > 3451, 12489 > 11743) — records
their paper doesn't know it holds; (b) the 2012 lower bounds 114/356/1132 are
heuristic claims whose codeword lists are not posted anywhere we could find, so
they are, by this repo's own rules, unverified claims — worth requesting from
the authors.

The novelty check worked exactly as designed: it killed the claim before we
made it. New targets, if the campaign continues: ≥115 / ≥357 / ≥1133.

UPDATE 2026-08-25: the n=6 target was met and passed — **E4(6,3) ≥ 116**, via
prescribed-symmetry search over the diagonal subgroups of S_4 × C_2 (top of this
file). The n=7 and n=8 targets (≥357, ≥1133) remain open, and the same method is
untried there.

---

## Binary campaign final ledger (2026-08-07)

Same-engine multi-seed ReduMIS on our independently built graphs, every result
verified by both verify.py and verify2.py in fresh processes:

| cell | published best | our best | verdict |
|---|---|---|---|
| (2,12) | 37 | 36 | stands |
| (2,13) | 56 | **56** | **matched, verified twice** |
| (2,14) | 87 | 86 | stands |
| (2,15) | 136 | **136** | **matched** — notable: their own ReduMIS got 135; ours reached the OnlineMIS-held 136 |
| (3,12) | 12 | **12** | matched |
| (3,13) | 16 | **16** | matched |
| quaternary n=8 | 1132 (lit.) / 1130 (W-H) | 1131 (twice: seeds 4 and 9) | 1132 stands; ours is the largest MACHINE-VERIFIED artifact |
| quaternary n=6 | 114 (lit.) / 113 (W-H) | **116** (2026-08-25, two distinct codes) | **RECORD BROKEN** — lower bound +2, cell now [116, 176] |

Across the full campaign (both alphabets): records matched at seven cells,
exceeded at zero, with every matched object now existing as a machine-verified
artifact — which none of the published incumbents do. The consistency of
"match-but-never-exceed" across 30+ independent seeds is itself evidence these
table values sit at or very near the true optima for this engine class.

### Novelty sweep for the 116 claim (2026-09-03, pre-submission)

Independent of the table check, a targeted literature sweep for any post-2012
improvement to E4(6,3). Nothing found. Evidence, strongest first:

1. Houghten's maintained table still reads **114-176** (confirmed on the live
   page 2026-09-03 by the repo owner; automated fetch is Cloudflare-blocked).
   This is the authority for the family and is demonstrably maintained -- it
   read 111-176 in the 2011 Wayback snapshot.
2. Weindel-Heckel, arXiv:2504.00613v2 (2025) -- a RECENT paper attacking this
   exact cell with LLM-guided search -- reports **113**, i.e. below the 2012
   value, and does not cite the edit-metric line at all. A 2025 paper on the
   same cell reporting less than 114 is strong evidence that nothing above 114
   is publicly known.
3. Web sweeps for 2012-2026 improvements return only the 2012 line
   (Ashlock-Houghten-Brown-Orth; Orth's thesis) and Houghten's table.
4. The 2021/22 follow-up in the same venue and group (Houghten et al.,
   "Effective decoders for DNA codes", BioSystems, online Dec 2021) is about
   DECODERS -- side-effect machines and decoding accuracy -- not new bounds.

RESIDUAL RISK, stated plainly: that 2021/22 paper is paywalled and we read only
its abstract and secondary summaries. If it carries an updated bounds table, we
have not seen it. Mitigating factor: its senior author maintains the web table,
which still reads 114. This is the weakest link in the novelty claim and should
be named as such in any submission or correspondence.

### Corrections from adversarial review (2026-09-04)

Three independent adversarial reviews were run against the paper and repo. The
central claim survived: an independent referee rebuilt the confusability graph
from scratch (392,358 edges, matching), re-derived both codes without using our
code, and additionally established that both are MAXIMAL (no 117th word can be
added) and that the criterion matches the literature's E4(n,d) definition,
checked against Brock/Sun's own worked example. **E4(6,3) >= 116 stands.**

Corrections made:

1. **S_4 x C_2 has 33 conjugacy classes of subgroups, not 30** (22 product-form
   + 11 diagonal, not 19 + 11). Cause: `all_subgroups()` in search/edit4_diag.py
   closed over PAIRS of generators under a docstring asserting all subgroups are
   2-generated. False -- C_2^3 needs three -- so it found 91 of 98 subgroups.
   The arithmetic that exposes it immediately: S_4 has 11 conjugacy classes of
   subgroups, so "H" plus "H x C_2" is already 22. FIXED; enumeration now asserts
   98 subgroups. The three previously missed classes were searched and top out at
   **104**, so the record is unaffected (results/edit4_missed3_n6.json).
2. **A third distinct 116-word code** was sitting unchecksummed in best/ from an
   overwritten filename. Verified and promoted to
   `best/edit4_n6_C4_diagonal_b_size116.txt` (sha256 c14fbb75...). Three distinct
   116-codes now exist.
3. **verify.py did not verify the edit-metric artifacts at all** -- it crashed on
   their headers, so CLAUDE.md rule 5 was false for that whole family, and the
   only edit verifier lived inside the search module. verify.py now carries a
   self-contained edit-metric path (no import from search/) plus accept/reject
   self-tests. GROUND TRUTH WAS MODIFIED: full suite re-run, 13/13 PASS, and all
   12 artifacts in best/ now verify through `python verify.py <artifact>`.
4. Paper corrections: "only two of 25 failures cost over an hour" was false (at
   least seven did, including a 6.9 h run and a 16 h certification); "the two
   highest rungs produced nothing usable" was false (SAT and ILP produced every
   certificate the campaign has, including the A_4-112 DRAT proof and 24/24
   certified LP upper bounds -- the inversion holds for LOWER bounds only); the
   100/92/8 cells were INSIDE the cap, not past it, and the results file DID
   record `proven_optimal_invariant: false` -- the flag was lost in
   summarisation, not in logging; "25 abandoned approaches" is 23 headed entries;
   "187 run records" counts gitignored .mis pairs and certificates (149 tracked).

## SUPERSEDED (2026-09-07): E4(6,3) >= 120, not 116

An adversarial forensic review of this repository ran a stronger local-search
operator (ARW-style 2-improvements with plateau moves) on OUR OWN 1021-orbit
graph for the free-acting V_4 class -- the class our sweep had recorded as
"112 on all four seeds" -- and reached 116 on every seed within a minute and
**120 on two seeds in six**.

`best/edit4_n6_V4_size120.txt`, sha256(codewords)
adeaa3604074f5f2b42ece623e45df96f9521650b98fb11a7182321fbc5c3576.

Verified five ways in fresh processes: from-scratch ball disjointness (0
colliding pairs), from-scratch Levenshtein (min distance exactly 3), verify.py
(valid/120), verify2.py (valid/120), and 0 violations as an independent set in
Edit4Graph(6). Stabiliser is exactly V_4 (Klein, no reversal). 120 <= 176 so no
upper bound is contradicted, and 120 is a multiple of 4 as free action requires.

**E4(6,3) is now bracketed [120, 176], improving the 2012 value of 114 by six.**

What this says about our own campaign, which is worse than the bound is good:

1. Our sweep gave every class four seeds. The two classes that returned 116 were
   later given ten more each; the V_4 class was never pushed. We then reported
   "112 on all four seeds" as though it were a property of the class.
2. Re-running V_4 with OUR operator at 10 seeds x 10 min gives
   [116,116,116,116,116,112,112,112,112,112] (results/repro_v4_120.json). So the
   112 was a seed-budget artifact, not a ceiling. Our operator still never finds
   the 120.
3. The paper had built an argument on that mis-reading -- "classes with identical
   arithmetic disagree, so the arithmetic is not what produces a 116" -- which is
   false. All three free-acting classes contain a 116.
4. The calibration cell PASSED the run that under-searched the winning class. The
   device catches gross malfunction and says nothing about insufficient effort.
5. The n=7 negative (355 vs 356) came from the same operator and should now be
   treated as a statement about our search, not about the cell.

This is the second time in one campaign that an unaudited recorded conclusion
pruned the search space and cost us the answer. The first was the subgroup
parameterisation. The paper's own thesis has now happened to the paper twice.

## NEW RESULT (2026-09-07): E4(7,3) >= 359, improving the published 356

**Best artifact: 359.** `best/edit4_n7_reversal_size359.txt`, sha256(codewords)
`c45d36c5cab37fff1947853beb5749dc609963fa73d68dffd462fdd8243e6012`.
Invariant under word reversal only (164 reversal pairs + 31 palindromes;
stabiliser in S_4 x C_2 is exactly {id, reversal}). Found in class 27 of
`search/edit4_cmis.py --list` (6844 usable orbits, 760,200 orbit edges) on
seed 1 of 2 at 600 s (seed 2: 355); results/edit4_cmis_n7_n7c12c27_20260907T163452.json.
Fresh-process checks: verify.py valid/359, verify2.py valid/359, verify_editq
and verify_editq2 accept at d=3 and reject at d=4, from-scratch minimum
pairwise Levenshtein = 3, 359 distinct words of length 7. 359 <= 614.
The cell is now **[359, 614]**. The 358 below was found first, in a different
class, and is kept as a second independent artifact.

### Earlier the same day: E4(7,3) >= 358 (order-6 class)

**Claim.** A quaternary code of length 7 with minimum edit distance 3 and
**358 codewords**. Houghten's table (snapshot 2026-08-09, and the owner's
browser check of 2026-08-26) reads 356-614 for this cell; Orth's thesis Table
9.5 gives 356 as the 2012 improvement. So this is +2 on a bound standing since
2012, and the cell becomes [358, 614].

**Artifact.** `best/edit4_n7_C3xC2_size358.txt`, sha256(codewords)
`c546ee08ff2326424cee82b92a68ede37c29699e6d13d4ba02b858432e023a5a`.
Invariant under the order-6 group C_3 x C_2 = <(1 2 3)> x <reversal> (the
3-cycle on symbols 1,2,3 fixing 0, and word reversal); its stabiliser in
S_4 x C_2 is exactly this group (checked by applying all 48 elements). It is a
union of 54 orbits of size 6, 11 orbits of size 3 (reversal-symmetric words
or 3-cycle-fixed words) and the fixed word 0000000: 324 + 33 + 1 = 358.

**Verification, every check in a fresh process (rule 6):**
1. `python verify.py best/edit4_n7_C3xC2_size358.txt` -> valid True, score 358.
2. `python verify2.py ...` -> valid True, score 358.
3. `verify_editq.py` and `verify_editq2.py` (general-d verifiers, new today,
   3463-check accept/reject suite in tests/test_editq_verify.py) at d = 3:
   both accept; at d = 4 both reject (parameter shift, rule 2).
4. From-scratch Levenshtein over all 63,903 pairs: minimum distance exactly 3.
5. 358 distinct words, all length 7 over {0,1,2,3}.
6. 358 <= 614 (table UB) and <= 853 (W-H ILP UB): no upper bound contradicted.

**How.** search/edit4_cmis.py + search/mis_ils.c: the SAME orbit graph the
2026-09-05 sweep searched (class 14 in `--list` order, 2218 usable orbits,
235,561 orbit edges), where the Python operator reached 355 in 4 seeds x 370 s.
The C engine (ARW (1,2)-swaps + plateau moves, ~500 iterations/ms) reached 358
on both of its first two seeds at 600 s (results/edit4_cmis_n7_n7c14_
20260907T154400.json). Three further seeds at 900 s returned 355, 355, 355
(results/edit4_cmis_n7_n7c14b_20260907T161627.json): the engine's long runs
sit in one basin, the same pattern seen at n = 6, so seed count matters more
than seed length with this operator.

**Consequence for the paper's n=7 negative.** "355 vs 356, a search shortfall"
is now confirmed as exactly that; the cell was not tight.

### Same day, n = 6: no progress past 120, and what was tried

- C engine reproduces 120 in the V_4 class (1 seed in 6 at 60 s) and finds a
  SECOND, distinct 120 in the C_4-product class (shares 2 words with the V_4
  code): `best/edit4_n6_C4_product_size120.txt`, sha256
  `70349cf811682484b02cad9d5d2ae32f82518316ad1910953e5f60c850a99a08`,
  verify.py / verify2.py valid 120.
- 124 (31 orbits) in the three free order-4 classes: 15 seeds x 600 s all
  returned 116 (results/edit4_cmis_n6_push124_*.json); 30 seeds x 60 s per class
  (results/edit4_cmis_n6_shortbatch_*.json): V_4 reached 120 on 4 seeds, C_4
  product 116 on all 30, C_4 diagonal 120 on 2 seeds
  (best/edit4_n6_C4_diagonal_size120.txt, verify.py/verify2.py valid 120). So
  ALL THREE free-acting classes contain a 120. No 124 anywhere.
- n = 7 follow-up in the reversal class (6 seeds x 300 s, results/edit4_cmis_n7_
  n7c27b_*.json): 357, 357, 358, 358, 356, 357. The 359 stands as best; a second
  358 in this class is kept as best/edit4_n7_reversal_size358.txt.
- 122 in the (0 1)(2 3) class (2046 orbits of size 2), seeded from the V_4
  120: 5 seeds x 300 s all stay at 120 (results/edit4_cmis_n6_c29seed_*.json).
- Lovasz theta of the V_4 orbit graph is 45.1 (float, SCS, 39 min;
  results/theta_edit4_n6_orbit_class20.json): far above 31, so it neither
  rules out 124 nor is worth certifying.

## NEW RESULT (2026-09-07): E4(6,4) >= 30, improving the table's 28 (UB 32)

**Claim.** A quaternary code of length 6 with minimum edit distance **4** and
**30 codewords**. Houghten's table cell (n=6, d=4) reads 28-32 in every
snapshot we hold (2011-09-30, 2025-08-22, 2026-08-09). The cell becomes
**[30, 32]**.

**Artifact.** `best/editq_n6_d4_size30_class30.txt`, sha256(codewords)
`fc6fb24268384a77f148646c7760bb3049b57258d28386310f392f89675a270f`.
Found in the ((2 3), rev) order-2 diagonal class, but its full stabiliser in
S_4 x C_2 is the order-4 DIAGONAL group {id, (0 1)(2 3), ((2 3), rev),
((0 1), rev)}. Three classes reach 30 (the order-4 diagonal class 23, the
(0 1)(2 3) product class 29, and class 30; results/editq_n6_d4_first_
20260907T165437.json); the first sweep saved all three to ONE filename, so only
the last survives (the class-23 code, sha d1723cdb..., was overwritten;
search/editq.py now tags artifact names with the class index). Exact MaxSAT
(search/editq_exact.py, RC2) proves 30 is the OPTIMUM of class 23's orbit graph
(results/editq_n6_d4_exact_class23.json): 31 needs a different class.

**Verification (fresh processes, on the surviving file).** verify_editq.py
valid/30; verify_editq2.py (Myers bit-parallel, no shared code) valid/30;
from-scratch Wagner-Fischer minimum distance 4 (reproduce.sh step 2b); 30
distinct words of length 6.

**Also: E4(7,4) >= 72** (table 65-128), `best/editq_n7_d4_size72_class12.txt`,
sha256 `85cb1ceb68de50ccc1566cddddc8240d28b328cbb2b688e8cc0520bc74aeb391`,
stabiliser the order-8 diagonal group; two order-8 classes reach 72
(results/editq_n7_d4_first_*.json). Fresh checks: both general-d verifiers
valid/72, from-scratch minimum distance 4, distance histogram
{4: 930, 5: 920, 6: 522, 7: 184} over 2556 pairs.

**Encoding gate for general d (rule 4 in the only form available, since no
lists are published).** The same pipeline reproduces the table's EXACT cells
without exceeding any: E4(5,3) = 44 (three classes reach 44, none higher),
E4(5,4) = 11, E4(6,5) = 8, E4(6,6) = 4 (results/editq_n5_d3_gate_*.json,
editq_n5_d4_gate_*.json, editq_n6_d5_gate_*.json, editq_n6_d6_gate_*.json).
Separately, the d = 3 conflict graph built by search/lev_graph.c has exactly the
392,358 edges of Edit4Graph(6), which was itself differentially tested against
the definition (results/encoding_check.json).

**Pipeline.** search/lev_graph.c (conflict graph lev <= d-1), search/editq.py
(same 33 conjugacy classes of S_4 x C_2, orbit reduction, C engine),
verify_editq.py + verify_editq2.py (tests/test_editq_verify.py, 3463 checks).

## NEW RESULT (2026-09-07, late): E4(7,5) >= 20, improving the table's 18 (UB 24)

**Claim.** A quaternary code of length 7 with minimum edit distance 5 and 20
codewords; Houghten's table reads 18-24. Cell becomes **[20, 24]**.

**Artifact.** `best/editq_n7_d5_size20_class20.txt`, sha256
`74f19f275251e84fc6cde1a1dc58efedefb314f53d0c830e5021537f7438a553`.
Invariant under V_4 (Klein, symbols, no reversal), which acts freely here, so
the code is 5 orbits of size 4 and the next attainable value in this class is
24 = the table's upper bound. Found on all 3 seeds x 30 s in class 20 of the
general-d sweep (results/editq_n7_d5_first_*.json).

**Verification (fresh processes).** verify_editq.py valid/20, verify_editq2.py
valid/20, from-scratch Wagner-Fischer distance histogram over the 190 pairs
{5: 144, 6: 26, 7: 20} (minimum exactly 5), 20 distinct words of length 7,
stabiliser in S_4 x C_2 exactly V_4.

Words: 0003332 0031211 0110012 0222202 0311320 1001103 1112223 1120300
1200231 1333313 2000020 2133102 2213033 2221110 2332230 3022013 3111131
3223321 3302122 3330001

## SUPERSEDED (2026-09-07, late): E4(7,3) >= 364, not 359

Six more two-minute seeds of the C engine on the order-6 class (C_3 x C_2,
2218 usable orbits) returned [361, 358, 355, 355, 355, 364]
(results/edit4_cmis_n7_n7round2_*.json); the order-8 diagonal class returned
a 360 in the same batch. **E4(7,3) is now [364, 614]**, +8 over the 2012 value.

Artifact `best/edit4_n7_C3xC2_size364.txt`, sha256(codewords)
`dadc469672cd5dfe21dbe551ca884c81851a9cc661440f323cf7c65d895e34e6`: 55 orbits
of size 6, 11 of size 3, and 0000000. Fresh-process checks: verify.py valid/364,
verify2.py valid/364, both general-d verifiers accept at d=3 and reject at d=4,
from-scratch minimum pairwise Levenshtein 3, stabiliser exactly C_3 x C_2.
Next attainable size in this class is 366 (sizes are 6a+3b+c, c in {0,1}).

Note on seed statistics: 364 came on seed 25 of the class; the spread across
the six seeds (355-364) is nine codewords, so a per-class verdict from fewer
than ~10 seeds of this engine means little at n=7.

## CONTROL (2026-09-07, late): ReduMIS on the orbit graph

The comparison the referees asked for. KaMIS/ReduMIS rebuilt from source
(tools/kamis; Apple clang, no libomp, four portability patches: two
`__gnu_cxx` namespaces, GCC parallel-mode `partial_sum`, an omp.h shim; single
threaded), run with --time_limit=600 on the METIS export of the V_4 orbit graph
(graphs/edit4_n6_orbit_class20.metis, unweighted, every orbit size 4):

| instance | seeds | ReduMIS result | C engine (search/mis_ils.c) |
|---|---|---|---|
| raw n=6 graph (4096) | 1 x 600 s | 112 | 114 (best raw) |
| V_4 orbit graph (1021) | 3 x 600 s | 29 orbits = 116, all three seeds | 120 on 4 of 30 one-minute seeds |

Each ReduMIS orbit set expands to a valid 116-word code (verify_edit4 True).
Reading: the instance reduction lets ReduMIS reach 116 in one seed where the
raw graph gives 112-113, but the 120 needs the ARW-style operator; symmetry
and operator each account for part of the gap, and neither alone reaches 120
here. Logs: results/redumis/.

## NEW RESULT (2026-09-08): E4(8,3) >= 1140, improving the published 1132

**Claim.** A quaternary code of length 8 with minimum edit distance 3 and
**1140 codewords**. Houghten's table reads 1132-2340 (the 1132 from
Ashlock-Houghten-Brown-Orth 2012); Weindel-Heckel Table 3 reports 1130; this
repo's best raw-graph result was 1131 (ReduMIS, 9 seeds). Cell becomes
**[1140, 2340]**, +8.

**Artifact.** `best/edit4_n8_V4_size1140.txt`, sha256(codewords)
`59cde1b2888f012e0b97cbae823f4f9c81885a9495cf66b2e417edcdb7d28e5d`.
Invariant under V_4 (Klein, symbols, no reversal), which acts freely at n=8:
285 orbits of size 4. Found on BOTH of two ten-minute seeds of the C engine
on the V_4 orbit graph (16,381 usable orbits, 2,951,295 orbit edges, built by
search/orbit_graph.c, which was validated edge-for-edge against the Python
builder at n=6); results/edit4_cmis_n8_n8c14c20_20260907T233640.json. The
order-6 class gave 1126 on both seeds; the reversal class 1115.

**Verification (fresh processes, rule 6).** verify.py valid/1140; verify2.py
valid/1140; verify_editq.py and verify_editq2.py accept at d=3 and reject at
d=4; from-scratch Wagner-Fischer over all 649,230 pairs: minimum 3; 1140
distinct words of length 8 over {0,1,2,3}; stabiliser in S_4 x C_2 exactly V_4;
1140 <= 2340. Re-derivable by reproduce.sh step 2b.

This is the first time prescribed symmetry was tried at n=8 in this repo (the
Python orbit builder could not build the instance). Next attainable size in
the class is 1144; six more seeds are running.

Also: a second, distinct 364 at n=7 (`best/edit4_n7_C3xC2_b_size364.txt`,
sha256 5bdd83e2..., shares 25 words with the first) from 12 more two-minute
seeds: [358,361,355,355,358,358,358,358,358,364,358,358].

### n = 8 follow-up (2026-09-08): all three free classes hold a 1140; no 1144

Klein class, six more ten-minute seeds: [1140, 1132, 1136, 1132, 1136, 1132]
(results/edit4_cmis_n8_n8c20b_*.json). C_4 product class: [1136, 1140]; C_4
diagonal class: [1128, 1140] (results/edit4_cmis_n8_n8c21c25_*.json). So 1140
appears in all three free-acting classes (3 of 12 seeds over the three), the
per-seed spread is 1128-1140, and nothing reached 1144. Second artifact:
`best/edit4_n8_C4_diagonal_size1140.txt` (verify.py/verify2.py valid 1140,
stabiliser exactly the diagonal C_4). E4(7,5): ten more seeds on each free
class give 20, 16, 20; 24 not found.

## NEW RESULTS (2026-09-08): E4(8,4) >= 180 (table 176-512) and E4(8,6) >= 16 (table 14-32)

First symmetric attacks on these cells, on the V_4 orbit graphs built by
search/orbit_graph.c (16,318 usable orbits / 23.1M edges at d=4; 12,055 /
72.4M at d=6), 2 seeds x 600 s each (search/edit4_cmis.py --n 8 --d D).

- **E4(8,4) >= 180**: `best/editq_n8_d4_size180_class20.txt`, sha256
  `a6aefa7ccccda2a3e31202d5bfff4fbf9ccda4322ea21ce3dc3c0474287dc313`. Both
  seeds 180 = 45 orbits of 4. Stabiliser exactly V_4. Fresh checks:
  verify_editq valid/180, verify_editq2 valid/180, from-scratch distance
  histogram over 16,110 pairs {4: 3164, 5: 5006, 6: 4984, 7: 2522, 8: 434}.
  Cell becomes [180, 512]. Log results/edit4_cmis_n8_n8d4_20260908T022441.json.
- **E4(8,6) >= 16**: `best/editq_n8_d6_size16_class20.txt`, sha256
  `82145007728143226f9a90dff0857d0658b8f1de6016b8a69e166c38151a0c88`. Both
  seeds 16 = 4 orbits of 4; the code's full stabiliser is V_4 x <reversal>
  (order 8). Fresh checks: both verifiers valid/16, distance histogram over
  120 pairs {6: 104, 7: 6, 8: 10}. Cell becomes [16, 32].
  Words: 00020222 01122110 02223331 03103203 10033001 11131333 12012312
  13332220 20001113 21321021 22202000 23300332 30230130 31110002 32211223
  33313111. Log results/edit4_cmis_n8_n8d6_20260908T022445.json.

## NEW RESULT (2026-09-08): E4(8,5) >= 44 (table 38-96)

`best/editq_n8_d5_size44_class20.txt`, V_4 orbit graph at n=8, d=5 (15,649
usable orbits, 83.1M edges); 2 seeds x 600 s gave 40 and 44 = 11 orbits of 4.
Fresh checks: verify_editq valid/44, verify_editq2 valid/44, from-scratch
minimum pairwise Levenshtein 5 over 946 pairs, 44 distinct words of length 8,
stabiliser exactly V_4. Cell becomes [44, 96]. Log
results/edit4_cmis_n8_n8d5_20260908T022453.json.

### n=8, d=4/6 follow-up (2026-09-08): plateaus in the V_4 class
E4(8,4): four more ten-minute seeds [180, 176, 180, 180]; 184 not found.
E4(8,6): six more five-minute seeds all 16; 20 not found. Both cells now being
tried in the two C_4 classes (results/editq_n8_d{4,6}_c21.log).

### n=8 general-d plateaus across the three free classes (2026-09-08)
E4(8,5): four more seeds in V_4 [40, 40, 40, 44]; 48 not found. E4(8,4): C_4
product [176, 176], C_4 diagonal [176, 180]. E4(8,6): 16 on every seed in all
three free classes (10 seeds). So 44 / 180 / 16 are what this operator finds
in the free order-4 classes at ten minutes per seed; non-free classes (order
6, reversal) and n=9 are the next places to look. Logs: results/editq_n8_d*_c2*.log.

## SUPERSEDED (2026-09-08): E4(8,4) >= 187, not 180

The order-6 class C_3 x C_2 (8996 usable orbits at d=4, weights 1/3/6) admits
odd sizes and beat the free classes: 2 seeds x 600 s gave [184, 187].
`best/editq_n8_d4_size187_class14.txt`, sha256
`05b00335dd555e3675a64be24fa4cb3cc0e6ebed657ce95c700916cc7e2292ee`: 26 orbits
of 6, 10 of 3, and one fixed word. Fresh checks: verify_editq valid/187,
verify_editq2 valid/187, from-scratch distance histogram over 17,391 pairs
{4: 3510, 5: 5457, 6: 5520, 7: 2442, 8: 462}, stabiliser exactly C_3 x C_2.
Cell becomes **[187, 512]**. Log results/edit4_cmis_n8_n8d4c14b_*.json. (An
earlier target-capped run stopped at 181 on both seeds within seconds.)
E4(8,4) order-6 follow-up: four more ten-minute seeds all 184; the 187 (seed
11) stands and the class's seed spread is 181-187 over eight seeds.

## SUPERSEDED (2026-09-08): E4(8,4) >= 188, not 187

The order-8 diagonal class (class 12; 6500 usable orbits at d=4, weights 4/8)
gave 188 on both of 2 ten-minute seeds: `best/editq_n8_d4_size188_class12.txt`,
sha256 `f86fe59c0490a594a14d94545d6532fc3f780d470282a476a3ca6fc434521a2a`,
16 orbits of 8 and 15 of 4. Fresh checks: verify_editq valid/188, verify_editq2
valid/188, distance histogram over 17,578 pairs {4: 3618, 5: 5256, 6: 5812,
7: 2386, 8: 506}, stabiliser exactly the order-8 diagonal group. Cell is now
**[188, 512]**. The same class held the 72 at (7,4). Log
results/edit4_cmis_n8_n8d4c12_*.json; four more seeds running.

## SUPERSEDED (2026-09-08): E4(8,4) >= 196, not 188

Four more ten-minute seeds on the order-8 diagonal class: [196, 192, 188, 188].
`best/editq_n8_d4_size196_class12.txt`, sha256
`9795dd0caeb09e76e41fe778c5e7a0496b2beecbda4866507421c35184deca53`, 16 orbits
of 8 and 17 of 4. Fresh checks: verify_editq valid/196, verify_editq2
valid/196, distance histogram over 19,110 pairs {4: 3926, 5: 6028, 6: 6276,
7: 2404, 8: 476}, stabiliser exactly the order-8 diagonal group. Cell is now
**[196, 512]**, +20 over the table. Seed spread in this class so far
188-196 over six seeds; six more running. Log
results/edit4_cmis_n8_n8d4c12b_*.json.

### E4(8,6) class optima (2026-09-08)
The order-8 diagonal class at (8,6) has only 265 usable orbits; RC2 MaxSAT
proves its optimum is 16 (results/editq_n8_d6_exact_class12.json), so the
16 found there and in all three free classes is not a search shortfall in
that class. The free V_4 class (12,055 orbits, 72M edges) is too large for
RC2; whether 20 exists there, or in the unrestricted graph (UB 32), is open.

### Correction (2026-09-08): two n=8 artifact files were overwritten by later same-size codes
The follow-up seed batches saved equal-size codes to the same class-tagged
filenames, replacing the originals: `best/editq_n8_d5_size44_class20.txt` now
has sha256 `b777627b5733905edf22910b860ee33ed56227f5f2e18ba73d090b27b20fdc99` and `best/editq_n8_d6_size16_class20.txt` now has sha256
`773280cf9dc013ce9337801d25485c142dfa0ca8d2f17d10faa71c9f575d903f` (both re-verified: valid 44 / valid 16, reproduce.sh step 2b). The
hashes quoted in the entries above for those two files describe codes that no
longer exist on disk; the driver now refuses to overwrite a different code
(`_alt` suffix). The bounds are unaffected.

## SUPERSEDED (2026-09-08): E4(8,4) >= 200, not 196

Six more ten-minute seeds on the order-8 diagonal class: [196, 196, 196, 196,
200, 196]. `best/editq_n8_d4_size200_class12.txt`, sha256
`2083125e16976bdc2c73a6de79fef82e798668ebb92c84ac32b9c3a4c148d251`, 17 orbits
of 8 and 16 of 4. Fresh checks: verify_editq valid/200, verify_editq2
valid/200, distance histogram over 19,900 pairs {4: 4182, 5: 6186, 6: 6218,
7: 2742, 8: 572}, stabiliser exactly the order-8 diagonal group. Cell is now
**[200, 512]**, +24 over the table. Twelve seeds in this class so far:
188-200. Log results/edit4_cmis_n8_n8d4c12c_*.json; six more running.

## NEW RESULT (2026-09-08): E4(9,3) >= 3697, improving the table's 3451

**Claim.** A quaternary code of length 9 with minimum edit distance 3 and
**3697 codewords**; Houghten's table reads 3451-9360. Cell becomes
**[3697, 9360]**, +246. First symmetric attack at n=9 (262,144 words); the
orbit graph of the order-6 class C_3 x C_2 (41,005 usable orbits, 8,942,596
edges) was built by search/orbit_graph.c in 96 minutes and searched by the C
engine for 2 seeds x 900 s: [3688, 3697].

**Artifact.** `best/edit4_n9_C3xC2_size3697.txt`, sha256(codewords)
`82378d8b7444b0a8054d34fa530c9b0bfa4ae7001c7691902daad489f30d6982`: 600
orbits of 6, 32 of 3, and 000000000.

**Verification (fresh processes, rule 6).** verify.py valid/3697 (16 s);
verify2.py valid/3697; verify_editq.py and verify_editq2.py accept at d=3
(6.8M pairs, 69 s for the bit-parallel one) and reject at d=4; 3697 distinct
words of length 9 over {0,1,2,3}; stabiliser in S_4 x C_2 exactly C_3 x C_2;
3697 <= 9360. Encoding at n=9 differentially tested: results/encoding_check_n9.txt
(0 mismatches on 200k pairs). Log results/edit4_cmis_n9_n9c14_20260908T074748.json.
The V_4 class returned 3708 on its first seed; verification in progress.

## SUPERSEDED (2026-09-08): E4(9,3) >= 3708, not 3697

The V_4 class at n=9 (65,533 usable orbits, 15,244,728 orbit edges; orbit
graph built in 3 h 27 min) gave [3708, 3684] on 2 seeds x 900 s.
`best/edit4_n9_V4_size3708.txt`, sha256(codewords)
`091abb6ea79cb3e8e2f773f27bd900df7d97b334827c31b7a7d0edc3b1cb55d6`, 927
orbits of size 4. Fresh checks: verify.py valid/3708, verify2.py valid/3708,
verify_editq/verify_editq2 accept at d=3 (66 s) and reject at d=4, 3708
distinct words of length 9, stabiliser exactly V_4. Cell is now
**[3708, 9360]**, +257 over the table. Log
results/edit4_cmis_n9_n9c20_20260908T075317.json; four more seeds running.
E4(8,4) order-8 class, six more seeds: [196, 196, 196, 200, 196, 196]; a second,
distinct 200 (best/editq_n8_d4_size200_class12_alt2.txt). Eighteen seeds in the
class: 188-200, with 200 twice. 204 not found.
E4(9,3) V_4 class, four more 15-minute seeds: [3684, 3668, 3676, 3700]; the
3708 stands. Six seeds in the class: 3668-3708. Under machine load these seeds
ran only ~180k iterations each, so longer or more seeds are worth trying.
E4(9,3) order-6 class, four more 15-minute seeds: [3685, 3682, 3667, 3688];
the class's 3697 stands. Six seeds: 3667-3697.
E4(8,4) reversal-only class (29,434 orbits, 37.2M edges): 2 seeds x 900 s gave
[182, 180]; the order-8 diagonal class's 200 stands. Class ranking at (8,4):
order-8 diagonal 200 > order-6 187 > reversal 182 > free order-4 180.
E4(9,3) order-8 diagonal class (30,006 orbits): 2 seeds x 900 s gave
[3668, 3676]. Class ranking at (9,3) so far: V_4 3708 > C_3 x C_2 3697 >
order-8 diagonal 3676; every class beats the table's 3451 on every seed.
Housekeeping (2026-09-08): `best/editq_n8_d4_size180_class20.txt` was also
replaced by a later same-size code before the overwrite guard existed (current
sha256 in the file header; re-verified valid/180). The superseded 188 and 196
files of the order-8 class were removed; 200 (two distinct codes) remains.
E4(9,3) V_4 class, two 30-minute seeds: [3700, 3688]; eight seeds in the class
now, 3668-3708, and the 3708 stands. Longer seeds did not help here either.

### n=9 general-d instances (2026-09-08)
Order-8 diagonal class orbit graphs built at d=4 (28,627 orbits, 53.4M edges),
d=5 (13,143 / 51.6M), d=6 (5,813 / 16.5M), d=7 (142 / 10k). At d=7 that class
holds only one usable orbit-pair-free set of 4 words (RC2 optimum 4,
results/editq_n9_d7_exact_class12.json) against a table of 11-20, so d=7
needs a smaller group; the order-6 class is being built. Searches at d=4,5,6
running.
(9,7): the order-6 class has 349 usable orbits and RC2 optimum 10
(results/editq_n9_d7_exact_class14.json), one below the table's 11; the
order-8 class optimum is 4. Both are exact within their classes, so E4(9,7)
needs a smaller group; the V_4 class is being built.

## NEW RESULT (2026-09-08): E4(9,4) >= 520, improving the table's 495 (UB 2048)

Order-8 diagonal class at n=9, d=4 (28,627 usable orbits, 53.4M edges; built
in 22 min): 2 seeds x 900 s both gave 520 = 60 orbits of 8 + 10 of 4.
`best/editq_n9_d4_size520_class12.txt`, sha256
`f48053ee681fc5791c5a4f367cecd8c5bbd3d32c45dd1ffe10c194e5f5746569`.
Fresh checks: verify_editq valid/520, verify_editq2 valid/520, from-scratch
distance histogram over 134,940 pairs {4: 12898, 5: 28224, 6: 44782, 7: 35528,
8: 11802, 9: 1706}, stabiliser exactly the order-8 diagonal group. Cell becomes
**[520, 2048]**, +25. Log results/edit4_cmis_n9_n9d4c12_20260908T093241.json;
two more seeds running.
(9,7), closed within three classes: the V_4 orbit graph (36,604 usable orbits)
is COMPLETE (670M edges, 0 compatible pairs; search/complement_edges.c), so
any V_4-invariant code has one orbit = 4 words; the order-6 class optimum is
10 and the order-8 class optimum 4 (RC2). E4(9,7) stays at the table's 11-20;
the order-2 and order-3 classes (87k-131k orbits) would take ~10 h to build
with the current O(orbits x words) scan and are not attempted.
(8,6) V_4 class CLOSED (2026-09-08): the orbit graph (12,055 orbits, 72.4M
edges) has only 269,826 compatible pairs; the exact maximum-weight clique of
that compatibility graph (networkx, search/compat_clique.py) is 4 orbits = 16
words (best/editq_n8_d6_class20_exact.txt). So 16 is the optimum of every
class tried at (8,6) (V_4 exact, order-8 exact, order-6 exact 13); 20 would
need a smaller group or no symmetry.
(7,5) V_4 class CLOSED (2026-09-08): 3568 orbits, 293,622 compatible pairs;
exact maximum-weight clique = 5 orbits = 20 words
(best/editq_n7_d5_class20_exact.txt). So 24 is impossible in the V_4 class;
E4(7,5) stays [20, 24].
(7,5) C_4 classes CLOSED: C_4 diagonal exact optimum 20
(best/editq_n7_d5_class25_exact.txt), C_4 product exact optimum 16
(best/editq_n7_d5_class21_exact.txt). All three free order-4 classes are
therefore exactly 20/20/16; 24 needs a class with mixed orbit sizes or no
symmetry.

## NEW RESULT (2026-09-08): E4(9,6) >= 32, improving the table's 27 (UB 128)

Order-8 diagonal class at n=9, d=6 (5813 usable orbits, 16.5M edges): 2 seeds
x 600 s both gave 32 = 3 orbits of 8 + 2 of 4.
`best/editq_n9_d6_size32_class12.txt`, sha256
`f034a20a9dac757148b2ea2653eca03b60a5f233aa60ddde872f667e20c8ef99`.
Fresh checks: verify_editq valid/32, verify_editq2 valid/32, from-scratch
distance histogram over 496 pairs {6: 302, 7: 96, 8: 54, 9: 44}, stabiliser
exactly the order-8 diagonal group. Cell becomes **[32, 128]**, +5.
Words: 000000322 000111000 001122112 001333333 021031021 022002233 030231113
033301202 110033003 110222222 111000111 111111233 121320002 122210313
130120130 133113322 200220011 203213203 211123020 212013331 222222100
222333222 223111111 223300330 300032131 303102220 311331100 312302312
332000000 332211221 333222333 333333011.
Log results/edit4_cmis_n9_n9d6c12_20260908T101211.json; exact closure of the
class (max clique of its compatibility graph) running.
(9,6) order-8 class CLOSED: exact maximum-weight clique of its compatibility
graph = 32 (best/editq_n9_d6_class12_exact.txt), so the 32 is that class's
optimum; 36+ needs another class (V_4 being built).

## NEW RESULT (2026-09-08): E4(9,5) >= 100, improving the table's 97 (UB 384)

Order-8 diagonal class at n=9, d=5 (13,143 usable orbits, 51.6M edges): 2
seeds x 900 s both gave 100 = 10 orbits of 8 + 5 of 4.
`best/editq_n9_d5_size100_class12.txt`, sha256
`f0647e8c174766334d8c526b705cc53cc77a5d21fb6a959e1f7ec85dd51eea4a`.
Fresh checks: verify_editq valid/100, verify_editq2 valid/100, from-scratch
distance histogram over 4950 pairs {5: 1456, 6: 1690, 7: 1220, 8: 456,
9: 128}, stabiliser exactly the order-8 diagonal group. Cell becomes
**[100, 384]**, +3. Log results/edit4_cmis_n9_n9d5c12_20260908T101213.json.
The class's compatibility graph (~35M pairs) is too large for the exact
clique closure.
E4(9,4) order-8 class, two more 15-minute seeds: [512, 516]; four seeds
512-520, the 520 stands.
E4(9,5) order-8 class, two more 15-minute seeds: [100, 100] (a second distinct
100, best/editq_n9_d5_size100_class12_alt2.txt). Four seeds all 100; 104 not
found; the order-6 class is being tried.
E4(9,5) order-6 class (20,893 orbits, 106M edges): 2 seeds x 900 s both 100
(best/editq_n9_d5_size100_class14.txt, a third distinct 100 in a different
class). Six seeds over two classes all 100; the cell looks settled at 100 for
this operator.
E4(9,4) order-6 class (39,793 orbits, 77.4M edges, built in 26 min): 2 seeds x
900 s gave [517, 514]; the order-8 class's 520 stands.
E4(8,5) V_4 class, four more ten-minute seeds all 40; ten seeds in the class
40-44 with 44 twice. 48 not found; 44 stands.

## SUPERSEDED (2026-09-08): E4(8,3) >= 1144, not 1140

Four more ten-minute seeds on the n=8 V_4 class: [1144, 1140, 1136, 1140].
`best/edit4_n8_V4_size1144.txt`, sha256(codewords)
`366ae67906745711600b3929d1ce1d546aa79f4b9215201f9de9aa6c37ec2de2`, 286 orbits
of size 4. Fresh checks: verify.py valid/1144, verify2.py valid/1144, both
general-d verifiers accept at d=3 and reject at d=4, from-scratch minimum
pairwise Levenshtein 3 over 653,796 pairs, stabiliser exactly V_4. Cell is now
**[1144, 2340]**, +12 over the table. Sixteen seeds in this class: 1128-1144.
Log results/edit4_cmis_n8_n8c20d_20260908T115622.json; four more seeds running.
E4(9,4) order-8 class, two more seeds: [520, 516]; a second distinct 520
(best/editq_n9_d4_size520_class12_alt2.txt). Six seeds 512-520 with 520 twice;
524 not found.
E4(9,3) V_4 class, four more 15-minute seeds: [3676, 3696, 3696, 3680]; twelve
seeds in the class 3668-3708, the 3708 stands.

## SUPERSEDED (2026-09-08): E4(8,3) >= 1148, not 1144

Four more ten-minute seeds on the n=8 V_4 class: [1136, 1148, 1136, 1124].
`best/edit4_n8_V4_size1148.txt`, sha256(codewords)
`e98fe9d121189c38208634f4b8978be329f6aeee0f2d2bedd97105bb87c57fb8`, 287 orbits
of size 4. Fresh checks: verify.py valid/1148, verify2.py valid/1148, both
general-d verifiers accept at d=3 and reject at d=4, from-scratch minimum
pairwise Levenshtein 3, stabiliser exactly V_4. Cell is now **[1148, 2340]**,
+16 over the table. Twenty seeds in this class: 1124-1148 (spread 24). Log
results/edit4_cmis_n8_n8c20e_20260908T121720.json; four more running.
n=8, d=3 follow-up batches (2026-09-08 afternoon): V_4 [1128, 1140, 1136, 1140]
(24 seeds in the class, 1124-1148; the 1148 stands); C_4 product [1132, 1128,
1132, 1144] (a 1144 in a second class, best/edit4_n8_C4_product_size1144.txt,
verify.py valid); C_4 diagonal [1108, 1108, 1136, 1136]. 1152 not found.
E4(8,4) order-8 class, six more ten-minute seeds: [196, 192, 192, 192, 192,
188]; 24 seeds in the class 188-200 with 200 twice; 204 not found.
E4(7,3) order-6 class, eight more two-minute seeds: [358, 358, 358, 358, 358,
355, 358, 358]; 26 seeds in the class 355-364 with 364 twice; 366 not found.
E4(9,3) V_4 class, six more 15-minute seeds: [3612, 3616, 3684, 3668, 3680,
3672] (the first two ran at ~80k iterations under heavy machine load);
eighteen seeds in the class 3612-3708, the 3708 stands.
E4(9,4) order-8 class, four more seeds: [512, 520, 512, 512] (a third distinct
520, best/editq_n9_d4_size520_class12_alt3.txt). Ten seeds 512-520 with 520
three times; 524 not found.
E4(6,3) V_4 class, sixty more one-minute seeds (2026-09-08): 120 three times,
116 otherwise, one 112; ninety seeds in the class now, 120 seven times, 124
never. The n=6 cell stays [120, 176].
n=8, order-4 diagonal class with fixed points (class 23, the (6,4) winner):
(8,4) 188 on both seeds (13,114 orbits), (8,5) 42 on both (4184 orbits),
(8,6) 16 on both (2594 orbits). Below or equal to the records 200 / 44 / 16.
n=8, d=3 fair trials: reversal class [1119, 1120, 1115, 1117] (4 x 600 s);
the class is genuinely weaker than V_4 (1148) at n=8.
(8,6) class 23 CLOSED: exact maximum-weight clique = 16 (7 orbits of sizes
1/2/4; best/editq_n8_d6_class23_exact.txt). Four classes now exactly 16 or
less at (8,6): V_4 16, order-8 16, class 23 16, order-6 13.
n=8, d=3, order-6 class fair trial (4 x 600 s): [1135, 1129, 1129, 1123]. Class
ranking at (8,3): V_4 1148 > C_4 product 1144 > C_4 diagonal 1140 > order-6
1135 > reversal 1120.

## SUPERSEDED (2026-09-08): E4(9,3) >= 3724, not 3708

The C_4-product class at n=9 (65,533 usable orbits, 15,244,600 edges; built in
90 min) gave [3688, 3724] on its first 2 seeds x 900 s.
`best/edit4_n9_C4_product_size3724.txt`, sha256(codewords)
`027e7b4ad5bcd2ff698d28d46bed74a61609e98cc6dd4edef9eb61b583e981f1`, 931
orbits of size 4. Fresh checks: verify.py valid/3724, verify2.py valid/3724,
verify_editq (155 s) and verify_editq2 (39 s) accept at d=3 and reject at d=4,
3724 distinct words of length 9, stabiliser exactly C_4. Cell is now
**[3724, 9360]**, +273 over the table. Log
results/edit4_cmis_n9_n9c21_20260908T173346.json; two more seeds and the C_4
diagonal class running.
E4(9,3) C_4-diagonal class (65,531 orbits; built in 38 min): 2 seeds x 900 s
gave [3692, 3668]. All three free order-4 classes at n=9 now searched: C_4
product 3724 (2 seeds), V_4 3708 (18 seeds), C_4 diagonal 3692 (2 seeds).
E4(9,3) C_4-product class, two more seeds: [3696, 3712]; four seeds 3688-3724,
the 3724 stands.
E4(9,3) C_4-diagonal class, two more seeds: [3696, 3688]; four seeds 3668-3696.
E4(9,3) C_4-product class, seeds 20-21: [3692, 3688]; six seeds 3688-3724, the
3724 (seed 2) stands; four more running.
E4(9,3) C_4-product class, seeds 30-33: [3680, 3696, 3684, 3700]; ten seeds
3680-3724, the 3724 stands. Class ranking at (9,3): C_4 product 3724 > V_4
3708 > C_3 x C_2 3697 > C_4 diagonal 3696 > order-8 diagonal 3676.

## SUPERSEDED (2026-09-08, evening): E4(9,3) >= 3744 by a DIFFERENT OPERATOR (ReduMIS)

ReduMIS (tools/kamis, single-threaded, --time_limit=900) on the METIS export
of the n=9 C_4-DIAGONAL orbit graph (65,531 orbits, all size 4): seeds 1, 2
gave 936 and 940 orbits = 3744 and 3760 words. On the same graph the ARW
engine's four seeds gave 3668-3696, and on the V_4 / C_4-product graphs
ReduMIS gave 929/929 and 928/930 (3712-3720) against ARW 3612-3724. So at
n=9 ReduMIS is the stronger operator, the reverse of n=6 (116 vs 120), and
the class it prefers differs from ARW's.

`best/edit4_n9_class25_redumis_s1.txt` (3744), sha256(codewords)
`7ffc2b48943f6bdd7a1e1c897f2d5a031b3d924353dc7c3b218292ac7f3cc87e`: verify.py
valid/3744, verify2.py valid/3744, verify_editq/verify_editq2 accept, 3744
distinct words of length 9, from-scratch minimum pairwise Levenshtein 3,
stabiliser exactly the diagonal C_4. The 3760 (seed 2) is being verified.
NuMVC (search/mis_numvc.c), the third operator tried, gave 3516-3524.
Logs: results/redumis/n9_c2*_s*.log, results/numvc_n9_c21_s*.txt.

## SUPERSEDED (2026-09-08, evening): E4(9,3) >= 3760

`best/edit4_n9_C4_diagonal_size3760.txt`, sha256(codewords)
`9c7c0ce7c75fbf30820c6af6863e58b28b867e0320ff2476abe864f29bcffb17`: ReduMIS
seed 2 on the C_4-diagonal orbit graph, 940 orbits of size 4. Fresh checks:
verify.py valid/3760, verify2.py valid/3760, both general-d verifiers accept
(mis_expand.py), 3760 distinct words of length 9, from-scratch minimum
pairwise Levenshtein 3, stabiliser exactly the diagonal C_4. Cell is now
**[3760, 9360]**, +309 over the table. Three more ReduMIS seeds on this class
are at 934-939 orbits (interim). At n=8, ReduMIS on the V_4 orbit graph gave
283 and 281 orbits (1132, 1124; both expand to valid codes) against the ARW
engine's 1148, so the operator ranking reverses between n=8 and n=9.

### ReduMIS follow-ups (2026-09-08, late)
n=9 C_4-diagonal orbit graph, ReduMIS seeds 3, 4 (900 s) and 5 (1800 s): 937,
939, 940 orbits = 3748, 3756, 3760 words (results/redumis/n9_c25_s{3,4,5}.log).
Five ReduMIS seeds in that class: 3744-3760, with 3760 twice; the ARW engine's
best there was 3696 over four seeds. n=8 comparisons on the same graphs: at
(8,3) ReduMIS gives 1132/1124 on V_4 and 1136/1140 on C_4 diagonal (ARW 1148
and 1140); at (8,4) ReduMIS gives 176/196 on V_4 where ARW gave 180 (the
order-8 class's 200 stands). So ReduMIS is the stronger operator at n=9 and
on the (8,4) V_4 graph, ARW at (8,3); neither dominates.
Second distinct 3760 at (9,3): ReduMIS seed 5 (1800 s) on the C_4-diagonal
orbit graph, `best/edit4_n9_C4_diagonal_b_size3760.txt`, sha256
`8dcf0bb29bbb5185e964ee8c65b1280740d0cdc18494a55b686032bda2eddc69`; shares 44
words with the first. Both general-d verifiers accept.
ReduMIS seeds 6, 7 on the n=9 C_4-diagonal orbit graph: 939, 938 orbits (3756,
3752). Seven seeds 3744-3760 with 3760 twice; 3764 not found. E4(9,3) stays
[3760, 9360].
ReduMIS on the n=8 V_4 orbit graphs at d=5 and d=6 (2 x 600 s each): 40, 40
and 16, 16 (results/redumis/n8d5_c20_s*.log, n8d6_c20_s*.log); below the ARW
44 at (8,5) and equal to the exact optimum 16 at (8,6). The order-8 diagonal
class is weighted, so ReduMIS cannot be run there. No cell moved.
