# OrbitReduction

A sorry-free Lean 4 (+ mathlib) formalization of the crux lemma behind this
repository's A₄-112 theorem (`proofs/a4-invariant-112.md`, Steps 1-2 of the
proof): **invariant sets are unions of orbits, and maximum invariant
independent sets reduce to the orbit graph**.

All content is in `OrbitReduction/Basic.lean`, namespace `OrbitReduction`.

## The abstract lemma

Setting: a finite group `G` (`[Group G] [Fintype G]`) acting on a type `V`
(`[MulAction G V] [DecidableEq V]`), and a conflict relation
`E : V → V → Prop`.

* `Independent E S` — no two **distinct** members of `S : Finset V` are
  `E`-related. Quantifying over distinct pairs only is the subtlety the
  informal proof flags: a word never conflicts with itself, so an orbit with
  an internal `E`-edge is excluded by *independence*, not by *invariance*.
* `Invariant G S` — `g • S = S` for every `g : G`.
* `Compatible E O₁ O₂` — no member of `O₁` is `E`-related to any member of
  `O₂`.
* `GoodFamily G E F` — `F : Finset (Finset V)` consists of `G`-orbits, each
  internally independent ("usable"), pairwise compatible.

Main results (all proved, no sorries, no axioms beyond mathlib's):

* **(a)** `Invariant.eq_biUnion_orbitF`, `Invariant.orbitF_subset_or_disjoint`,
  `invariant_iff_orbitF_subset` — an invariant set is literally the union of
  the orbits of its members, and each orbit is either contained in it or
  disjoint from it. Packaged as `part_a`.
* **(b)** `orbitF_subset_invariant_independent_iff`,
  `isOrbit_subset_invariant_independent_iff` — an orbit is contained in *some*
  invariant independent set iff it is internally independent (the witness for
  the hard direction is the orbit itself). Packaged as `part_b`.
* **(c)** `invIndepSizes_eq_orbitFamilySizes` — the set of sizes attained by
  invariant independent sets **equals** the set of weights `∑ O ∈ F, O.card`
  attained by good families; hence `sSup_invIndepSizes_eq` /
  `isGreatest_invIndepSizes_iff` — the two maxima agree, and both are attained
  (`isGreatest_sSup_invIndepSizes`, `isGreatest_sSup_orbitFamilySizes`, using
  `[Fintype V]`). Packaged as `part_c`.

Two further lemmas are the ones that genuinely use `G`-invariance of `E`
(`hE : ∀ g x y, E (g • x) (g • y) ↔ E x y`):

* `independent_orbitF_iff` — an orbit is internally independent iff no
  nontrivial translate of one representative conflicts with that
  representative (an `|O|²` test becomes a `|G|` test);
* `compatible_orbitF_iff` — two orbits are compatible iff one fixed
  representative of the first conflicts with no translate of a representative
  of the second.

These are exactly what makes the concrete orbit-graph computation cheap.

Note on hypotheses: parts (a), (b), (c) are proved **without** assuming
symmetry, irreflexivity, or `G`-invariance of `E`; the packaged theorems
`part_a`/`part_b`/`part_c` carry those hypotheses anyway for fidelity to the
informal statement. This is a strengthening of the informal lemma, not a
weakening.

## Correspondence to the concrete A₄ instance

In `proofs/a4-invariant-112.md` the lemma is instantiated with:

| abstract | concrete |
|---|---|
| `V` | quaternary words of length 6, i.e. `Fin 6 → Fin 4` (4096 words) |
| `G` | the alternating group A₄ ≤ S₄, acting symbol-wise: `σ • (x₁…x₆) = σ(x₁)…σ(x₆)` |
| `E x y` | `x ≠ y ∧ Levenshtein distance d(x,y) ≤ 2` (symmetric, irreflexive, A₄-invariant because A₄ maps edit scripts to edit scripts) |
| invariant independent `S` | A₄-invariant single-edit-correcting code `C ⊆ Σ⁶` |
| internally independent orbit | "usable" orbit |
| `GoodFamily` | family of usable, pairwise-compatible orbits |
| `part_c` | max invariant code size = max weight of a compatible usable-orbit family = 4 + 12·m |

Step 1 of the informal proof is `part_a`; the "usable orbits" reduction of
Step 2 is `part_b` + `part_c`; the "test one representative" trick implicit in
building the 270-vertex orbit graph is
`independent_orbitF_iff` / `compatible_orbitF_iff`.

## What is NOT formalized

* The concrete instance computation: the orbit census (342 orbits: 341 of
  size 12 and one of size 4), the identification of the 270 usable orbits, and
  the 17,992-edge orbit compatibility graph. These are computed twice by
  independent code paths in this repo (`search/edit4.py` and the adversarial
  referee's reimplementation) and cross-checked.
* The SAT step (`m ≤ 9`, i.e. no 10 pairwise-compatible usable size-12
  orbits): covered by the DRUP certificate
  `best/certificates/no10orbits.drup`, checked with `drat-trim`
  (`s VERIFIED`), not by Lean.
* The lower bound artifact (the 112-word code itself): verified by
  `verify.py` / `verify2.py`.

So Lean covers the *reduction* ("it suffices to solve the orbit graph
problem"), and the certificate + two independent verifiers cover the
*instance*.

## Building

```
cd lean/OrbitReduction
lake exe cache get   # fetch mathlib olean cache
lake build           # expect: zero errors, zero sorries
```

Toolchain: `leanprover/lean4:v4.32.2`, mathlib `v4.32.2`.

## Build verification (2026-09-03)

`lake build` on toolchain `leanprover/lean4:v4.32.2`: **completed successfully,
8657 jobs, exit 0** — only style-linter warnings (unused section variables in the
`Packaged` section; a duplicated `[Fintype V]` instance in `part_c`'s signature,
cosmetic and not a soundness issue).

Axiom audit via `#print axioms` on `part_a`, `part_b`, `part_c`,
`invIndepSizes_eq_orbitFamilySizes`, `independent_orbitF_iff` and
`compatible_orbitF_iff` — every one reports exactly:

    [propext, Classical.choice, Quot.sound]

i.e. the three standard mathlib axioms, with **no `sorryAx`** and no
`native_decide`. `grep` for `sorry`/`admit`/`native_decide` over the sources
returns nothing. Before this date the "sorry-free" claim in this README had
never actually been checked by a build.
