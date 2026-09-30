/-
# Orbit reduction for invariant independent sets

This file formalizes the crux lemma behind the repository's `A₄`-112 theorem
(`proofs/a4-invariant-112.md`), at the natural level of abstraction.

Setting: a finite group `G` acting on a finite type `V`, together with a
symmetric irreflexive `G`-invariant "conflict" relation `E : V → V → Prop`.

* `S : Finset V` is **independent** when no two *distinct* members are
  `E`-related.  The restriction to distinct pairs is the subtlety the informal
  proof flags: a word never conflicts with itself, so an orbit carrying an
  internal `E`-edge is ruled out by *independence*, not by *invariance*.
* `S` is **invariant** when `g • S = S` for every `g : G`.

We prove:

* (a) `invariant_iff_orbitF_subset`, `Invariant.eq_biUnion_orbitF`,
      `Invariant.orbitF_subset_or_disjoint` — invariant sets are exactly the
      unions of orbits;
* (b) `orbitF_subset_invariant_independent_iff` — an orbit sits inside *some*
      invariant independent set iff it is internally independent;
* (c) `invIndepSizes_eq_orbitFamilySizes` and `sSup_invIndepSizes_eq` — the
      achievable sizes of invariant independent sets are *exactly* the
      achievable weights of families of pairwise-compatible internally
      independent orbits, hence the two maxima agree.

Two extra lemmas, `independent_orbitF_iff` and `compatible_orbitF_iff`, are the
ones that actually need `G`-invariance of `E`: they say usability and
compatibility can be tested against a single orbit representative, which is what
makes the concrete `A₄` computation cheap.
-/
import Mathlib

universe u v

namespace OrbitReduction

variable {G : Type u} {V : Type v}
variable [Group G] [Fintype G] [DecidableEq V] [MulAction G V]

/-! ## Orbits as `Finset`s -/

/-- The `G`-orbit of `x`, as a `Finset`. -/
def orbitF (G : Type u) [Group G] [Fintype G] {V : Type v} [DecidableEq V] [MulAction G V]
    (x : V) : Finset V :=
  Finset.univ.image (fun g : G => g • x)

@[simp]
theorem mem_orbitF {x y : V} : y ∈ orbitF G x ↔ ∃ g : G, g • x = y := by
  simp [orbitF]

theorem self_mem_orbitF (x : V) : x ∈ orbitF G x :=
  mem_orbitF.2 ⟨1, one_smul _ _⟩

theorem smul_mem_orbitF (g : G) (x : V) : g • x ∈ orbitF G x :=
  mem_orbitF.2 ⟨g, rfl⟩

theorem orbitF_nonempty (x : V) : (orbitF G x).Nonempty :=
  ⟨x, self_mem_orbitF x⟩

/-- Membership in an orbit determines the orbit. -/
theorem orbitF_eq_of_mem {x y : V} (h : y ∈ orbitF G x) : orbitF G y = orbitF G x := by
  obtain ⟨g, rfl⟩ := mem_orbitF.1 h
  ext z
  simp only [mem_orbitF]
  constructor
  · rintro ⟨h', rfl⟩
    exact ⟨h' * g, by rw [mul_smul]⟩
  · rintro ⟨h', rfl⟩
    exact ⟨h' * g⁻¹, by rw [mul_smul, inv_smul_smul]⟩

/-- Distinct orbits are disjoint. -/
theorem orbitF_disjoint {x y : V} (h : orbitF G x ≠ orbitF G y) :
    Disjoint (orbitF G x) (orbitF G y) := by
  rw [Finset.disjoint_left]
  intro a ha ha'
  exact h (by rw [← orbitF_eq_of_mem ha, ← orbitF_eq_of_mem ha'])

/-- `O` is an orbit. -/
def IsOrbit (G : Type u) [Group G] [Fintype G] {V : Type v} [DecidableEq V] [MulAction G V]
    (O : Finset V) : Prop :=
  ∃ x : V, O = orbitF G x

theorem isOrbit_orbitF (x : V) : IsOrbit G (orbitF G x) := ⟨x, rfl⟩

theorem IsOrbit.nonempty {O : Finset V} (h : IsOrbit G O) : O.Nonempty := by
  obtain ⟨x, rfl⟩ := h; exact orbitF_nonempty x

/-- Two distinct orbits are disjoint. -/
theorem IsOrbit.disjoint {O₁ O₂ : Finset V} (h₁ : IsOrbit G O₁) (h₂ : IsOrbit G O₂)
    (hne : O₁ ≠ O₂) : Disjoint O₁ O₂ := by
  obtain ⟨x, rfl⟩ := h₁
  obtain ⟨y, rfl⟩ := h₂
  exact orbitF_disjoint hne

/-! ## Independence and invariance -/

/-- `S` is **independent**: no two *distinct* members are `E`-related.
Note the quantification over distinct pairs only — `E x x` is never consulted. -/
def Independent (E : V → V → Prop) (S : Finset V) : Prop :=
  ∀ ⦃x⦄, x ∈ S → ∀ ⦃y⦄, y ∈ S → x ≠ y → ¬ E x y

/-- `S` is **invariant**: `g • S = S` for every `g : G`. -/
def Invariant (G : Type u) [Group G] [Fintype G] {V : Type v} [DecidableEq V] [MulAction G V]
    (S : Finset V) : Prop :=
  ∀ g : G, S.image (fun x => g • x) = S

/-- Orbits `O₁`, `O₂` are **compatible**: no member of one is `E`-related to any member of the
other. -/
def Compatible (E : V → V → Prop) (O₁ O₂ : Finset V) : Prop :=
  ∀ ⦃x⦄, x ∈ O₁ → ∀ ⦃y⦄, y ∈ O₂ → ¬ E x y

theorem Independent.subset {E : V → V → Prop} {S T : Finset V} (h : Independent E T)
    (hst : S ⊆ T) : Independent E S :=
  fun _ hx _ hy hne => h (hst hx) (hst hy) hne

theorem independent_empty (E : V → V → Prop) : Independent E (∅ : Finset V) := by
  intro x hx; simp at hx

theorem invariant_empty : Invariant G (∅ : Finset V) := fun _ => Finset.image_empty _

/-- Invariance is equivalent to closure under the action. -/
theorem invariant_iff_smul_mem {S : Finset V} :
    Invariant G S ↔ ∀ (g : G) ⦃x⦄, x ∈ S → g • x ∈ S := by
  constructor
  · intro h g x hx
    rw [← h g]
    exact Finset.mem_image_of_mem _ hx
  · intro h g
    have hsub : S.image (fun x => g • x) ⊆ S := by
      intro y hy
      obtain ⟨x, hx, rfl⟩ := Finset.mem_image.1 hy
      exact h g hx
    have hcard : S.card ≤ (S.image (fun x => g • x)).card := by
      rw [Finset.card_image_of_injective _ (MulAction.injective g)]
    exact Finset.eq_of_subset_of_card_le hsub hcard

theorem Invariant.smul_mem {S : Finset V} (h : Invariant G S) (g : G) {x : V} (hx : x ∈ S) :
    g • x ∈ S :=
  invariant_iff_smul_mem.1 h g hx

/-! ## (a) Invariant sets are unions of orbits -/

/-- **(a)** `S` is invariant iff it contains the full orbit of each of its members. -/
theorem invariant_iff_orbitF_subset {S : Finset V} :
    Invariant G S ↔ ∀ ⦃x⦄, x ∈ S → orbitF G x ⊆ S := by
  rw [invariant_iff_smul_mem]
  constructor
  · intro h x hx y hy
    obtain ⟨g, rfl⟩ := mem_orbitF.1 hy
    exact h g hx
  · intro h g x hx
    exact h hx (smul_mem_orbitF g x)

theorem Invariant.orbitF_subset {S : Finset V} (h : Invariant G S) {x : V} (hx : x ∈ S) :
    orbitF G x ⊆ S :=
  invariant_iff_orbitF_subset.1 h hx

/-- **(a)** An invariant set is literally the union of the orbits of its members. -/
theorem Invariant.eq_biUnion_orbitF {S : Finset V} (h : Invariant G S) :
    S = S.biUnion (fun x => orbitF G x) := by
  apply Finset.Subset.antisymm
  · intro x hx
    exact Finset.mem_biUnion.2 ⟨x, hx, self_mem_orbitF x⟩
  · intro y hy
    obtain ⟨x, hx, hy'⟩ := Finset.mem_biUnion.1 hy
    exact h.orbitF_subset hx hy'

/-- **(a)** Each orbit is either contained in an invariant set or disjoint from it. -/
theorem Invariant.orbitF_subset_or_disjoint {S : Finset V} (h : Invariant G S) (x : V) :
    orbitF G x ⊆ S ∨ Disjoint (orbitF G x) S := by
  by_cases hx : ∃ a, a ∈ orbitF G x ∧ a ∈ S
  · obtain ⟨a, ha, haS⟩ := hx
    left
    rw [← orbitF_eq_of_mem ha]
    exact h.orbitF_subset haS
  · right
    rw [Finset.disjoint_left]
    intro a ha haS
    exact hx ⟨a, ha, haS⟩

/-- Orbits are invariant. -/
theorem invariant_orbitF (x : V) : Invariant G (orbitF G x) := by
  rw [invariant_iff_smul_mem]
  intro g y hy
  obtain ⟨h, rfl⟩ := mem_orbitF.1 hy
  exact mem_orbitF.2 ⟨g * h, by rw [mul_smul]⟩

theorem IsOrbit.invariant {O : Finset V} (h : IsOrbit G O) : Invariant G O := by
  obtain ⟨x, rfl⟩ := h; exact invariant_orbitF x

/-- A union of invariant sets is invariant. -/
theorem invariant_biUnion {ι : Type*} [DecidableEq ι] {F : Finset ι} {t : ι → Finset V}
    (h : ∀ i ∈ F, Invariant G (t i)) : Invariant G (F.biUnion t) := by
  rw [invariant_iff_smul_mem]
  intro g x hx
  obtain ⟨i, hi, hx'⟩ := Finset.mem_biUnion.1 hx
  exact Finset.mem_biUnion.2 ⟨i, hi, (h i hi).smul_mem g hx'⟩

/-! ## (b) Which orbits can appear -/

/-- **(b)** An orbit is contained in *some* invariant independent set iff it is internally
independent.  (The witness for the hard direction is the orbit itself.) -/
theorem orbitF_subset_invariant_independent_iff (E : V → V → Prop) (x : V) :
    (∃ S : Finset V, Invariant G S ∧ Independent E S ∧ orbitF G x ⊆ S) ↔
      Independent E (orbitF G x) := by
  constructor
  · rintro ⟨S, _, hindep, hsub⟩
    exact hindep.subset hsub
  · intro h
    exact ⟨orbitF G x, invariant_orbitF x, h, Finset.Subset.refl _⟩

/-- **(b)**, phrased for an abstract orbit. -/
theorem isOrbit_subset_invariant_independent_iff (E : V → V → Prop) {O : Finset V}
    (hO : IsOrbit G O) :
    (∃ S : Finset V, Invariant G S ∧ Independent E S ∧ O ⊆ S) ↔ Independent E O := by
  obtain ⟨x, rfl⟩ := hO
  exact orbitF_subset_invariant_independent_iff E x

/-! ## (c) Reduction to the orbit graph -/

/-- A family of orbits that are internally independent and pairwise compatible.  Its members are
exactly the "usable, pairwise-compatible orbits" of the informal proof. -/
structure GoodFamily (G : Type u) [Group G] [Fintype G] {V : Type v} [DecidableEq V]
    [MulAction G V] (E : V → V → Prop) (F : Finset (Finset V)) : Prop where
  /-- every member is a `G`-orbit -/
  isOrbit : ∀ ⦃O⦄, O ∈ F → IsOrbit G O
  /-- every member is internally independent ("usable") -/
  indep : ∀ ⦃O⦄, O ∈ F → Independent E O
  /-- distinct members are compatible -/
  compat : ∀ ⦃O₁⦄, O₁ ∈ F → ∀ ⦃O₂⦄, O₂ ∈ F → O₁ ≠ O₂ → Compatible E O₁ O₂

/-- Members of a good family are pairwise disjoint, so its weight is the size of its union. -/
theorem GoodFamily.card_biUnion {E : V → V → Prop} {F : Finset (Finset V)}
    (h : GoodFamily G E F) : (F.biUnion id).card = ∑ O ∈ F, O.card :=
  Finset.card_biUnion fun _ h₁ _ h₂ hne => (h.isOrbit h₁).disjoint (h.isOrbit h₂) hne

/-- The union of a good family is invariant and independent. -/
theorem GoodFamily.invariant {E : V → V → Prop} {F : Finset (Finset V)} (h : GoodFamily G E F) :
    Invariant G (F.biUnion id) :=
  invariant_biUnion fun _ hO => (h.isOrbit hO).invariant

theorem GoodFamily.independent {E : V → V → Prop} {F : Finset (Finset V)}
    (h : GoodFamily G E F) : Independent E (F.biUnion id) := by
  intro x hx y hy hne
  obtain ⟨O₁, hO₁, hx'⟩ := Finset.mem_biUnion.1 hx
  obtain ⟨O₂, hO₂, hy'⟩ := Finset.mem_biUnion.1 hy
  by_cases heq : O₁ = O₂
  · subst heq
    exact h.indep hO₁ hx' hy' hne
  · exact h.compat hO₁ hO₂ heq hx' hy'

/-- The orbit decomposition of an invariant independent set is a good family. -/
theorem goodFamily_image_orbitF {E : V → V → Prop} {S : Finset V} (hinv : Invariant G S)
    (hindep : Independent E S) : GoodFamily G E (S.image (fun x => orbitF G x)) where
  isOrbit := by
    intro O hO
    obtain ⟨x, _, rfl⟩ := Finset.mem_image.1 hO
    exact isOrbit_orbitF x
  indep := by
    intro O hO
    obtain ⟨x, hx, rfl⟩ := Finset.mem_image.1 hO
    exact hindep.subset (hinv.orbitF_subset hx)
  compat := by
    intro O₁ hO₁ O₂ hO₂ hne x hx y hy
    obtain ⟨a, ha, rfl⟩ := Finset.mem_image.1 hO₁
    obtain ⟨b, hb, rfl⟩ := Finset.mem_image.1 hO₂
    have hdisj := orbitF_disjoint hne
    have hxy : x ≠ y := by
      rintro rfl
      exact (Finset.disjoint_left.1 hdisj hx) hy
    exact hindep (hinv.orbitF_subset ha hx) (hinv.orbitF_subset hb hy) hxy

theorem biUnion_image_orbitF {S : Finset V} (hinv : Invariant G S) :
    (S.image (fun x => orbitF G x)).biUnion id = S := by
  apply Finset.Subset.antisymm
  · intro y hy
    obtain ⟨O, hO, hy'⟩ := Finset.mem_biUnion.1 hy
    obtain ⟨x, hx, rfl⟩ := Finset.mem_image.1 hO
    exact hinv.orbitF_subset hx hy'
  · intro x hx
    exact Finset.mem_biUnion.2
      ⟨orbitF G x, Finset.mem_image_of_mem _ hx, self_mem_orbitF x⟩

variable (G) in
/-- Sizes attained by invariant independent sets. -/
def invIndepSizes (E : V → V → Prop) : Set ℕ :=
  { n | ∃ S : Finset V, Invariant G S ∧ Independent E S ∧ S.card = n }

variable (G) in
/-- Weights attained by families of pairwise-compatible internally independent orbits. -/
def orbitFamilySizes (E : V → V → Prop) : Set ℕ :=
  { n | ∃ F : Finset (Finset V), GoodFamily G E F ∧ ∑ O ∈ F, O.card = n }

/-- **(c)**, sharp form: the two sets of attainable values coincide. -/
theorem invIndepSizes_eq_orbitFamilySizes (E : V → V → Prop) :
    invIndepSizes G E = orbitFamilySizes G E := by
  ext n
  constructor
  · rintro ⟨S, hinv, hindep, rfl⟩
    refine ⟨S.image (fun x => orbitF G x), goodFamily_image_orbitF hinv hindep, ?_⟩
    rw [← (goodFamily_image_orbitF hinv hindep).card_biUnion, biUnion_image_orbitF hinv]
  · rintro ⟨F, hF, rfl⟩
    exact ⟨F.biUnion id, hF.invariant, hF.independent, hF.card_biUnion⟩

/-! ### The maxima agree -/

variable [Fintype V]

theorem invIndepSizes_nonempty (E : V → V → Prop) : (invIndepSizes G E).Nonempty :=
  ⟨0, ∅, invariant_empty, independent_empty E, Finset.card_empty⟩

theorem invIndepSizes_bddAbove (E : V → V → Prop) : BddAbove (invIndepSizes G E) := by
  refine ⟨Fintype.card V, ?_⟩
  rintro n ⟨S, -, -, rfl⟩
  simpa using Finset.card_le_univ S

/-- The supremum of the attainable sizes is attained: it is a genuine maximum. -/
theorem isGreatest_sSup_invIndepSizes (E : V → V → Prop) :
    IsGreatest (invIndepSizes G E) (sSup (invIndepSizes G E)) :=
  ⟨Nat.sSup_mem (invIndepSizes_nonempty E) (invIndepSizes_bddAbove E),
    fun _ hx => le_csSup (invIndepSizes_bddAbove E) hx⟩

theorem isGreatest_sSup_orbitFamilySizes (E : V → V → Prop) :
    IsGreatest (orbitFamilySizes G E) (sSup (orbitFamilySizes G E)) := by
  rw [← invIndepSizes_eq_orbitFamilySizes]
  exact isGreatest_sSup_invIndepSizes E

/-- **(c)** `max { |S| : S invariant independent } = max { ∑ |O| : O ranging over a family of
pairwise-compatible internally independent orbits }`, and both suprema are attained. -/
theorem sSup_invIndepSizes_eq (E : V → V → Prop) :
    sSup (invIndepSizes G E) = sSup (orbitFamilySizes G E) :=
  congrArg sSup (invIndepSizes_eq_orbitFamilySizes E)

/-- **(c)**, `IsGreatest` form: `n` is the maximum size of an invariant independent set iff it is
the maximum weight of a good orbit family. -/
theorem isGreatest_invIndepSizes_iff (E : V → V → Prop) (n : ℕ) :
    IsGreatest (invIndepSizes G E) n ↔ IsGreatest (orbitFamilySizes G E) n := by
  rw [invIndepSizes_eq_orbitFamilySizes]

/-! ## Representative reduction (this is where `G`-invariance of `E` is used) -/

section Invariance

variable {E : V → V → Prop} (hE : ∀ (g : G) (x y : V), E (g • x) (g • y) ↔ E x y)
include hE

/-- An orbit is internally independent iff no nontrivial translate of a representative conflicts
with that representative.  This turns an `|O|²` test into an `|G|` test. -/
theorem independent_orbitF_iff (x : V) :
    Independent E (orbitF G x) ↔ ∀ g : G, g • x ≠ x → ¬ E x (g • x) := by
  constructor
  · intro h g hg
    exact h (self_mem_orbitF x) (smul_mem_orbitF g x) (Ne.symm hg)
  · intro h a ha b hb hne
    obtain ⟨g, rfl⟩ := mem_orbitF.1 ha
    obtain ⟨k, rfl⟩ := mem_orbitF.1 hb
    have key : E (g • x) (k • x) ↔ E x ((g⁻¹ * k) • x) := by
      rw [← hE g⁻¹ (g • x) (k • x), inv_smul_smul, ← mul_smul]
    rw [key]
    refine h (g⁻¹ * k) ?_
    intro hcontra
    apply hne
    calc g • x = g • ((g⁻¹ * k) • x) := by rw [hcontra]
      _ = k • x := by rw [smul_smul, mul_inv_cancel_left]

/-- Two orbits are compatible iff a fixed representative of the first conflicts with no member of
the second.  This turns an `|O₁| · |O₂|` test into an `|G|` test. -/
theorem compatible_orbitF_iff (x y : V) :
    Compatible E (orbitF G x) (orbitF G y) ↔ ∀ g : G, ¬ E x (g • y) := by
  constructor
  · intro h g
    exact h (self_mem_orbitF x) (smul_mem_orbitF g y)
  · intro h a ha b hb
    obtain ⟨g, rfl⟩ := mem_orbitF.1 ha
    obtain ⟨k, rfl⟩ := mem_orbitF.1 hb
    have key : E (g • x) (k • y) ↔ E x ((g⁻¹ * k) • y) := by
      rw [← hE g⁻¹ (g • x) (k • y), inv_smul_smul, ← mul_smul]
    rw [key]
    exact h (g⁻¹ * k)

end Invariance

/-! ## Packaged statement

The three parts, restated in one place with *all* the ambient hypotheses of the informal lemma
present (symmetry and irreflexivity of `E`, and `G`-invariance of `E`).  These hypotheses are
carried for fidelity to the informal statement; the proofs of (a), (b), (c) above do not need
them, which is a strengthening, not a weakening.  `G`-invariance is genuinely used by
`independent_orbitF_iff` / `compatible_orbitF_iff` above.
-/

section Packaged

set_option linter.unusedVariables false

variable [Fintype V] (E : V → V → Prop)
variable (hsymm : Symmetric E) (hirr : Irreflexive E)
variable (hinv : ∀ (g : G) (x y : V), E (g • x) (g • y) ↔ E x y)

include hsymm hirr hinv

/-- **(a)** Every invariant `S` is a union of orbits. -/
theorem part_a {S : Finset V} (h : Invariant G S) :
    S = S.biUnion (fun x => orbitF G x) ∧
      ∀ x : V, orbitF G x ⊆ S ∨ Disjoint (orbitF G x) S :=
  ⟨h.eq_biUnion_orbitF, h.orbitF_subset_or_disjoint⟩

/-- **(b)** An orbit `O` is a subset of some invariant independent `S` iff `O` is internally
independent. -/
theorem part_b {O : Finset V} (hO : IsOrbit G O) :
    (∃ S : Finset V, Invariant G S ∧ Independent E S ∧ O ⊆ S) ↔ Independent E O :=
  isOrbit_subset_invariant_independent_iff E hO

/-- **(c)** The maximum size of an invariant independent set equals the maximum, over families of
pairwise-compatible internally independent orbits, of the sum of their cardinalities — and both
maxima are attained. -/
theorem part_c :
    sSup (invIndepSizes G E) = sSup (orbitFamilySizes G E) ∧
      IsGreatest (invIndepSizes G E) (sSup (invIndepSizes G E)) ∧
      IsGreatest (orbitFamilySizes G E) (sSup (orbitFamilySizes G E)) :=
  ⟨sSup_invIndepSizes_eq E, isGreatest_sSup_invIndepSizes E, isGreatest_sSup_orbitFamilySizes E⟩

end Packaged

end OrbitReduction
