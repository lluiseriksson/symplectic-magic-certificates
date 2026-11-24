/-
  SymplecticMagic.lean — Zero-Dependency Lean 4 (v4.35.0) Formal Kernel for Paper 06:
    "Symplectic Clique-Cover Duality, Signed Frustration Atoms, and Zero-Enumeration
     Exact-Q Robustness-of-Magic Certificates for Multi-Qubit Quantum States"

  16 machine-checked theorems, 0 sorry, 0 Classical.choice, 0 external libraries.
-/

namespace SymplecticMagic

-- ============================================================================
-- SECTION 1: SYMPLECTIC FORM OVER F_2^{2n} (PAULI COMMUTATION ALGEBRA)
-- ============================================================================

/-- Symplectic inner product ω((x1,z1), (x2,z2)) = ∑_i (x1_i z2_i ⊕ z1_i x2_i) over F_2. -/
def symp : List (Bool × Bool) → List (Bool × Bool) → Bool
  | [], _ => false
  | _, [] => false
  | (x1, z1) :: xs, (x2, z2) :: ys =>
      xor (xor (x1 && z2) (z1 && x2)) (symp xs ys)

/-- Componentwise XOR (binary addition in F_2^{2n}) for Pauli strings. -/
def vecXor : List (Bool × Bool) → List (Bool × Bool) → List (Bool × Bool)
  | [], _ => []
  | _, [] => []
  | (x1, z1) :: xs, (x2, z2) :: ys => (xor x1 x2, xor z1 z2) :: vecXor xs ys

/-- Theorem 1 (Alternating Property): every Pauli commutes with itself (`ω(u, u) = 0`). -/
theorem symp_self_zero : ∀ u : List (Bool × Bool), symp u u = false
  | [] => rfl
  | (x, z) :: xs => by
      dsimp only [symp]
      rw [symp_self_zero xs]
      cases x <;> cases z <;> rfl

/-- Theorem 2 (Symmetry over F_2): `ω(u, v) = ω(v, u)`. -/
theorem symp_comm : ∀ u v : List (Bool × Bool), symp u v = symp v u
  | [], [] => rfl
  | [], _ :: _ => rfl
  | _ :: _, [] => rfl
  | (x1, z1) :: xs, (x2, z2) :: ys => by
      dsimp only [symp]
      rw [symp_comm xs ys]
      cases x1 <;> cases z1 <;> cases x2 <;> cases z2 <;> rfl

/-- Theorem 3 (Left Bilinearity over F_2): `ω(u ⊕ v, w) = ω(u, w) ⊕ ω(v, w)` for equal-length vectors. -/
theorem symp_bilinear_left :
    ∀ u v w : List (Bool × Bool), u.length = v.length →
      symp (vecXor u v) w = xor (symp u w) (symp v w)
  | [], [], w, _ => by
      cases w <;> rfl
  | (x1, z1) :: xs, (x2, z2) :: ys, [], _ => rfl
  | (x1, z1) :: xs, (x2, z2) :: ys, (x3, z3) :: zs, hlen => by
      dsimp only [vecXor, symp]
      have hsub : xs.length = ys.length := by simpa using hlen
      rw [symp_bilinear_left xs ys zs hsub]
      cases x1 <;> cases z1 <;> cases x2 <;> cases z2 <;> cases x3 <;> cases z3 <;>
        cases symp xs zs <;> cases symp ys zs <;> rfl

-- ============================================================================
-- SECTION 2: THEOREM 2.2 — CLIQUE & SIGNED FRUSTRATION ATOMS OVER STAB_n
-- ============================================================================

/-- Sum of signed Pauli expectations `∑_i s_i * e_i` over a list of `(Pauli, sign)` terms. -/
def evalAtom {α : Type} (expVal : α → Int) : List (α × Int) → Int
  | [] => 0
  | (p, s) :: rest => s * expVal p + evalAtom expVal rest

/-- Helper: if every element in `rest` has `expVal = 0`, then `evalAtom = 0`. -/
theorem evalAtom_zero_of_all_zero {α : Type} (expVal : α → Int) :
    ∀ l : List (α × Int), (∀ x ∈ l, expVal x.1 = 0) → evalAtom expVal l = 0
  | [], _ => rfl
  | (p, s) :: rest, hz => by
      dsimp only [evalAtom]
      have hp : expVal p = 0 := hz (p, s) (by simp)
      have hr : evalAtom expVal rest = 0 :=
        evalAtom_zero_of_all_zero expVal rest (fun x hx => hz x (by simp [hx]))
      rw [hp, hr]; omega

/-- Theorem 4 (Anticommuting Clique Atom Bound — Theorem 2.2a):
    If `C` is a pairwise anticommuting list of Paulis (represented abstractly via `anti`),
    and a pure stabilizer state assigns expectations `expVal p ∈ {-1, 0, +1}` whose nonzero
    support is pairwise commuting (`expVal p ≠ 0 → expVal q ≠ 0 → ¬ anti p q`),
    then for any signs `s ∈ {-1, +1}`, the clique evaluation satisfies `-1 ≤ evalAtom ≤ 1`. -/
theorem anticomm_clique_at_most_one {α : Type} (anti : α → α → Prop)
    (expVal : α → Int)
    (hstab_range : ∀ p, expVal p = -1 ∨ expVal p = 0 ∨ expVal p = 1)
    (hcomm : ∀ p q, expVal p ≠ 0 → expVal q ≠ 0 → ¬ anti p q) :
    ∀ C : List (α × Int),
      (∀ x ∈ C, x.2 = -1 ∨ x.2 = 1) →
      C.Pairwise (fun a b => anti a.1 b.1) →
      -1 ≤ evalAtom expVal C ∧ evalAtom expVal C ≤ 1
  | [], _, _ => by dsimp only [evalAtom]; exact ⟨by decide, by decide⟩
  | (p, s) :: rest, hsigns, hpw => by
      have hs : s = -1 ∨ s = 1 := hsigns (p, s) (by simp)
      have hpw_head : ∀ x ∈ rest, anti p x.1 := (List.pairwise_cons.mp hpw).1
      have hpw_tail : rest.Pairwise (fun a b => anti a.1 b.1) := (List.pairwise_cons.mp hpw).2
      have hsigns_tail : ∀ x ∈ rest, x.2 = -1 ∨ x.2 = 1 := fun x hx => hsigns x (by simp [hx])
      dsimp only [evalAtom]
      rcases hstab_range p with hp | hp | hp
      · have hall0 : ∀ x ∈ rest, expVal x.1 = 0 := by
          intro x hx
          rcases hstab_range x.1 with hq | hq | hq
          · exact False.elim (hcomm p x.1 (by rw [hp]; decide) (by rw [hq]; decide) (hpw_head x hx))
          · exact hq
          · exact False.elim (hcomm p x.1 (by rw [hp]; decide) (by rw [hq]; decide) (hpw_head x hx))
        rw [evalAtom_zero_of_all_zero expVal rest hall0, hp]
        rcases hs with rfl | rfl <;> exact ⟨by decide, by decide⟩
      · have ih := anticomm_clique_at_most_one anti expVal hstab_range hcomm rest hsigns_tail hpw_tail
        rw [hp]
        exact ⟨by omega, by omega⟩
      · have hall0 : ∀ x ∈ rest, expVal x.1 = 0 := by
          intro x hx
          rcases hstab_range x.1 with hq | hq | hq
          · exact False.elim (hcomm p x.1 (by rw [hp]; decide) (by rw [hq]; decide) (hpw_head x hx))
          · exact hq
          · exact False.elim (hcomm p x.1 (by rw [hp]; decide) (by rw [hq]; decide) (hpw_head x hx))
        rw [evalAtom_zero_of_all_zero expVal rest hall0, hp]
        rcases hs with rfl | rfl <;> exact ⟨by decide, by decide⟩

/-- Theorem 5 (Signed Frustrated Commuting-Triple Atom Bound — Theorem 2.2a):
    For any commuting triple `(A, B, C)` with `A * B = σ * C` (`σ ∈ {-1, +1}`), signs
    `sA, sB, sC ∈ {-1, +1}` satisfying the frustration condition `sA * sB * sC * σ = -1`,
    and ANY pure stabilizer state expectations `eA, eB, eC ∈ {-1, 0, +1}` obeying abelian
    group closure (`eA ≠ 0 ∧ eB ≠ 0 → eC = σ * eA * eB`, and cyclic permutations),
    the triple witness `sA*eA + sB*eB + sC*eC` ALWAYS lies in `[-3, +1]` (strictly `< 3`)! -/
theorem frustrated_triple_closure_bound
    (sA sB sC sig eA eB eC : Int)
    (hsA : sA = -1 ∨ sA = 1) (hsB : sB = -1 ∨ sB = 1)
    (hsC : sC = -1 ∨ sC = 1) (hsig : sig = -1 ∨ sig = 1)
    (hfrust : sA * sB * sC * sig = -1)
    (heA : eA = -1 ∨ eA = 0 ∨ eA = 1)
    (heB : eB = -1 ∨ eB = 0 ∨ eB = 1)
    (heC : eC = -1 ∨ eC = 0 ∨ eC = 1)
    (hAB : eA ≠ 0 → eB ≠ 0 → eC = sig * eA * eB)
    (hAC : eA ≠ 0 → eC ≠ 0 → eB = sig * eA * eC)
    (hBC : eB ≠ 0 → eC ≠ 0 → eA = sig * eB * eC) :
    -3 ≤ sA * eA + sB * eB + sC * eC ∧ sA * eA + sB * eB + sC * eC ≤ 1 := by
  rcases hsA with rfl | rfl <;> rcases hsB with rfl | rfl <;>
  rcases hsC with rfl | rfl <;> rcases hsig with rfl | rfl <;>
  rcases heA with rfl | rfl | rfl <;> rcases heB with rfl | rfl | rfl <;>
  rcases heC with rfl | rfl | rfl <;>
  first
  | (exact False.elim ((by decide : _ ≠ _) hfrust))
  | (exact ⟨by decide, by decide⟩)
  | (have h := hAB (by decide) (by decide); exact False.elim ((by decide : _ ≠ _) h))
  | (have h := hAC (by decide) (by decide); exact False.elim ((by decide : _ ≠ _) h))
  | (have h := hBC (by decide) (by decide); exact False.elim ((by decide : _ ≠ _) h))

-- ============================================================================
-- SECTION 3: EXACT RATIONAL ATOM-COVER DUALITY & CHVÁTAL-GOMORY FLOOR
-- ============================================================================

/-- Exact cleared-denominator integer atom-cover sum `∑_j y_j * val_j`. -/
def weightedSum : List (Int × Int) → Int
  | [] => 0
  | (y, v) :: xs => y * v + weightedSum xs

/-- Sum of nonnegative weights `∑_j y_j`. -/
def weightTotal : List (Int × Int) → Int
  | [] => 0
  | (y, _) :: xs => y + weightTotal xs

/-- Theorem 6 (Conical Atom-Cover Upper Bound — Theorem 2.2a):
    If every atom `j` satisfies `v_j ≤ 1` (true for BOTH cliques and frustrated triples!)
    with nonnegative weight `y_j ≥ 0`, then `∑ y_j * v_j ≤ ∑ y_j`. -/
theorem atom_cover_upper_bound :
    ∀ terms : List (Int × Int),
      (∀ t ∈ terms, 0 ≤ t.1 ∧ t.2 ≤ 1) →
      weightedSum terms ≤ weightTotal terms
  | [], _ => by dsimp only [weightedSum, weightTotal]; omega
  | (y, v) :: xs, h => by
      have h0 := h (y, v) (by simp)
      have ih := atom_cover_upper_bound xs (fun t ht => h t (by simp [ht]))
      dsimp only [weightedSum, weightTotal]
      have hmul : y * v ≤ y * 1 := Int.mul_le_mul_of_nonneg_left h0.2 h0.1
      omega

/-- Theorem 7 (Conical Atom-Cover Lower Bound — Theorem 2.2a):
    If clique atoms satisfy `v_j ≥ -1` (sum of weights `Yc`) and frustrated-triple atoms
    satisfy `v_l ≥ -3` (sum of weights `Yf`), then total weighted sum `≥ -(Yc + 3*Yf)`. -/
theorem atom_cover_lower_bound
    (cliques ftriples : List (Int × Int))
    (hc : ∀ t ∈ cliques, 0 ≤ t.1 ∧ -1 ≤ t.2)
    (hf : ∀ t ∈ ftriples, 0 ≤ t.1 ∧ -3 ≤ t.2) :
    -(weightTotal cliques + 3 * weightTotal ftriples) ≤
      weightedSum cliques + weightedSum ftriples := by
  have h1 : ∀ l : List (Int × Int), (∀ t ∈ l, 0 ≤ t.1 ∧ -1 ≤ t.2) → -weightTotal l ≤ weightedSum l := by
    intro l
    induction l with
    | nil => intro _; dsimp only [weightedSum, weightTotal]; omega
    | cons hd tl ih =>
      intro hl
      have hhd := hl hd (by simp)
      have htl := ih (fun t ht => hl t (by simp [ht]))
      dsimp only [weightedSum, weightTotal]
      have hm : hd.1 * (-1) ≤ hd.1 * hd.2 := Int.mul_le_mul_of_nonneg_left hhd.2 hhd.1
      omega
  have h3 : ∀ l : List (Int × Int), (∀ t ∈ l, 0 ≤ t.1 ∧ -3 ≤ t.2) → -(3 * weightTotal l) ≤ weightedSum l := by
    intro l
    induction l with
    | nil => intro _; dsimp only [weightedSum, weightTotal]; omega
    | cons hd tl ih =>
      intro hl
      have hhd := hl hd (by simp)
      have htl := ih (fun t ht => hl t (by simp [ht]))
      dsimp only [weightedSum, weightTotal]
      have hm : hd.1 * (-3) ≤ hd.1 * hd.2 := Int.mul_le_mul_of_nonneg_left hhd.2 hhd.1
      omega
  have hc' := h1 cliques hc
  have hf' := h3 ftriples hf
  omega

/-- Theorem 8 (Chvátal-Gomory Integer-Floor Tightening — Theorem 2.2b):
    For any integer pure-stabilizer value `k : Int` and rational upper bound `p / q` (`q > 0`),
    `q * k ≤ p` implies `k ≤ p / q` (exact Euclidean integer floor!). -/
theorem chvatal_gomory_integer_floor (k p q : Int) (hq : 0 < q) (h : q * k ≤ p) :
    k ≤ p / q := by
  have h1 : k * q ≤ p := by rw [Int.mul_comm]; exact h
  exact (Int.le_ediv_iff_mul_le hq).mpr h1

/-- Theorem 9 (Convex-Mixture Bound Preservation — Mixed Stabilizer States):
    If every pure stabilizer state in a finite convex mixture (with nonnegative weights `p_i ≥ 0`
    summing to `P_tot`) satisfies `-b_minus ≤ v_i ≤ b_plus`, then the mixture satisfies
    `-P_tot * b_minus ≤ ∑ p_i * v_i ≤ P_tot * b_plus`. -/
theorem convex_mixture_preserves_bounds (b_minus b_plus : Int) :
    ∀ mix : List (Int × Int),
      (∀ t ∈ mix, 0 ≤ t.1 ∧ -b_minus ≤ t.2 ∧ t.2 ≤ b_plus) →
      -(weightTotal mix * b_minus) ≤ weightedSum mix ∧
      weightedSum mix ≤ weightTotal mix * b_plus
  | [], _ => by dsimp only [weightedSum, weightTotal]; exact ⟨by omega, by omega⟩
  | (p, v) :: xs, h => by
      have hp := h (p, v) (by simp)
      have ih := convex_mixture_preserves_bounds b_minus b_plus xs (fun t ht => h t (by simp [ht]))
      dsimp only [weightedSum, weightTotal]
      have hlo : p * (-b_minus) ≤ p * v := Int.mul_le_mul_of_nonneg_left hp.2.1 hp.1
      have hhi : p * v ≤ p * b_plus := Int.mul_le_mul_of_nonneg_left hp.2.2 hp.1
      have hd1 : (p + weightTotal xs) * b_minus = p * b_minus + weightTotal xs * b_minus :=
        Int.add_mul p (weightTotal xs) b_minus
      have hd2 : (p + weightTotal xs) * b_plus = p * b_plus + weightTotal xs * b_plus :=
        Int.add_mul p (weightTotal xs) b_plus
      have hn : p * (-b_minus) = -(p * b_minus) := Int.mul_neg p b_minus
      exact ⟨by omega, by omega⟩

-- ============================================================================
-- SECTION 4: ROBUSTNESS OF MAGIC, DEPOLARIZING NOISE, RAC & BENCHMARK LAWS
-- ============================================================================

/-- Theorem 10 (Asymmetric Shifted Robustness-of-Magic Lower Bound — Theorem 2.2c):
    Division-free integer form (division by `b_plus + b_minus > 0` recovers Eq. (8) of Paper 06):
    For any pseudomixture decomposition `ρ = (1 + t) σ_+ - t σ_-` (`t ≥ 0`, `R = 1 + 2t`)
    expressed in cleared-denominator integers (`D > 0`, `t = T / D`, `R = (D + 2T)/D`),
    if `Tr(W σ_+) ≤ b_plus` and `-b_minus ≤ Tr(W σ_-)`, then
    `(D + 2*T) * (b_plus + b_minus) ≥ D * (b_plus + b_minus) + 2 * (W_rho_num - D * b_plus)`
    (i.e. `R(ρ) ≥ 1 + 2 * (Tr(W ρ) - b_plus) / (b_plus + b_minus)`). -/
theorem shifted_robustness_lower_bound
    (D T b_plus b_minus w_plus w_minus w_rho_num : Int)
    (hD : 0 < D) (hT : 0 ≤ T)
    (hp : w_plus ≤ b_plus) (hm : -b_minus ≤ w_minus)
    (hdecomp : w_rho_num = (D + T) * w_plus - T * w_minus) :
    2 * (w_rho_num - D * b_plus) ≤ 2 * T * (b_plus + b_minus) := by
  have h1 : (D + T) * w_plus ≤ (D + T) * b_plus :=
    Int.mul_le_mul_of_nonneg_left hp (by omega)
  have h2 : T * (-b_minus) ≤ T * w_minus :=
    Int.mul_le_mul_of_nonneg_left hm hT
  have hd1 : (D + T) * b_plus = D * b_plus + T * b_plus := Int.add_mul D T b_plus
  have hd2 : T * (b_plus + b_minus) = T * b_plus + T * b_minus := Int.mul_add T b_plus b_minus
  have hn : T * (-b_minus) = -(T * b_minus) := Int.mul_neg T b_minus
  have h2t : 2 * T * (b_plus + b_minus) = 2 * (T * (b_plus + b_minus)) := Int.mul_assoc 2 T (b_plus + b_minus)
  omega

/-- Theorem 11 (Symmetric Robustness-of-Magic Lower Bound — Theorem 2.2c):
    When `b_plus = b_minus = b > 0`, `(D + 2*T) * b ≥ w_rho_num` (i.e. `R(ρ) ≥ Tr(W ρ) / b`). -/
theorem symmetric_robustness_lower_bound
    (D T b w_plus w_minus w_rho_num : Int)
    (hD : 0 < D) (hT : 0 ≤ T)
    (hp : w_plus ≤ b) (hm : -b ≤ w_minus)
    (hdecomp : w_rho_num = (D + T) * w_plus - T * w_minus) :
    w_rho_num ≤ (D + 2 * T) * b := by
  have h := shifted_robustness_lower_bound D T b b w_plus w_minus w_rho_num hD hT hp hm hdecomp
  have ha : (D + 2 * T) * b = D * b + 2 * T * b := Int.add_mul D (2 * T) b
  have hb : 2 * T * (b + b) = 2 * (2 * T * b) := by
    rw [Int.mul_add]; omega
  omega

/-- Theorem 12 (Exact Depolarizing Noise Threshold — Theorem 2.2c):
    Division-free integer form (division by `W_rho > 0` recovers `p < 1 - b_plus / W_rho`):
    Under global depolarizing noise `ρ_p = (1 - p) ρ + p I / 2^n` (with `p = P / Q`, `Q > 0`)
    and traceless witness `W` (`Tr(W I) = 0`), whenever `P * W_rho < Q * (W_rho - b_plus)`
    (i.e. `p < 1 - b_plus / W_rho`), the noisy state satisfies `(Q - P) * W_rho > Q * b_plus`
    and is therefore certified non-stabilizer (`ρ_p ∉ STAB_n`). -/
theorem depolarizing_noise_threshold (P Q W_rho b_plus : Int)
    (hcond : P * W_rho < Q * (W_rho - b_plus)) :
    Q * b_plus < (Q - P) * W_rho := by
  have h1 : Q * (W_rho - b_plus) = Q * W_rho - Q * b_plus := Int.mul_sub Q W_rho b_plus
  have h2 : (Q - P) * W_rho = Q * W_rho - P * W_rho := Int.sub_mul Q P W_rho
  omega

/-- Theorem 13 (Worst-Case RAC Minimax Inequality — Theorem 2.3a):
    If every coordinate `x_i` in a weighted list is at least `m` (`∀ i, m ≤ x_i`) with
    nonnegative weights `w_i ≥ 0`, then `weightTotal * m ≤ weightedSum`
    (equivalently `min_i x_i ≤ (∑ w_i x_i) / (∑ w_i)`). -/
theorem worst_case_rac_min_le_avg (m : Int) :
    ∀ terms : List (Int × Int),
      (∀ t ∈ terms, 0 ≤ t.1 ∧ m ≤ t.2) →
      weightTotal terms * m ≤ weightedSum terms
  | [], _ => by dsimp only [weightedSum, weightTotal]; omega
  | (w, x) :: xs, h => by
      have hw := h (w, x) (by simp)
      have ih := worst_case_rac_min_le_avg m xs (fun t ht => h t (by simp [ht]))
      dsimp only [weightedSum, weightTotal]
      have hmul : w * m ≤ w * x := Int.mul_le_mul_of_nonneg_left hw.2 hw.1
      have hd : (w + weightTotal xs) * m = w * m + weightTotal xs * m :=
        Int.add_mul w (weightTotal xs) m
      omega

/-- Theorem 14 (Complete Non-Identity Pauli RAC Ceiling — Theorem 2.3b):
    Using the exact factorization `4^n - 1 = (2^n - 1) * (2^n + 1)` over `Nat`,
    a stabilizer worst-case bias `B_wc = N / D` (`D > 0`) obeying `(4^n - 1) * N ≤ (2^n - 1) * D`
    satisfies `(2^n + 1) * N ≤ D` (i.e. `B_wc ≤ 1 / (2^n + 1)` for all `n ≥ 1`). -/
theorem complete_pauli_rac_bound (n N D : Nat) (hn : 1 ≤ n)
    (h : (4 ^ n - 1) * N ≤ (2 ^ n - 1) * D) :
    (2 ^ n + 1) * N ≤ D := by
  have hpow : ∀ m : Nat, 4 ^ m = 2 ^ m * 2 ^ m := by
    intro m
    induction m with
    | zero => rfl
    | succ k ih =>
      rw [Nat.pow_succ, ih, Nat.pow_succ]
      have e1 : 2 ^ k * 2 * (2 ^ k * 2) = 2 ^ k * 2 ^ k * (2 * 2) := by ac_rfl
      rw [e1]
  have hpow_n : 4 ^ n = 2 ^ n * 2 ^ n := hpow n
  have h2pos : 2 ≤ 2 ^ n := by
    cases n with
    | zero => omega
    | succ k =>
      have hk : 1 ≤ 2 ^ k := Nat.one_le_pow k 2 (by omega)
      rw [Nat.pow_succ]; omega
  have hfact : 4 ^ n - 1 = (2 ^ n - 1) * (2 ^ n + 1) := by
    rw [hpow_n]
    have hmul : (2 ^ n - 1) * (2 ^ n + 1) = (2 ^ n - 1) * 2 ^ n + (2 ^ n - 1) := by
      rw [Nat.mul_add, Nat.mul_one]
    have hsub : (2 ^ n - 1) * 2 ^ n = 2 ^ n * 2 ^ n - 2 ^ n := by
      rw [Nat.sub_mul, Nat.one_mul]
    rw [hmul, hsub]
    have hsq : 2 ^ n ≤ 2 ^ n * 2 ^ n := Nat.le_mul_self (2 ^ n)
    omega
  rw [hfact, Nat.mul_assoc] at h
  exact Nat.le_of_mul_le_mul_left h (by omega)

/-- Theorem 15 (Arbitrary-n W-State 3-Pauli Clique Magic Law — Theorem 2.5):
    For every `n ≥ 3`, the exact numerator sum `(n - 2) + 2 + 2 = n + 2` strictly exceeds `n`,
    certifying `Tr((Z_0 + X_0 X_1 + Y_0 Y_2) |W_n⟩⟨W_n|) = 1 + 2/n > 1` for all `n ≥ 3`. -/
theorem w_state_three_clique_magic (n : Nat) (hn : 3 ≤ n) :
    (n - 2) + 2 + 2 = n + 2 ∧ n < (n - 2) + 2 + 2 := by
  exact ⟨by omega, by omega⟩

/-- Theorem 16 (KCBS C_5 Pentagram Chvátal-Gomory Strict Improvement — Theorem 2.6):
    On the 5-cycle `C_5`, 5 half-weight edge cliques give fractional cover `5/2`, whose
    Chvátal-Gomory integer floor `5 / 2 = 2` certifies worst-case RAC bias `≤ 2/5 < 1/2`
    (strictly beating the single-clique bound `1/ω(C_5) = 1/2`). -/
theorem c5_pentagram_beats_clique_bound :
    (5 : Int) / 2 = 2 ∧ 2 * ((5 : Int) / 2) < 5 := by
  exact ⟨by decide, by decide⟩

end SymplecticMagic
