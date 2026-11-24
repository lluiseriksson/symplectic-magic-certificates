# Symplectic Clique-Cover Duality, Signed Frustration Atoms, and Zero-Enumeration Exact-$\mathbb{Q}$ Robustness-of-Magic Certificates for Multi-Qubit Quantum States

[![Lean 4 v4.35.0](https://img.shields.io/badge/Lean_4-v4.35.0_(16%2F16_Theorems,_0_sorry)-blue)](lean/SymplecticMagic.lean)
[![Zero-Float Verifier](https://img.shields.io/badge/Exact_Verifier-Q_Rational_(0_floats,_9%2F9_certs)-brightgreen)](src/magic_certifier.py)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Complete formal verification artifact, zero-floating-point-trust exact-$\mathbb{Q}$ magic certificate suite (`SYMPLECTIC-MAGIC-CERT-1.0`), multi-qubit benchmarks up to $n = 32$ and $n = 50$ qubits, and publication manuscript advancing `arXiv:2611.26705` (Zamora, Nov 2026) and Howard & Campbell (PRL 2017).

## Main Results

1. **Symplectic Atom-Cover & Chvátal–Gomory Rounding Theorem (Theorem 2.2):**
   Any conical decomposition $W = \sum_j y_j C_j + \sum_\ell z_\ell F_\ell$ into anticommuting Pauli cliques ($C_j \in [-1, +1]$) and signed frustrated commuting triples ($F_\ell \in [-3, +1]$ via $\mathbb{F}_2$-subgroup closure) yields zero-enumeration stabilizer bounds $[-\hat{\beta}_-, +\hat{\beta}_+]$ (tightened by exact integer floors $\lfloor \beta_\pm \rfloor$ whenever $W$ has integer weights), certifying exact Robustness of Magic $\mathcal{R}(\rho) \ge 1 + 2(\operatorname{Tr}(W\rho) - \hat{\beta}_+)/(\hat{\beta}_+ + \hat{\beta}_-)$ and depolarizing noise tolerance $p_{\mathrm{cert}} = 1 - \hat{\beta}_+ / \operatorname{Tr}(W\rho)$ without ever enumerating $\mathrm{STAB}_n$.
2. **Worst-Case RAC Ceilings (Theorem 2.3) & Exponential Separation via Majorana Cliques $K_{2n+1}$ (Theorem 2.4):**
   On all non-identity Paulis $\mathcal{P}_n^*$, every $\sigma \in \mathrm{STAB}_n$ satisfies $B_{\mathrm{wc}}(\sigma) \le 1/(2^n + 1)$. For every $n \ge 1$, the rank-$2^{n-1}$ mixed state $\rho_{K_{2n+1}}$ has Pauli $\ell_1$-norm bound $(1+\sqrt{2n+1})/2^n \to 0$ (blind for all $n \ge 2$) and purity $2^{1-n} \to 0$, yet our $(2n+1)$-Pauli clique certificate proves $\mathcal{R}(\rho_{K_{2n+1}}) \ge \sqrt{2n+1}$ and $p_{\mathrm{cert}} = 1 - 1/\sqrt{2n+1}$ (reaching $\mathcal{R} \ge 8.062$ and $p_{\mathrm{cert}} = 87.6\%$ at $n = 32$ qubits in `< 2 ms`, bypassing $|\mathrm{STAB}_{32}| \approx 10^{168.95}$).
3. **Three-Pauli Stabilizer-Marginal Magic Theorem for Arbitrary $n$-Qubit $W$-States (Theorem 2.5):**
   Every 2-qubit marginal of $|W_n\rangle$ is an exact stabilizer mixture $\rho^{(2)}_{W_n} \in \mathrm{STAB}_2$, yet the 3-Pauli weight-$\le 2$ clique $\{Z_0, X_0 X_1, Y_0 Y_2\}$ certifies $\mathcal{R}(|W_n\rangle) \ge \mathcal{R}(\rho^{(3)}_{W_n}) \ge 1 + 2/n$ and $p_{\mathrm{cert}} = 2/(n+2)$ for every $n \ge 3$.
4. **Strict 3-Tier Hierarchy Separation (Theorem 2.6):**
   - **Tier 1 (Single Clique) $>$ Pauli $\ell_1$:** $K_{2n+1}$ ($\sqrt{2n+1}$ vs $\to 0$) and $|\mathrm{CCZ}\rangle$ (4-Pauli clique $W_{\mathrm{CCZ}} = XII + ZXI - ZZX + YXY$ gives $\mathcal{R} \ge 2$ vs 63-Pauli $\ell_1 = 15/8 = 1.875$).
   - **Tier 2 (Chvátal–Gomory Clique Cover) $>$ Tier 1:** 2-qubit KCBS $C_5$ pentagram ($\lfloor 5/2 \rfloor = 2 \implies B_{\mathrm{wc}} \le 2/5 < 1/2$, certifying depolarized $B_{\mathrm{wc}} = 9/20$ with $\mathcal{R} \ge 9/8$ where every single-clique bound fails).
   - **Tier 3 (Signed Frustration Atoms) $>$ All Unsigned Graph Bounds $\alpha_w(G_{\mathrm{anti}})$:** On `W_FT6 = (XX+YY+ZZ) + (XI+YI+ZI)` for $|H_3\rangle^{\otimes 2}$, unsigned $\alpha_w = 3 > 1+\sqrt{3} \approx 2.732$ fails completely, whereas 1 frustrated triple + 1 clique certifies $\beta_+ = 2, \beta_- = 4 \implies \mathcal{R} \ge (2+\sqrt{3})/3 \approx 1.244$ ($p < 26.8\%$).
5. **Zero-Dependency Lean 4 Kernel (`lean/SymplecticMagic.lean`, 16/16 Theorems):**
   Machine-checked on `leanprover/lean4:v4.35.0` with **0 `sorry`**, **0 `Classical.choice`**, and **0 external dependencies**.

## One-Command End-to-End Verification

```bash
python scripts/run_all_experiments.py
python scripts/make_figures.py
python scripts/verify_all.py
```

`scripts/verify_all.py` independently compiles the Lean 4 kernel, audits all 16 theorem axioms (`AuditAxioms.lean`), replays all 9 exact-$\mathbb{Q}$ certificates (`results/certificates/*.cert.json`), runs the 13-test pytest suite (including **28 live adversarial forgery, canonical-parser, and type-injection attacks**), and verifies all SHA-256 checksums in `SHA256SUMS.json`.
