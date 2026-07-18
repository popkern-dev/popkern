# Catmull-Clark characteristic ring — independent verification scripts

First-principles Python rebuild of the CC characteristic map machinery,
used during the 2026-07-02 KP20 recovery review to verify PR08
(Peters & Reif, *Subdivision Surfaces*, Springer 2008, §5–§6.1) and
Stam 1998 (`docs/papers/sig98.pdf`) against mesh connectivity, and to
cross-check the two sources against each other.

These scripts are the reference oracle for the Rust
characteristic-ring substrate (see
`docs/architecture/kp20-paper-faithful-recovery-plan.md`, WP1).
Requires numpy (and sympy for `sym_eig.py`).

- `verify_pr08.py` — builds the 13n-coefficient CC subdivision operator
  from mesh connectivity (stencils applied to decoded q-label grid
  coordinates) and checks the printed PR08 eq (6.4)/(6.5) DFT blocks.
  **Found erratum E1 (unpublished): PR08 p.112 row k′=4 first/third
  entries are swapped — correct row is `q₁δ, q₁+q₂w, q₂ | q₁, q₁w, 0, 0`.**
- `verify_final.py` / `sym_eig.py` — closed-form/symbolic subdominant
  eigenvector; **found erratum E2 (unpublished): PR08 p.113
  v̂_im[4], v̂_im[6] = ±84λ²(8λ+1), not the printed ±140λ²(8λ−1)
  (they agree only at n = 4).** Also applies the published erratum for
  the f(1,1,0) normalization matrix ((3,3) entry is v̂¹¹).
- `verify_stam.py` — Stam extended subdivision matrices S, S11/S12,
  S21/S22 from connectivity, the closed-form frequency-1 characteristic
  net, the outer-7 extension solve `(S12 − λI)x = −S11·ring`, one Ā
  subdivision, and the three exact picking vectors q¹/q²/q³ (with the
  N=3 fold q¹₁=2).
- `final_xcheck.py` — proves Stam's three picked 4×4 grids equal the
  PR08 characteristic-ring patch grids (A/B/C ↔ Ω₃/Ω₂/Ω₁) up to one
  global complex constant, ≤ 4e−15 for n ∈ {3..8, 12}.
- `verify_eig.py` — eigenvalue/subdominance guards.

Known consumer-facing consequence: the kernel's
`kp20_bi6_cc_ring.rs` transcribed the misprinted PR08 row k′=4
(erratum E1), so every characteristic-ring value it produced for
q₄..q₁₃ is wrong at every valence. See the recovery plan, RC5.

## PK-134 resolution (2026-07-18): the printed µ0 table reproduced

The KP20 eq (23) µ0 table IS the eq (21) normalization of the
F3-solved paper tau. Proof lives in the Rust gate
`printed_mu0_table_is_the_eq21_normalized_solved_tau`
(`kp20_paper_tau_eq21_fit.rs`): with µ FIXED at the printed values, a
9-parameter sigma-equivariant domain 3-jet psi maps the solved tau's
cubic jet at o exactly onto the eq (22) template — residual
1e-10..1e-9 at every printed valence n in {3, 5..10}, with 11
equations of over-determination beyond the fitted gauge.

Falsified en route (do not re-derive):
- free-µ fits match the family at 2.9e-10 with µ0 ~ 8.9-9.6 — the
  eq (22) template family is CLOSED under domain 3-jets, so a free µ0
  is gauge, not data (this is why `verify_mu0.py`'s closed-form
  candidates could never work: µ0 is only pinned relative to the
  printed eq (21) chart convention);
- quadratic-psi + eq (23) couplings strands Gauss-Newton in a wrong
  basin (residual 0.46).
