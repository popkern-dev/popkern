//! KP09a § 6 low-valence special cases (Deferred #4).
//!
//! KP09a § 6 ("Practical matters & special cases") describes what
//! to do at valences `n ∈ {3, 4}` where the KP13 circulant system
//! loses full rank. The paper's two prescriptions are:
//!
//! * **`n = 4`** — the circulant retains solutions via
//!   **guided surfacing**: the two extra free DoFs (two DFT modes
//!   `j ∈ {1, 3}` have zero denominator) are set by minimising
//!   `F_5` over the guide surface. All KP09a § 5 constraints
//!   (star quadratic, C¹/C², shared edges) remain enforceable.
//!
//! * **`n = 3`** — the paper recommends a **structural change**:
//!   either (a) remove the valence-3 EV via one Loop subdivision
//!   step (which converts 3-valent ↔ 5-valent for tetrahedra,
//!   Figure 15), or (b) use a dedicated ternary-sector
//!   construction. Option (a) is preferred for most meshes.
//!
//! # What this module delivers
//!
//! The [`solve_g21_circulant`](super::kp13_bi5_circulant::solve_g21_circulant)
//! entry point now accepts `n ∈ {3, 4}` via **min-norm
//! pseudo-inverse**, which matches exactly what an LSQ-based
//! guided-surfacing solver would produce when the free DoFs are
//! unconstrained (no `F_5` preference). In other words:
//!
//! - The paper's guided-surfacing choice biases the free DoFs
//!   toward the guide-surface shape — minimising `F_5(g△ − guide)`.
//! - Our min-norm pseudo-inverse chooses the free DoFs to have
//!   zero `L²` norm — minimising `‖g^i_{21}‖²` in the null-space
//!   of the circulant.
//!
//! Both are valid, well-defined KP09a § 6-compatible resolutions
//! of the free DoFs. They agree whenever the guide-surface data
//! is zero on the null space (which is the common case when the
//! guide is built from the boundary data via [`crate::geometry::gregory::kp13_bi5_lsq_cc`]).
//!
//! For deterministic, regression-stable results this module
//! standardises on **min-norm pseudo-inverse** as the default
//! low-valence resolver. Callers wanting a non-default
//! `F_5`-guided choice of the free DoFs can solve a follow-up
//! quadratic programme over the two-dimensional null-space — the
//! same machinery used by [`crate::geometry::gregory::kp15_bi6_cap_f5_solve`]
//! at bi-6 polar valence.
//!
//! # Null-space structure (for reference)
//!
//! The KP09a eq-8 circulant `M · g_{21} = 2R` has rank
//! `n − k(n)`, where `k(n)` counts DFT modes `j ∈ 0..n` with
//! `2c + cos(jα) = 0`, `c = cos(2π/n)`, `α = 2π/n`.
//!
//! | `n` | `c`     | singular modes `j`     | `k(n)` |
//! |-----|---------|------------------------|--------|
//! | 3   | `−1/2`  | `{0}` (constant)       | 1      |
//! | 4   | `0`     | `{1, 3}` (alternating) | 2      |
//! | 5   | `cos 72°` | none                 | 0      |
//! | 6   | `1/2`   | `{3}` (±1 alternating) | 1      |
//! | 7+  | varies  | none                   | 0      |
//!
//! At `n = 4` the two singular modes correspond to the kernel
//! vectors `(1, 0, −1, 0)` and `(0, 1, 0, −1)` (DFT at `j = 1`
//! and `j = 3` for a 4-cyclic recurrence).
//!
//! # Further reading
//!
//! * KP09a § 6 — the paper text on special cases.
//! * [`super::kp13_bi5_circulant::solve_g21_circulant`] — unified
//!   pseudo-inverse solver covering all `n ≥ 3`.
//! * `docs/journal/polar-kp-math.md` §4.1 — singular-mode analysis
//!   for the KP13 circulant.

#![allow(dead_code)]

use popkern_math::PortableMath;
use std::f64::consts::PI;

/// Returns the list of singular DFT-mode indices `j ∈ 0..n` at
/// which the KP13 circulant denominator `2c + cos(jα)` vanishes
/// (to within tolerance `1e-12`). Equivalent to the dimension of
/// the circulant's null-space is `singular_modes(n).len()`.
///
/// Matches the table in this module's header — smoke-test fodder
/// for valence-specific invariants.
pub fn singular_modes(n: usize) -> Vec<usize> {
    assert!(n >= 3, "singular_modes: need n ≥ 3, got {n}");
    let c = (2.0 * PI / n as f64).pcos();
    let alpha = 2.0 * PI / n as f64;
    let mut out = Vec::new();
    for j in 0..n {
        let denom = 2.0 * c + (j as f64 * alpha).pcos();
        if denom.abs() < 1e-12 {
            out.push(j);
        }
    }
    out
}

/// Dimension of the KP09a eq-8 circulant's null-space at valence
/// `n`. Alias for `singular_modes(n).len()`. Useful when a caller
/// only cares about the rank deficiency.
pub fn null_space_dim(n: usize) -> usize {
    singular_modes(n).len()
}

#[cfg(test)]
#[path = "kp13_bi5_low_valence_tests.rs"]
mod tests;
