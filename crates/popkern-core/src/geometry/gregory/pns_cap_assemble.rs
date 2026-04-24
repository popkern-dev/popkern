//! Phase IV PnS — End-to-end bi-4 cap assembly.
//!
//! Composes the PnS pipeline primitives into a single function
//! [`assemble_pns_bi4_cap`] that turns a boundary ring + pole into
//! `N` bi-4 sector patches:
//!
//! ```text
//! (Bi2TensorBorderRing, β, pole) → N × [[Point3; 5]; 5]
//! ```
//!
//! # Architecture
//!
//! Each sector `k` has a `5 × 5` grid of BB-coefs `P[k][i][j]`.
//! Their sources:
//!
//! | Source                   | CPs per sector                             | Count |
//! |--------------------------|--------------------------------------------|-------|
//! | Pole (constant)          | `(0, 0)`                                   | 1     |
//! | Central `a∘σ^k`          | `(1, 0), (0, 1), (2, 0), (1, 1), (0, 2)`   | 5     |
//! | Eq (2) pair-sum closed   | `(2, 1)`                                   | 1     |
//! | Boundary `b∘β`           | `(3, 0..4), (4, 0..4)`                     | 10    |
//! | Prev sector (j=4 edge)   | `(0..4, 4)` (duplicates)                   | 5     |
//! | Free (F_κ min)           | `(0, 3), (1, 2), (1, 3), (2, 2), (2, 3)`   | 5     |
//! | **Total**                |                                            | 27¹   |
//!
//! ¹ 27 = 25 CPs per sector + 2 duplicated via j=4 sharing with prev;
//!   the shared CPs collapse.
//!
//! The cap has `1 + 19·N` unique CPs overall (pole + per-seam
//! non-pole + per-sector interior).
//!
//! # DOFs and the F_κ solve
//!
//! The truly free variables, determined jointly by the F_κ
//! minimisation (paper eq 6), are:
//!
//! - `a`'s 5 Point3 components `(a_x, a_y, a_xx, a_xy, a_yy)` — 15
//!   scalars (shared across sectors; `a_0` pinned to pole)
//! - `5 × N` free interior CPs per sector — `15·N` scalars
//!
//! Total: `15 + 15·N` free scalars. The problem is separable per
//! coordinate (x, y, z), so we solve three `(5 + 5·N)`-dimensional
//! symmetric PSD systems.
//!
//! # Eq (2) nxt_{2,1} closed form
//!
//! Eq (2) gives `P[k+1][2, 1] + P[k][2, 1] = b_k` where `b_k`
//! depends on `P[k][2, 0]` (central, linear in `a`) and `P[k][3, 0]`
//! (boundary, constant). For odd `n` the pair-sum system has a
//! unique closed-form solution (alternating sum) and the
//! `(2, 1)` CPs are eliminated from the F_κ DOFs. For even `n`
//! the system is rank-`n − 1` and requires one extra free scalar;
//! **only odd `n ∈ {3, 5, 7, 9, ...}` is supported in this first
//! cut**; even-`n` support is a follow-up.
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` §4 (eq 2, 4, 5, 6), Fig. 8c (cap BB
//! labels, shared seam convention).
//!
//! # Convention: indexing within a sector patch
//!
//! Following paper's Fig. 8c: each sector is a bi-4 patch with
//! local `(i, j) ∈ [0, 4]²`:
//! - `i` runs ALONG the shared seam with nxt, from pole (`i=0`) to
//!   outer end (`i=4`)
//! - `j` runs TRANSVERSE to the seam, into the sector interior
//!   (`j=0` on seam with nxt, `j=4` on seam with prev)
//! - Cap's outer boundary is at `i=4` (per wip.md and
//!   [`super::pns_boundary`])
//!
//! Shared-edge identification across sectors (per clean-room
//! clarification): `P[k][i, 0] = P[k+1][i, 4]` (i.e., sector k's
//! seam-with-nxt at j=0 equals sector (k+1)'s seam-with-prev at
//! j=4). Both sides refer to the SAME physical Point3.

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::pns_beta::BetaReparam;
use super::pns_boundary::Bi2TensorBorderRing;
use super::pns_central_quad::QuadMapJet;
use super::pns_sigma::SigmaCenterJet;
use crate::geometry::Point3;

/// Inputs for the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleInput {
    /// Per-sector bi-2 tensor-border strips (cap's outer boundary
    /// and first transversal derivative).
    pub boundary: Bi2TensorBorderRing,
    /// β reparameterisation for `n` sectors.
    pub beta: BetaReparam,
    /// Canonical σ centre jet (sector-0). Sector `k`'s σ is
    /// obtained by rotating this jet in image space by `2π·k/n`.
    pub sigma: SigmaCenterJet,
    /// Pole = central point (paper's eq 5).
    pub pole: Point3,
    /// F_κ quality functional order. Typical: `κ = 3` for convex
    /// caps (balanced smoothness vs fit tightness).
    pub kappa: usize,
}

/// Result of the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleResult {
    /// Fitted quadratic map `a: R² → R³` (paper's eq 6
    /// minimisation). `a.a_0 = pole`.
    pub a: QuadMapJet,
    /// `N` bi-4 sector patches. `patches[k][i][j]` is `P[k][i][j]`.
    pub patches: Vec<[[Point3; 5]; 5]>,
}

/// Cap assembly error cases.
#[derive(Debug, thiserror::Error)]
pub enum CapAssembleError {
    /// Valence `n < 3` or `n > 12` (upper bound is practical
    /// — numerical conditioning of the F_κ solve degrades at
    /// very high n).
    #[error("PnS cap assemble: n must be in [3, 12], got {0}")]
    InvalidValence(usize),
    /// Even `n` is not supported in this first-cut implementation
    /// (eq (2) pair-sum system for p_{2,1} is rank-deficient for
    /// even `n`; needs extra free-scalar handling).
    #[error("PnS cap assemble: even n not yet supported (got {0})")]
    EvenValenceUnsupported(usize),
    /// `κ` out of valid range `1..=4`.
    #[error("PnS cap assemble: kappa must be 1..=4, got {0}")]
    InvalidKappa(usize),
    /// Boundary ring size mismatches `n`.
    #[error("PnS cap assemble: boundary ring has {actual} sectors, expected {expected}")]
    BoundaryValenceMismatch { expected: usize, actual: usize },
    /// The `(5 + 5·n)`-dimensional F_κ system is singular or
    /// ill-conditioned (should not happen in practice but reported
    /// defensively).
    #[error("PnS cap assemble: F_κ linear system is singular")]
    SingularSystem,
}

/// Assemble a PnS bi-4 cap from boundary + pole + σ jet.
///
/// Performs the joint F_κ minimisation over `a` and per-sector
/// free CPs, subject to the eq (2) pair-sum closed-form on
/// `P[k][2, 1]`.
///
/// # Errors
///
/// See [`CapAssembleError`] for failure modes.
pub fn assemble_pns_bi4_cap(
    input: &CapAssembleInput,
) -> Result<CapAssembleResult, CapAssembleError> {
    let n = input.boundary.n();
    if !(3..=12).contains(&n) {
        return Err(CapAssembleError::InvalidValence(n));
    }
    if n.is_multiple_of(2) {
        return Err(CapAssembleError::EvenValenceUnsupported(n));
    }
    if !(1..=4).contains(&input.kappa) {
        return Err(CapAssembleError::InvalidKappa(input.kappa));
    }
    // Scaffold — full solve forthcoming in subsequent commits.
    // For now: return an N-sector cap with all zeros except the
    // pinned boundary rows, as a structural placeholder.
    let placeholder = build_scaffold(input);
    Ok(placeholder)
}

/// Build a first-cut placeholder result: pole pinned at `(0, 0)`
/// per sector, boundary rows (`i = 3, 4`) filled from `b∘β`, all
/// other CPs zero. Used until the full F_κ solve lands.
fn build_scaffold(input: &CapAssembleInput) -> CapAssembleResult {
    use super::pns_boundary::compose_bi4_outer_layers;

    let n = input.boundary.n();
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut patches = vec![[[zero; 5]; 5]; n];

    for k in 0..n {
        let strip = &input.boundary.sectors[k];
        let (row0, row1) = compose_bi4_outer_layers(strip, &input.beta);
        // Cap convention: pns_boundary's row0 = cap's i=4 row
        // (outer boundary, indexed by j=0..4); row1 gives i=3 via
        // p_{3, j} = p_{4, j} + row1[j]/4 (see mod docs).
        for j in 0..5 {
            patches[k][4][j] = row0[j];
            patches[k][3][j] = Point3::new(
                row0[j].x + 0.25 * row1[j].x,
                row0[j].y + 0.25 * row1[j].y,
                row0[j].z + 0.25 * row1[j].z,
            );
        }
        // Pole at (0, 0)
        patches[k][0][0] = input.pole;
    }

    // Placeholder `a`: linear part zero, Hessian zero. Real F_κ
    // solve will populate.
    let a = QuadMapJet {
        a_0: input.pole,
        a_x: zero,
        a_y: zero,
        a_xx: zero,
        a_xy: zero,
        a_yy: zero,
    };

    CapAssembleResult { a, patches }
}

#[cfg(test)]
#[path = "pns_cap_assemble_tests.rs"]
mod tests;
