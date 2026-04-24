//! Phase IV PnS — End-to-end bi-4 cap assembly.
//!
//! Composes the PnS pipeline primitives into a single function
//! [`assemble_pns_bi4_cap`] that turns a boundary ring + pole into
//! `N` bi-4 sector patches:
//!
//! ```text
//! (Bi2TensorBorderRing, β, σ, pole) → N × [[Point3; 5]; 5]
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
//! | Nxt sector (j=4 edge)    | `(0..4, 4)`                                | 5     |
//! | Free (F_κ min)           | `(0, 3), (1, 2), (1, 3), (2, 2), (2, 3)`   | 5     |
//!
//! The j=4 row is SHARED with nxt sector's j=0 row (cap topology):
//! `P[k][i, 4] = P[(k+1) mod N][i, 0]`. In our representation each
//! sector reads its j=4 row via that alias (no new DOFs).
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
//! Total: `15 + 15·N` scalars. Problem is separable per coordinate
//! (x, y, z), so we solve three `(5 + 5·N)`-dimensional symmetric
//! PSD systems.
//!
//! # Eq (2) nxt_{2,1} closed form
//!
//! Eq (2) gives `P[k+1][2, 1] + P[k][2, 1] = b_k` where
//! ```text
//!   b_k = (2 − c)·P[k][2, 0] + c·P[k][3, 0]
//! ```
//! depends on `P[k][2, 0]` (central, linear in `a`) and
//! `P[k][3, 0]` (boundary, constant). For odd `n` the pair-sum
//! system has closed-form solution
//! ```text
//!   P[k][2, 1] = (1/2) · Σ_{j=0}^{n−1} (−1)^j · b_{(k+j) mod n}
//! ```
//! (same alternating-sum structure as [`super::pns_bi3_twists`]).
//! For even `n` the system is rank-`n − 1` and needs extra free
//! scalars; this first-cut supports **odd `n` only**.
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` §4 (eqs 2, 4, 5, 6), Fig. 8c.
//!
//! # Convention
//!
//! Per Fig. 8c: each sector is a bi-4 patch with local
//! `(i, j) ∈ [0, 4]²`:
//! - `i` runs along the seam with nxt (pole at `i=0`, outer at `i=4`)
//! - `j` runs transverse to the seam (`j=0` on seam with nxt)
//! - Cap's outer boundary is at `i=4` (per wip.md convention)

#![allow(dead_code)]
#![allow(clippy::needless_range_loop)]

use super::direct_linalg::solve_dense;
use super::pns_beta::BetaReparam;
use super::pns_boundary::{Bi2TensorBorderRing, compose_bi4_outer_layers};
use super::pns_central_quad::{QuadMapJet, central_bb_from_a_and_sigma};
use super::pns_gram::{bi4_flat_index, gram_bi4_f_kappa};
use super::pns_sigma::SigmaCenterJet;
use crate::geometry::Point3;

/// Inputs for the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleInput {
    /// Per-sector bi-2 tensor-border strips.
    pub boundary: Bi2TensorBorderRing,
    /// β reparameterisation for `n` sectors.
    pub beta: BetaReparam,
    /// Canonical σ centre jet.
    pub sigma: SigmaCenterJet,
    /// Pole (central point, eq 5).
    pub pole: Point3,
    /// F_κ order, `κ ∈ 1..=4`.
    pub kappa: usize,
}

/// Result of the bi-4 cap assembly.
#[derive(Debug, Clone)]
pub struct CapAssembleResult {
    /// Fitted quadratic map `a` (with `a.a_0 = pole`).
    pub a: QuadMapJet,
    /// `N` bi-4 sector patches.
    pub patches: Vec<[[Point3; 5]; 5]>,
}

/// Cap assembly error cases.
#[derive(Debug, thiserror::Error)]
pub enum CapAssembleError {
    #[error("PnS cap assemble: n must be in [3, 12], got {0}")]
    InvalidValence(usize),
    #[error("PnS cap assemble: even n not yet supported (got {0})")]
    EvenValenceUnsupported(usize),
    #[error("PnS cap assemble: kappa must be 1..=4, got {0}")]
    InvalidKappa(usize),
    #[error("PnS cap assemble: F_κ linear system is singular")]
    SingularSystem,
}

/// The 5 "free" CPs per sector, in fixed order for DOF indexing.
/// These are the deep-interior CPs not pinned by boundary, central,
/// closed-form, or shared-edge aliasing.
const FREE_CPS: [(usize, usize); 5] = [(0, 3), (1, 2), (1, 3), (2, 2), (2, 3)];

/// Index within [`FREE_CPS`] for a given `(i, j)`. Returns `None`
/// if the CP is not free.
fn free_cp_index(i: usize, j: usize) -> Option<usize> {
    FREE_CPS.iter().position(|&(a, b)| a == i && b == j)
}

/// Total scalar DOFs per coordinate: 5 (a's 5 Point3 coefs) + 5·n
/// (free CPs per sector).
fn dof_count(n: usize) -> usize {
    5 + 5 * n
}

/// DOF index for a's `a_x` (0), `a_y` (1), `a_xx` (2), `a_xy` (3),
/// `a_yy` (4). 0..=4.
fn dof_a(a_idx: usize) -> usize {
    a_idx
}

/// DOF index for sector `k`'s free CP #`f` (0..5 per `FREE_CPS`).
fn dof_free(k: usize, f: usize) -> usize {
    5 + 5 * k + f
}

/// Assemble a PnS bi-4 cap from boundary + pole + σ jet.
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

    // Precompute per-sector σ^k and boundary (row0, row1).
    //
    // σ convention: per the clean-room clarification (KK-Peters
    // reference implementation author), adjacent sector patches do
    // NOT relate by rotation — the "common recipe" applies the
    // SAME σ formulas in each sector. So σ^k = σ^0 for all k here.
    // The earlier image-space `rotate_sigma_jet` per sector was
    // removed once the paper's construction was clarified.
    let mut sigma_k = Vec::with_capacity(n);
    let mut row0 = Vec::with_capacity(n);
    let mut row1 = Vec::with_capacity(n);
    for k in 0..n {
        sigma_k.push(input.sigma);
        let (r0, r1) = compose_bi4_outer_layers(&input.boundary.sectors[k], &input.beta);
        row0.push(r0);
        row1.push(r1);
    }

    // Build the coefficient map: for each sector k and each CP
    // (i, j), express it as a linear combination of DOFs plus a
    // constant Point3.
    let coef_map = build_coef_map(input, &sigma_k, &row0, &row1);

    // Assemble the DOF-space Gram matrix and RHS, solve per coord,
    // reconstruct patches.
    let dofs = solve_dof_system(&coef_map, input.kappa, n)?;
    let patches = reconstruct_patches(&coef_map, &dofs, n);
    let a = reconstruct_a(&dofs, input.pole);

    Ok(CapAssembleResult { a, patches })
}

/// Linear representation of a single BB-coef: `cp = Σ_d coef[d]·dof[d] + constant`.
#[derive(Debug, Clone)]
struct CpRepr {
    /// Scalar coefficients (same across x/y/z by separability).
    /// Length `dof_count(n)`.
    coef: Vec<f64>,
    /// Constant Point3 (contribution from pole + boundary).
    constant: Point3,
}

impl CpRepr {
    fn zeros(n: usize) -> Self {
        Self {
            coef: vec![0.0; dof_count(n)],
            constant: Point3::new(0.0, 0.0, 0.0),
        }
    }
}

/// Full coefficient map: `sector_cps[k][5·i + j]` gives the
/// representation for `P[k][i, j]`.
#[derive(Debug, Clone)]
struct CoefMap {
    n: usize,
    sector_cps: Vec<Vec<CpRepr>>,
}

fn build_coef_map(
    input: &CapAssembleInput,
    sigma_k: &[SigmaCenterJet],
    row0: &[[Point3; 5]],
    row1: &[[Point3; 5]],
) -> CoefMap {
    let n = input.boundary.n();
    let mut sector_cps: Vec<Vec<CpRepr>> = (0..n).map(|_| vec![CpRepr::zeros(n); 25]).collect();

    // Pass 1: fill in pole, central (a∘σ^k), boundary, free CPs.
    for k in 0..n {
        for i in 0..5 {
            for j in 0..5 {
                let repr = compute_non_shared_cp(input, sigma_k, row0, row1, k, i, j);
                if let Some(r) = repr {
                    sector_cps[k][bi4_flat_index(i, j)] = r;
                }
            }
        }
    }

    // Pass 2: j=4 row aliases to NEXT sector's j=0 row (cap
    // topology: cur's v=4 seam is shared with nxt's v=0 seam).
    for k in 0..n {
        let nxt = (k + 1) % n;
        for i in 0..5 {
            let src = sector_cps[nxt][bi4_flat_index(i, 0)].clone();
            sector_cps[k][bi4_flat_index(i, 4)] = src;
        }
    }

    // Pass 3: closed-form p_{2, 1} via pair-sum over all sectors.
    // P[k][2, 1] = (1/2) Σ_j (-1)^j · b_{(k+j) mod n}
    // where b_m = (2-c)·P[m][2, 0] + c·P[m][3, 0].
    let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();
    let two_minus_c = 2.0 - c;
    for k in 0..n {
        let mut combined = CpRepr::zeros(n);
        for j in 0..n {
            let sign = if j % 2 == 0 { 0.5 } else { -0.5 };
            let m = (k + j) % n;
            let p20_m = &sector_cps[m][bi4_flat_index(2, 0)];
            let p30_m = &sector_cps[m][bi4_flat_index(3, 0)];
            for d in 0..combined.coef.len() {
                combined.coef[d] += sign * (two_minus_c * p20_m.coef[d] + c * p30_m.coef[d]);
            }
            combined.constant.x += sign * (two_minus_c * p20_m.constant.x + c * p30_m.constant.x);
            combined.constant.y += sign * (two_minus_c * p20_m.constant.y + c * p30_m.constant.y);
            combined.constant.z += sign * (two_minus_c * p20_m.constant.z + c * p30_m.constant.z);
        }
        sector_cps[k][bi4_flat_index(2, 1)] = combined;
    }

    // Also update the j=4 aliases for (2, 4) now that (2, 0) is
    // finalised. Actually (2, 0) doesn't depend on p_{2, 1}, so
    // pass 2 copy is still valid. No re-sync needed.

    CoefMap { n, sector_cps }
}

/// Compute the representation for `P[k][i, j]` EXCLUDING shared
/// j=4 row entries and closed-form (2, 1). Returns `None` for
/// those (handled by later passes).
fn compute_non_shared_cp(
    input: &CapAssembleInput,
    sigma_k: &[SigmaCenterJet],
    row0: &[[Point3; 5]],
    row1: &[[Point3; 5]],
    k: usize,
    i: usize,
    j: usize,
) -> Option<CpRepr> {
    let n = input.boundary.n();

    // Deferred: j=4 row (pass 2), closed-form (2, 1) (pass 3).
    if j == 4 {
        return None;
    }
    if (i, j) == (2, 1) {
        return None;
    }

    let mut repr = CpRepr::zeros(n);

    // Pole.
    if (i, j) == (0, 0) {
        repr.constant = input.pole;
        return Some(repr);
    }

    // Central (a∘σ^k): i+j ≤ 2 and (i, j) ≠ (0, 0).
    if i + j <= 2 {
        central_cp_repr(&mut repr, &sigma_k[k], input.pole, i, j);
        return Some(repr);
    }

    // Boundary i=3 or i=4.
    if i == 4 {
        repr.constant = row0[k][j];
        return Some(repr);
    }
    if i == 3 {
        repr.constant = Point3::new(
            row0[k][j].x + 0.25 * row1[k][j].x,
            row0[k][j].y + 0.25 * row1[k][j].y,
            row0[k][j].z + 0.25 * row1[k][j].z,
        );
        return Some(repr);
    }

    // Free CP.
    if let Some(f) = free_cp_index(i, j) {
        repr.coef[dof_free(k, f)] = 1.0;
        return Some(repr);
    }

    // Any remaining CP (i ∈ {1, 2}, j ∈ {2, 3} not in FREE_CPS,
    // or i=0 at j=3 already in FREE_CPS, etc.) should have been
    // covered above. Should be unreachable.
    unreachable!(
        "build_coef_map: unexpected CP ({i}, {j}) in sector {k} (i+j={}, free? {})",
        i + j,
        free_cp_index(i, j).is_some()
    );
}

/// Fill `repr` for a central CP `(i, j)` with `i+j ≤ 2` (excluding
/// the pole). Coefficients depend on `sigma_k` via the chain rule
/// in [`central_bb_from_a_and_sigma`]; constants = pole.
fn central_cp_repr(repr: &mut CpRepr, sigma: &SigmaCenterJet, pole: Point3, i: usize, j: usize) {
    // Use unit-basis evaluation: for each a-DOF, set that
    // component of `a` to (1, 0, 0), others zero, evaluate
    // central_bb_from_a_and_sigma; the resulting Point3.x = coef.
    //
    // (By separability the coef on a's x-component equals the coef
    // on y- and z-components, so we need only the x-component.)
    let zero = Point3::new(0.0, 0.0, 0.0);
    let basis = Point3::new(1.0, 0.0, 0.0);

    // Baseline (a = zero, a_0 = pole): central CP = pole.
    // Already handled as constant.
    repr.constant = pole;

    // a_x
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: basis,
            a_y: zero,
            a_xx: zero,
            a_xy: zero,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(0)] = cp.x - pole.x;
    }
    // a_y
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: basis,
            a_xx: zero,
            a_xy: zero,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(1)] = cp.x - pole.x;
    }
    // a_xx
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: zero,
            a_xx: basis,
            a_xy: zero,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(2)] = cp.x - pole.x;
    }
    // a_xy
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: zero,
            a_xx: zero,
            a_xy: basis,
            a_yy: zero,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(3)] = cp.x - pole.x;
    }
    // a_yy
    {
        let a = QuadMapJet {
            a_0: pole,
            a_x: zero,
            a_y: zero,
            a_xx: zero,
            a_xy: zero,
            a_yy: basis,
        };
        let cp = cp_at_ij(central_bb_from_a_and_sigma(&a, sigma), i, j);
        repr.coef[dof_a(4)] = cp.x - pole.x;
    }
}

/// Extract central BB coef `(i, j)` from a `CentralSixBB`.
fn cp_at_ij(cbb: super::pns_central_quad::CentralSixBB, i: usize, j: usize) -> Point3 {
    match (i, j) {
        (0, 0) => cbb.p00,
        (1, 0) => cbb.p10,
        (0, 1) => cbb.p01,
        (2, 0) => cbb.p20,
        (1, 1) => cbb.p11,
        (0, 2) => cbb.p02,
        _ => panic!("cp_at_ij: not a central CP ({i}, {j})"),
    }
}

/// Solve the joint F_κ system in DOF space, per coordinate.
fn solve_dof_system(
    coef_map: &CoefMap,
    kappa: usize,
    n: usize,
) -> Result<Vec<[f64; 3]>, CapAssembleError> {
    let gram = gram_bi4_f_kappa(kappa);
    let n_dofs = dof_count(n);

    // Build G_dof = Σ_k M_k^T · G · M_k  (size n_dofs × n_dofs).
    // M_k is the (25) × n_dofs matrix of scalar coefficients for
    // sector k (using coef_map.sector_cps[k][idx].coef).
    //
    // Also build b_dof = Σ_k M_k^T · G · c_k (three components:
    // b_dof_x, b_dof_y, b_dof_z) where c_k is sector k's constants.
    let mut g_dof = vec![vec![0.0_f64; n_dofs]; n_dofs];
    let mut b_dof_x = vec![0.0_f64; n_dofs];
    let mut b_dof_y = vec![0.0_f64; n_dofs];
    let mut b_dof_z = vec![0.0_f64; n_dofs];

    for k in 0..coef_map.n {
        let cps = &coef_map.sector_cps[k];
        for a in 0..25 {
            for b in 0..25 {
                let gab = gram[a][b];
                if gab == 0.0 {
                    continue;
                }
                let ma = &cps[a];
                let mb = &cps[b];
                // M_k^T · G · M_k contribution:
                //   for each (d1, d2): G_dof[d1][d2] += ma.coef[d1] · gab · mb.coef[d2]
                for d1 in 0..n_dofs {
                    let ca = ma.coef[d1];
                    if ca == 0.0 {
                        continue;
                    }
                    for d2 in 0..n_dofs {
                        let cb = mb.coef[d2];
                        if cb == 0.0 {
                            continue;
                        }
                        g_dof[d1][d2] += ca * gab * cb;
                    }
                }
                // M_k^T · G · c_k contribution:
                //   b_dof[d] += ma.coef[d] · gab · mb.constant.?
                // Accumulate for each coord.
                for d in 0..n_dofs {
                    let ca = ma.coef[d];
                    if ca == 0.0 {
                        continue;
                    }
                    b_dof_x[d] += ca * gab * mb.constant.x;
                    b_dof_y[d] += ca * gab * mb.constant.y;
                    b_dof_z[d] += ca * gab * mb.constant.z;
                }
            }
        }
    }

    // Solve G_dof · x = -b_dof per coord. Gram is symmetric PSD;
    // if the free-DOF subspace has no kernel, solve_dense works.
    let mut dofs = vec![[0.0_f64; 3]; n_dofs];

    let xs =
        solve_dense(g_dof.clone(), negate(&b_dof_x)).ok_or(CapAssembleError::SingularSystem)?;
    let ys =
        solve_dense(g_dof.clone(), negate(&b_dof_y)).ok_or(CapAssembleError::SingularSystem)?;
    let zs = solve_dense(g_dof, negate(&b_dof_z)).ok_or(CapAssembleError::SingularSystem)?;
    for d in 0..n_dofs {
        dofs[d] = [xs[d], ys[d], zs[d]];
    }
    Ok(dofs)
}

fn negate(v: &[f64]) -> Vec<f64> {
    v.iter().map(|x| -x).collect()
}

/// Reconstruct the cap patches from solved DOF values.
fn reconstruct_patches(coef_map: &CoefMap, dofs: &[[f64; 3]], n: usize) -> Vec<[[Point3; 5]; 5]> {
    let mut patches = vec![[[Point3::new(0.0, 0.0, 0.0); 5]; 5]; n];
    for k in 0..n {
        for i in 0..5 {
            for j in 0..5 {
                let repr = &coef_map.sector_cps[k][bi4_flat_index(i, j)];
                let mut p = repr.constant;
                for d in 0..repr.coef.len() {
                    let c = repr.coef[d];
                    if c == 0.0 {
                        continue;
                    }
                    p.x += c * dofs[d][0];
                    p.y += c * dofs[d][1];
                    p.z += c * dofs[d][2];
                }
                patches[k][i][j] = p;
            }
        }
    }
    patches
}

/// Reconstruct `a` from the DOF solution (first 5 DOFs, per coord).
fn reconstruct_a(dofs: &[[f64; 3]], pole: Point3) -> QuadMapJet {
    let dof_to_p3 = |d: usize| Point3::new(dofs[d][0], dofs[d][1], dofs[d][2]);
    QuadMapJet {
        a_0: pole,
        a_x: dof_to_p3(dof_a(0)),
        a_y: dof_to_p3(dof_a(1)),
        a_xx: dof_to_p3(dof_a(2)),
        a_xy: dof_to_p3(dof_a(3)),
        a_yy: dof_to_p3(dof_a(4)),
    }
}

#[cfg(test)]
#[path = "pns_cap_assemble_tests.rs"]
mod tests;
