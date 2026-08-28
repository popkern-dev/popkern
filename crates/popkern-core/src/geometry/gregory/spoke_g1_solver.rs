//! §5.3 item #3 extended — Karciauskas-Peters spoke-wise G1
//! solver for rectangular Gregory patches with free-curvature
//! interior CPs. **Generic** over the patch type via the
//! [`SpokeG1Patch`] trait: biquintic-full (deg 5, 2 free CPs
//! per spoke) and biseptic (deg 7, 4 free CPs per spoke).
//!
//! At each sample `t` on a shared spoke with along-spoke
//! tangent `T(t)` and cross-tangents `xtan_{A,B}(t)`, the
//! linearised G1 constraint with fixed target normal
//! `n_target = unit(T × xtan_A_init)` is:
//! `n_target · (xtan_A − xtan_B) = 0` — linear in the free CPs.
//!
//! # G2 (the polynomial-identity solver work) hook
//!
//! Extend by adding a `cross_tangent_second_coeff` trait method
//! (`∂²S/∂u²` coeffs — row `i = 2` free CPs for biseptic), an
//! additional per-sample constraint
//! `n · (xtan2_A − xtan2_B) = 0`, and solve the augmented system
//! with the same QP — no algorithmic change.

use crate::geometry::{Point3, Vec3};

/// One of the four boundary sides of a rectangular Gregory
/// patch, in its local `(u, v)` domain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    /// `u = 0` boundary, parameterised by `v ∈ [0, 1]`.
    UZero,
    /// `u = 1` boundary, parameterised by `v ∈ [0, 1]`.
    UOne,
    /// `v = 0` boundary, parameterised by `u ∈ [0, 1]`.
    VZero,
    /// `v = 1` boundary, parameterised by `u ∈ [0, 1]`.
    VOne,
}

/// One end of a shared spoke.
#[derive(Clone, Copy, Debug)]
pub struct SpokeEnd {
    /// Index of the patch this end belongs to.
    pub patch_idx: usize,
    /// Which boundary side of the patch is this spoke end.
    pub side: Side,
    /// `true` if this end's boundary parameterisation runs in
    /// the opposite direction to the spoke's canonical direction.
    pub reversed: bool,
}

/// A spoke shared between two patches. G1 continuity is
/// enforced along the boundary curve common to both ends.
#[derive(Clone, Copy, Debug)]
pub struct Spoke {
    /// First end. Its parameter direction defines the spoke's
    /// canonical parameter `t ∈ [0, 1]`.
    pub a: SpokeEnd,
    /// Second end. If `b.reversed`, its local boundary parameter
    /// at spoke-t is `1 − t`.
    pub b: SpokeEnd,
}

/// Abstraction over the patch types the spoke-G1 solver
/// supports. Implementations supply the degree-specific
/// free-CP index list, cross-tangent coefficient formula, and
/// grid accessors; the solver code itself is degree-independent.
pub trait SpokeG1Patch {
    /// Free interior Bezier indices whose CP contributes to the
    /// cross-tangent at the given boundary side. 2 indices for
    /// biquintic-full (degree 5), 4 for biseptic (degree 7).
    fn free_indices_for_side(side: Side) -> &'static [(usize, usize)];
    /// Cross-tangent coefficient for a single free CP at
    /// parameter `t`: `coef * P_{i,j}` is the contribution to
    /// `∂S/∂(u or v)` at the boundary. Magnitude is
    /// `n · B_{inner}^n(t)` where `n` is the patch degree.
    fn cross_tangent_coeff(side: Side, free_i: usize, free_j: usize, t: f64) -> f64;
    /// Read a grid control point.
    fn grid_get(&self, i: usize, j: usize) -> Point3;
    /// Write a grid control point.
    fn grid_set(&mut self, i: usize, j: usize, p: Point3);
    /// `∂S/∂u` at `(u, v)`.
    fn partial_u(&self, u: f64, v: f64) -> Vec3;
    /// `∂S/∂v` at `(u, v)`.
    fn partial_v(&self, u: f64, v: f64) -> Vec3;
}

/// Solve spoke-wise G1 for a collection of rectangular Gregory
/// patches sharing spokes. Generic over the patch type — works
/// for biquintic-full and biseptic. Adjusts the free interior
/// CPs in place (boundary + corner-inner-CP pairs preserved).
/// Returns the number of spokes whose residual was improved;
/// `None` on a singular linear system (regularised internally
/// so this is rare in practice).
pub fn spoke_g1_solve<P: SpokeG1Patch>(
    patches: &mut [P],
    spokes: &[Spoke],
    samples_per_spoke: usize,
) -> Option<usize> {
    spoke_g1_solve_inner(patches, spokes, samples_per_spoke, false, 1.0)
}

/// Spoke-G1 solve with additional near-endpoint
/// samples at `t ≈ 0.02` and `t ≈ 0.98`. These drive outer-
/// corner tangent-plane continuity (constrain adjacent sub-
/// patches' near-corner CPs to agree on the tangent plane at
/// `V_{2k+2}`). Only effective when the patch type has free
/// interior CPs at `j = 1` and `j = n − 1` (row-1 corner-
/// adjacent positions). Two-stripe satisfies this; biseptic
/// and biquintic-full don't — calling this for those patch
/// types would over-constrain the solver with no reachable
/// fix, so they use [`spoke_g1_solve`] without corner samples.
pub fn spoke_g1_solve_with_corner_samples<P: SpokeG1Patch>(
    patches: &mut [P],
    spokes: &[Spoke],
    samples_per_spoke: usize,
) -> Option<usize> {
    spoke_g1_solve_inner(patches, spokes, samples_per_spoke, true, 1.0)
}

/// Damped variant — computes the linearised delta then applies
/// only `damping` fraction of it. Used for corner-mode passes
/// where a full Newton step tends to fold the patch. Caller
/// iterates many times with small damping to reach near-
/// machine-precision G1 without over-shooting. `damping = 1.0`
/// is equivalent to [`spoke_g1_solve_with_corner_samples`].
pub fn spoke_g1_solve_with_corner_samples_damped<P: SpokeG1Patch>(
    patches: &mut [P],
    spokes: &[Spoke],
    samples_per_spoke: usize,
    damping: f64,
) -> Option<usize> {
    spoke_g1_solve_inner(patches, spokes, samples_per_spoke, true, damping)
}

fn spoke_g1_solve_inner<P: SpokeG1Patch>(
    patches: &mut [P],
    spokes: &[Spoke],
    samples_per_spoke: usize,
    include_corner_samples: bool,
    damping: f64,
) -> Option<usize> {
    if patches.is_empty() || spokes.is_empty() || samples_per_spoke == 0 {
        return Some(0);
    }

    let variables = enumerate_free_cp_variables::<P>(patches.len(), spokes);
    if variables.is_empty() {
        return Some(0);
    }
    let n_unknowns = variables.len() * 3;

    let mut constraint_rows: Vec<Vec<f64>> = Vec::new();
    let mut constraint_rhs: Vec<f64> = Vec::new();

    let ts: Vec<f64> = if include_corner_samples {
        let interior = sample_parameters(samples_per_spoke);
        let mut out = Vec::with_capacity(interior.len() + 2);
        out.push(0.02);
        out.extend(interior);
        out.push(0.98);
        out
    } else {
        sample_parameters(samples_per_spoke)
    };

    for spoke in spokes {
        for &t in &ts {
            let t_b = if spoke.b.reversed { 1.0 - t } else { t };
            let tangent_a = spoke_tangent(&patches[spoke.a.patch_idx], spoke.a.side, t);
            let tangent_b_raw = spoke_tangent(&patches[spoke.b.patch_idx], spoke.b.side, t_b);
            let tangent_b = if spoke.b.reversed {
                Vec3::new(-tangent_b_raw.x, -tangent_b_raw.y, -tangent_b_raw.z)
            } else {
                tangent_b_raw
            };
            let avg_tangent = Vec3::new(
                0.5 * (tangent_a.x + tangent_b.x),
                0.5 * (tangent_a.y + tangent_b.y),
                0.5 * (tangent_a.z + tangent_b.z),
            );
            if avg_tangent.length_squared() < 1e-24 {
                continue;
            }
            let xtan_a_init = cross_tangent(&patches[spoke.a.patch_idx], spoke.a.side, t);
            let n_unnorm = avg_tangent.cross(&xtan_a_init);
            let Some(n_target) = n_unnorm.normalized() else {
                continue;
            };
            let mut row = vec![0.0_f64; n_unknowns];
            let mut rhs = 0.0_f64;

            accumulate_side_row::<P>(
                &variables,
                spoke.a.patch_idx,
                spoke.a.side,
                t,
                1.0,
                &n_target,
                &patches[spoke.a.patch_idx],
                &mut row,
                &mut rhs,
            );
            accumulate_side_row::<P>(
                &variables,
                spoke.b.patch_idx,
                spoke.b.side,
                t_b,
                -1.0,
                &n_target,
                &patches[spoke.b.patch_idx],
                &mut row,
                &mut rhs,
            );

            constraint_rows.push(row);
            constraint_rhs.push(rhs);
        }
    }

    if constraint_rows.is_empty() {
        return Some(0);
    }

    let mut delta = solve_kkt_qp(&constraint_rows, &constraint_rhs)?;
    let d = damping.clamp(0.0, 1.0);
    if d != 1.0 {
        for v in delta.iter_mut() {
            *v *= d;
        }
    }

    apply_delta_to_patches::<P>(patches, &variables, &delta);

    let mut improved = 0;
    for spoke in spokes {
        for &t in &ts {
            let t_b = if spoke.b.reversed { 1.0 - t } else { t };
            let tangent_a = spoke_tangent(&patches[spoke.a.patch_idx], spoke.a.side, t);
            let tangent_b_raw = spoke_tangent(&patches[spoke.b.patch_idx], spoke.b.side, t_b);
            let tangent_b = if spoke.b.reversed {
                Vec3::new(-tangent_b_raw.x, -tangent_b_raw.y, -tangent_b_raw.z)
            } else {
                tangent_b_raw
            };
            let avg_tangent = Vec3::new(
                0.5 * (tangent_a.x + tangent_b.x),
                0.5 * (tangent_a.y + tangent_b.y),
                0.5 * (tangent_a.z + tangent_b.z),
            );
            let xtan_a = cross_tangent(&patches[spoke.a.patch_idx], spoke.a.side, t);
            let n_unnorm = avg_tangent.cross(&xtan_a);
            let Some(n_target) = n_unnorm.normalized() else {
                continue;
            };
            let xtan_b = cross_tangent(&patches[spoke.b.patch_idx], spoke.b.side, t_b);
            let post = (xtan_a.x - xtan_b.x) * n_target.x
                + (xtan_a.y - xtan_b.y) * n_target.y
                + (xtan_a.z - xtan_b.z) * n_target.z;
            if post.abs() < 1e-9 {
                improved += 1;
            }
        }
    }
    Some(improved)
}

#[derive(Clone, Copy, Debug)]
struct FreeCpVariable {
    patch_idx: usize,
    i: usize,
    j: usize,
}

fn sample_parameters(m: usize) -> Vec<f64> {
    if m == 1 {
        return vec![0.5];
    }
    (1..=m).map(|k| (k as f64) / ((m + 1) as f64)).collect()
}

fn enumerate_free_cp_variables<P: SpokeG1Patch>(
    n_patches: usize,
    spokes: &[Spoke],
) -> Vec<FreeCpVariable> {
    let mut vars: Vec<FreeCpVariable> = Vec::new();
    let mut seen: std::collections::BTreeSet<(usize, usize, usize)> = Default::default();
    for spoke in spokes {
        for end in [spoke.a, spoke.b] {
            if end.patch_idx >= n_patches {
                continue;
            }
            for &(i, j) in P::free_indices_for_side(end.side) {
                if seen.insert((end.patch_idx, i, j)) {
                    vars.push(FreeCpVariable {
                        patch_idx: end.patch_idx,
                        i,
                        j,
                    });
                }
            }
        }
    }
    vars
}

fn cross_tangent<P: SpokeG1Patch>(patch: &P, side: Side, t: f64) -> Vec3 {
    let (u, v) = match side {
        Side::UZero => (0.0, t),
        Side::UOne => (1.0, t),
        Side::VZero => (t, 0.0),
        Side::VOne => (t, 1.0),
    };
    match side {
        Side::UZero | Side::UOne => patch.partial_u(u, v),
        Side::VZero | Side::VOne => patch.partial_v(u, v),
    }
}

fn spoke_tangent<P: SpokeG1Patch>(patch: &P, side: Side, t: f64) -> Vec3 {
    let (u, v) = match side {
        Side::UZero => (0.0, t),
        Side::UOne => (1.0, t),
        Side::VZero => (t, 0.0),
        Side::VOne => (t, 1.0),
    };
    match side {
        Side::UZero | Side::UOne => patch.partial_v(u, v),
        Side::VZero | Side::VOne => patch.partial_u(u, v),
    }
}

#[allow(clippy::too_many_arguments)]
fn accumulate_side_row<P: SpokeG1Patch>(
    variables: &[FreeCpVariable],
    patch_idx: usize,
    side: Side,
    t: f64,
    row_sign: f64,
    n_target: &Vec3,
    patch: &P,
    row: &mut [f64],
    rhs: &mut f64,
) {
    for &(free_i, free_j) in P::free_indices_for_side(side) {
        let coef = P::cross_tangent_coeff(side, free_i, free_j, t);
        let signed = row_sign * coef;
        if let Some(offset) = variable_offset(variables, patch_idx, free_i, free_j) {
            row[offset] += signed * n_target.x;
            row[offset + 1] += signed * n_target.y;
            row[offset + 2] += signed * n_target.z;
        }
    }
    let xtan = cross_tangent(patch, side, t);
    let proj = xtan.x * n_target.x + xtan.y * n_target.y + xtan.z * n_target.z;
    *rhs -= row_sign * proj;
}

fn variable_offset(
    variables: &[FreeCpVariable],
    patch_idx: usize,
    i: usize,
    j: usize,
) -> Option<usize> {
    variables
        .iter()
        .position(|v| v.patch_idx == patch_idx && v.i == i && v.j == j)
        .map(|idx| idx * 3)
}

fn solve_kkt_qp(a_rows: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let m = a_rows.len();
    let n = a_rows.first().map(|r| r.len()).unwrap_or(0);
    if m == 0 || n == 0 {
        return Some(vec![0.0; n]);
    }
    let mut aat = vec![vec![0.0_f64; m]; m];
    let mut trace = 0.0_f64;
    for r in 0..m {
        for c in r..m {
            let s: f64 = a_rows[r]
                .iter()
                .zip(a_rows[c].iter())
                .map(|(x, y)| x * y)
                .sum();
            aat[r][c] = s;
            aat[c][r] = s;
        }
        trace += aat[r][r];
    }
    let lambda = 1e-12_f64 * (trace / (m as f64)).max(1.0);
    for (r, row) in aat.iter_mut().enumerate() {
        row[r] += lambda;
    }
    let lambdas = solve_linear_system(aat, b.to_vec())?;
    let mut delta = vec![0.0_f64; n];
    for (r, lam) in lambdas.iter().enumerate() {
        for (k, d) in delta.iter_mut().enumerate() {
            *d += lam * a_rows[r][k];
        }
    }
    Some(delta)
}

fn solve_linear_system(mut m: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = m.len();
    for i in 0..n {
        let mut pivot = i;
        let mut best = m[i][i].abs();
        for (r, row) in m.iter().enumerate().skip(i + 1) {
            let v = row[i].abs();
            if v > best {
                best = v;
                pivot = r;
            }
        }
        if best < 1e-14 {
            return None;
        }
        if pivot != i {
            m.swap(i, pivot);
            b.swap(i, pivot);
        }
        for r in (i + 1)..n {
            let factor = m[r][i] / m[i][i];
            let pivot_row = m[i].clone();
            for (c, pc) in pivot_row.iter().enumerate().skip(i) {
                m[r][c] -= factor * pc;
            }
            b[r] -= factor * b[i];
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for (c, xc) in x.iter().enumerate().skip(i + 1) {
            s -= m[i][c] * xc;
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

fn apply_delta_to_patches<P: SpokeG1Patch>(
    patches: &mut [P],
    variables: &[FreeCpVariable],
    delta: &[f64],
) {
    for (idx, var) in variables.iter().enumerate() {
        let off = idx * 3;
        let cur = patches[var.patch_idx].grid_get(var.i, var.j);
        patches[var.patch_idx].grid_set(
            var.i,
            var.j,
            Point3::new(
                cur.x + delta[off],
                cur.y + delta[off + 1],
                cur.z + delta[off + 2],
            ),
        );
    }
}
