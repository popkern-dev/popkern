//! Bidegree-(9, 9) Bezier patch with Chiyokura-Kimura dual-twist
//! corners.
//!
//! Identical to [`RectangularGregoryBiNonic`] except the four
//! interior corner controls `(1,1)`, `(1,8)`, `(8,1)`, `(8,8)`
//! are evaluation-time rational blends of per-boundary candidate
//! pairs. A plain tensor product has ONE twist control per
//! corner demanded by BOTH boundary cross-derivative bands
//! through that corner; when the demands differ (curvature ×
//! control spacing makes them differ on real fillet corners), any
//! single value damages at least one boundary's band. The blend
//! gives each boundary its own candidate: along `v = 0` the
//! blended control reduces to `P_v` exactly (and the rational
//! term's `∂/∂v` contribution is killed by the `B_1^9(0) = 0`
//! factor), so each boundary band realizes its own derived
//! coefficient with zero interference from the other.
//!
//! Design: `docs/decisions/2026-06-11-kp20-binonic-dual-twist-corners.md`.
//! Precedent at degrees 3/5/7: `rectangular.rs`,
//! `rectangular_deg5.rs`, `rectangular_biquintic_full.rs`.

use crate::geometry::gregory::RectangularGregoryBiNonic;
use crate::geometry::gregory::bernstein_binonic::{
    bernstein_binonic_all, bernstein_binonic_deriv_all,
};
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Bi-nonic patch with rational dual-twist corner controls.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularGregoryBiNonicDualTwist {
    /// 10×10 control-point grid, row-major (`grid[i][j]`, `i`
    /// indexes `u`). Entries `(1,1)`, `(1,8)`, `(8,1)`, `(8,8)`
    /// are placeholders — evaluation replaces them with the
    /// rational blend from [`inner_cps`](Self::inner_cps).
    pub grid: [[Point3; 10]; 10],
    /// One candidate pair per corner, CCW from the parameter
    /// origin: `[(0,0), (0,9), (9,9), (9,0)]`. Pair = `(P_u,
    /// P_v)`: `P_u` drives the cross-derivative of the
    /// u-constant boundary through the corner (`u = 0` or
    /// `u = 1`), `P_v` the v-constant boundary (`v = 0` or
    /// `v = 1`) — on each boundary the blend reduces to that
    /// boundary's own candidate.
    pub inner_cps: [(Point3, Point3); 4],
}

/// The four blended grid indices.
#[inline]
fn is_dual_twist(i: usize, j: usize) -> bool {
    matches!((i, j), (1, 1) | (1, 8) | (8, 1) | (8, 8))
}

/// Map a blended index to (corner pair index, wu, wv): `wu`
/// vanishes on the u-constant boundary through the corner, `wv`
/// on the v-constant boundary, so each boundary sees only its
/// own candidate.
#[inline]
fn quadrant_weights(i: usize, j: usize, u: f64, v: f64) -> (usize, f64, f64) {
    match (i, j) {
        (1, 1) => (0, u, v),
        (1, 8) => (1, u, 1.0 - v),
        (8, 8) => (2, 1.0 - u, 1.0 - v),
        (8, 1) => (3, 1.0 - u, v),
        _ => unreachable!("quadrant_weights called with non-dual-twist index ({i}, {j})"),
    }
}

#[inline]
fn blended_cp(inner_cps: &[(Point3, Point3); 4], i: usize, j: usize, u: f64, v: f64) -> Point3 {
    let (k, wu, wv) = quadrant_weights(i, j, u, v);
    let (p_u, p_v) = inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
        // Corner limit convention (midpoint): the blended control
        // multiplies B_1(u)·B_1(v), which vanishes quadratically
        // at the corner, so this choice cannot affect any surface
        // derivative up to second order.
        return Point3::new(
            0.5 * (p_u.x + p_v.x),
            0.5 * (p_u.y + p_v.y),
            0.5 * (p_u.z + p_v.z),
        );
    }
    let c_u = wv / denom;
    let c_v = wu / denom;
    Point3::new(
        c_u * p_u.x + c_v * p_v.x,
        c_u * p_u.y + c_v * p_v.y,
        c_u * p_u.z + c_v * p_v.z,
    )
}

/// `∂(blended_cp)/∂u` via quotient rule.
#[inline]
fn blended_cp_partial_u(
    inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv) = quadrant_weights(i, j, u, v);
    let du_sign = match k {
        0 | 1 => 1.0,
        _ => -1.0,
    };
    let (p_u, p_v) = inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let factor = du_sign * wv / (denom * denom);
    Vec3::new(
        factor * (p_v.x - p_u.x),
        factor * (p_v.y - p_u.y),
        factor * (p_v.z - p_u.z),
    )
}

/// `∂(blended_cp)/∂v` via quotient rule.
#[inline]
fn blended_cp_partial_v(
    inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv) = quadrant_weights(i, j, u, v);
    let dv_sign = match k {
        0 | 3 => 1.0,
        _ => -1.0,
    };
    let (p_u, p_v) = inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let factor = dv_sign * wu / (denom * denom);
    Vec3::new(
        factor * (p_u.x - p_v.x),
        factor * (p_u.y - p_v.y),
        factor * (p_u.z - p_v.z),
    )
}

impl RectangularGregoryBiNonicDualTwist {
    /// Lift a plain bi-nonic: both candidates of every pair take
    /// the plain twist control, so evaluation is identical until
    /// the pairs are written separately.
    pub fn from_plain(plain: &RectangularGregoryBiNonic) -> Self {
        let grid = plain.grid;
        let inner_cps = [
            (grid[1][1], grid[1][1]),
            (grid[1][8], grid[1][8]),
            (grid[8][8], grid[8][8]),
            (grid[8][1], grid[8][1]),
        ];
        Self { grid, inner_cps }
    }

    /// Evaluate the patch at `(u, v) ∈ [0, 1]²`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_all(v);
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for (i, &bu_i) in bu.iter().enumerate() {
            for (j, &bv_j) in bv.iter().enumerate() {
                let w = bu_i * bv_j;
                let p = if is_dual_twist(i, j) {
                    blended_cp(&self.inner_cps, i, j, u, v)
                } else {
                    self.grid[i][j]
                };
                out = Point3::new(out.x + w * p.x, out.y + w * p.y, out.z + w * p.z);
            }
        }
        out
    }

    /// `∂S/∂u`: product rule over the tensor basis plus the
    /// rational blend's own `u`-derivative.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bu_prime = bernstein_binonic_deriv_all(u);
        let bv = bernstein_binonic_all(v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for (i, (&bu_i, &bu_prime_i)) in bu.iter().zip(bu_prime.iter()).enumerate() {
            for (j, &bv_j) in bv.iter().enumerate() {
                let (p, dp_du) = if is_dual_twist(i, j) {
                    (
                        blended_cp(&self.inner_cps, i, j, u, v),
                        blended_cp_partial_u(&self.inner_cps, i, j, u, v),
                    )
                } else {
                    (self.grid[i][j], Vec3::new(0.0, 0.0, 0.0))
                };
                let w_prime = bu_prime_i * bv_j;
                let w = bu_i * bv_j;
                out = out
                    + Vec3::new(
                        w_prime * p.x + w * dp_du.x,
                        w_prime * p.y + w * dp_du.y,
                        w_prime * p.z + w * dp_du.z,
                    );
            }
        }
        out
    }

    /// `∂S/∂v`: symmetric to [`partial_u`](Self::partial_u).
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let bu = bernstein_binonic_all(u);
        let bv = bernstein_binonic_all(v);
        let bv_prime = bernstein_binonic_deriv_all(v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for (i, &bu_i) in bu.iter().enumerate() {
            for (j, (&bv_j, &bv_prime_j)) in bv.iter().zip(bv_prime.iter()).enumerate() {
                let (p, dp_dv) = if is_dual_twist(i, j) {
                    (
                        blended_cp(&self.inner_cps, i, j, u, v),
                        blended_cp_partial_v(&self.inner_cps, i, j, u, v),
                    )
                } else {
                    (self.grid[i][j], Vec3::new(0.0, 0.0, 0.0))
                };
                let w_prime = bu_i * bv_prime_j;
                let w = bu_i * bv_j;
                out = out
                    + Vec3::new(
                        w_prime * p.x + w * dp_dv.x,
                        w_prime * p.y + w * dp_dv.y,
                        w_prime * p.z + w * dp_dv.z,
                    );
            }
        }
        out
    }

    /// `∂²S/∂u²` — central difference on the analytic first
    /// partial (the rational contribution's closed-form second
    /// derivative is deferred, matching the degree-5 precedent;
    /// the blended controls carry `B_1` weights so the FD error
    /// is far below the G2 gate's resolution).
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-4_f64;
        let u_plus = (u + h).min(1.0);
        let u_minus = (u - h).max(0.0);
        let (fwd, bwd) = (self.partial_u(u_plus, v), self.partial_u(u_minus, v));
        finite_difference(fwd, bwd, u_plus - u_minus)
    }

    /// `∂²S/∂u∂v` — central difference on `partial_u` along `v`.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-4_f64;
        let v_plus = (v + h).min(1.0);
        let v_minus = (v - h).max(0.0);
        let (fwd, bwd) = (self.partial_u(u, v_plus), self.partial_u(u, v_minus));
        finite_difference(fwd, bwd, v_plus - v_minus)
    }

    /// `∂²S/∂v²` — central difference on `partial_v` along `v`.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-4_f64;
        let v_plus = (v + h).min(1.0);
        let v_minus = (v - h).max(0.0);
        let (fwd, bwd) = (self.partial_v(u, v_plus), self.partial_v(u, v_minus));
        finite_difference(fwd, bwd, v_plus - v_minus)
    }

    /// Parameter range, `[0,1] × [0,1]`.
    pub fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, 1.0), (0.0, 1.0))
    }

    /// TEMPORARY CONTAINMENT (recorded in the dual-twist design
    /// doc and PK-085): convert via the plain bi-nonic whose four
    /// twist controls take the pair midpoints. The blended control
    /// deviates from the midpoint by at most `|P_u − P_v| / 2`,
    /// weighted by `B_1^9(u)·B_1^9(v) ≤ 0.153`, so the surface
    /// error is bounded by `0.077 · max|P_u − P_v|` — zero when
    /// the pairs coincide. The exact rational conversion is the
    /// degree-(13, 13) denominator-clearing construction used by
    /// the degree-5 variant (`rectangular_deg5_to_nurbs.rs`),
    /// staged as Stage B follow-up before any producer writes
    /// distinct pairs into committed geometry.
    pub fn to_nurbs_midpoint_containment(&self) -> crate::geometry::nurbs_surface::NurbsSurface {
        let mut grid = self.grid;
        let indices = [(1usize, 1usize, 0usize), (1, 8, 1), (8, 8, 2), (8, 1, 3)];
        for (i, j, k) in indices {
            let (p_u, p_v) = self.inner_cps[k];
            grid[i][j] = Point3::new(
                0.5 * (p_u.x + p_v.x),
                0.5 * (p_u.y + p_v.y),
                0.5 * (p_u.z + p_v.z),
            );
        }
        RectangularGregoryBiNonic { grid }.to_nurbs_exact()
    }
}

#[inline]
fn finite_difference(forward: Vec3, backward: Vec3, dh: f64) -> Vec3 {
    if dh < 1e-15 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    Vec3::new(
        (forward.x - backward.x) / dh,
        (forward.y - backward.y) / dh,
        (forward.z - backward.z) / dh,
    )
}
