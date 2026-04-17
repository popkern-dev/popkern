//! Chiyokura-Kimura rectangular Gregory patch evaluator.
//!
//! See the `gregory` module's top-level docs for background on
//! the twist-compatibility problem and the rational inner-CP
//! blending scheme.

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// A rectangular (N=4) rational Gregory patch on the parameter
/// domain `[0, 1] × [0, 1]`.
///
/// Stores a 4×4 grid of Bezier control points whose **corner**
/// entries (`boundary_cps[i][j]` for `(i, j) ∈ {(0,0), (0,3),
/// (3,0), (3,3)}`) are the four patch corners, and whose
/// **boundary** entries along each edge trace out a cubic Bezier
/// boundary curve. The four interior positions `(1,1), (1,2),
/// (2,1), (2,2)` are replaced at evaluation time by a rational
/// blend of a pair of candidate control points stored in
/// [`inner_cps`](Self::inner_cps) — one candidate per incident
/// boundary at that corner.
///
/// The pair at `inner_cps[k]` corresponds to the corner reachable
/// from patch corner `k` (ordered (0,0) → (0,3) → (3,3) → (3,0)
/// counter-clockwise). Each pair is `(P_u, P_v)` where `P_u`
/// matches the cross-tangent of the boundary in the `u` direction
/// at that corner and `P_v` matches the boundary in the `v`
/// direction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularGregory {
    /// 4×4 control-point grid stored row-major: `boundary_cps[i][j]`
    /// is the control point at Bezier index `(i, j)` with `i` the
    /// u-index and `j` the v-index. The four interior entries
    /// `(1,1)`, `(1,2)`, `(2,1)`, `(2,2)` are placeholders —
    /// evaluation replaces them with the rational blends from
    /// [`inner_cps`](Self::inner_cps).
    pub boundary_cps: [[Point3; 4]; 4],
    /// Four pairs of inner control points — one pair per corner,
    /// ordered (0,0) → (0,3) → (3,3) → (3,0) counter-clockwise.
    /// Each pair is `(P_u, P_v)` where `P_u` drives the
    /// u-direction boundary's cross-tangent at the corner and
    /// `P_v` drives the v-direction boundary's cross-tangent.
    pub inner_cps: [(Point3, Point3); 4],
}

/// Cubic Bernstein polynomials `B_i^3(t)` for `i ∈ {0, 1, 2, 3}`.
#[inline]
fn bernstein_cubic(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => s * s * s,
        1 => 3.0 * t * s * s,
        2 => 3.0 * t * t * s,
        3 => t * t * t,
        _ => 0.0,
    }
}

/// Derivative of the cubic Bernstein `dB_i^3/dt`.
#[inline]
fn bernstein_cubic_deriv(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => -3.0 * s * s,
        1 => 3.0 * s * s - 6.0 * t * s,
        2 => 6.0 * t * s - 3.0 * t * t,
        3 => 3.0 * t * t,
        _ => 0.0,
    }
}

#[inline]
fn is_inner(i: usize, j: usize) -> bool {
    matches!((i, j), (1, 1) | (1, 2) | (2, 1) | (2, 2))
}

/// Rationally-blended inner control points at parameter `(u, v)`.
///
/// The four corners are ordered (0,0) → (0,3) → (3,3) → (3,0)
/// counter-clockwise, matching [`RectangularGregory::inner_cps`].
/// The blending weight at each corner vanishes on the boundary
/// opposite its cross-tangent direction — so the corner's
/// contribution is wholly `P_u` on the v-boundary through it and
/// wholly `P_v` on the u-boundary through it.
#[inline]
fn inner_cp(inner_cps: &[(Point3, Point3); 4], i: usize, j: usize, u: f64, v: f64) -> Point3 {
    let (k, wu, wv) = match (i, j) {
        (1, 1) => (0, u, v),
        (1, 2) => (1, u, 1.0 - v),
        (2, 2) => (2, 1.0 - u, 1.0 - v),
        (2, 1) => (3, 1.0 - u, v),
        _ => unreachable!("inner_cp called with non-inner index"),
    };
    let (p_u, p_v) = inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
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

/// ∂/∂u of the rational inner-CP blend.
///
/// Let `D = wu + wv` and `P = (wv·P_u + wu·P_v) / D`. Quotient rule:
/// `∂P/∂wu = wv · (P_v - P_u) / D²`, and `∂wu/∂u = du_sign ∈ {±1}`.
fn inner_cp_partial_u(
    inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv, du_sign) = match (i, j) {
        (1, 1) => (0, u, v, 1.0),
        (1, 2) => (1, u, 1.0 - v, 1.0),
        (2, 2) => (2, 1.0 - u, 1.0 - v, -1.0),
        (2, 1) => (3, 1.0 - u, v, -1.0),
        _ => unreachable!(),
    };
    let (p_u, p_v) = inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let d2 = denom * denom;
    let factor = du_sign * wv / d2;
    Vec3::new(
        factor * (p_v.x - p_u.x),
        factor * (p_v.y - p_u.y),
        factor * (p_v.z - p_u.z),
    )
}

/// ∂/∂v of the rational inner-CP blend. `∂P/∂wv = wu · (P_u - P_v) / D²`
/// with `∂wv/∂v = dv_sign ∈ {±1}`.
fn inner_cp_partial_v(
    inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv, dv_sign) = match (i, j) {
        (1, 1) => (0, u, v, 1.0),
        (1, 2) => (1, u, 1.0 - v, -1.0),
        (2, 2) => (2, 1.0 - u, 1.0 - v, -1.0),
        (2, 1) => (3, 1.0 - u, v, 1.0),
        _ => unreachable!(),
    };
    let (p_u, p_v) = inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let d2 = denom * denom;
    let factor = dv_sign * wu / d2;
    Vec3::new(
        factor * (p_u.x - p_v.x),
        factor * (p_u.y - p_v.y),
        factor * (p_u.z - p_v.z),
    )
}

impl RectangularGregory {
    /// Evaluate the patch at `(u, v) ∈ [0, 1] × [0, 1]`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..4 {
            let bu = bernstein_cubic(i, u);
            for j in 0..4 {
                let bv = bernstein_cubic(j, v);
                let w = bu * bv;
                let p = if is_inner(i, j) {
                    inner_cp(&self.inner_cps, i, j, u, v)
                } else {
                    self.boundary_cps[i][j]
                };
                x += w * p.x;
                y += w * p.y;
                z += w * p.z;
            }
        }
        Point3::new(x, y, z)
    }

    /// Partial derivative ∂S/∂u at `(u, v)`.
    ///
    /// The inner control points depend on `(u, v)` via the rational
    /// blend, so the full derivative is
    /// `Σ_ij ∂B_i^3/∂u · B_j^3 · P_ij + B_i^3 · B_j^3 · ∂P_ij/∂u`.
    /// The second term vanishes on all four boundaries (where one
    /// Bernstein factor is zero at the inner indices) — exactly why
    /// the rational construction delivers G1 there without further
    /// adjustment. Inside the domain we compute `∂P_ij/∂u` by the
    /// quotient rule.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..4 {
            let bu_d = bernstein_cubic_deriv(i, u);
            let bu = bernstein_cubic(i, u);
            for j in 0..4 {
                let bv = bernstein_cubic(j, v);
                let (p, dp_du) = if is_inner(i, j) {
                    (
                        inner_cp(&self.inner_cps, i, j, u, v),
                        inner_cp_partial_u(&self.inner_cps, i, j, u, v),
                    )
                } else {
                    (self.boundary_cps[i][j], Vec3::new(0.0, 0.0, 0.0))
                };
                let w1 = bu_d * bv;
                let w2 = bu * bv;
                x += w1 * p.x + w2 * dp_du.x;
                y += w1 * p.y + w2 * dp_du.y;
                z += w1 * p.z + w2 * dp_du.z;
            }
        }
        Vec3::new(x, y, z)
    }

    /// Partial derivative ∂S/∂v at `(u, v)`. Symmetric to
    /// [`partial_u`](Self::partial_u).
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..4 {
            let bu = bernstein_cubic(i, u);
            for j in 0..4 {
                let bv_d = bernstein_cubic_deriv(j, v);
                let bv = bernstein_cubic(j, v);
                let (p, dp_dv) = if is_inner(i, j) {
                    (
                        inner_cp(&self.inner_cps, i, j, u, v),
                        inner_cp_partial_v(&self.inner_cps, i, j, u, v),
                    )
                } else {
                    (self.boundary_cps[i][j], Vec3::new(0.0, 0.0, 0.0))
                };
                let w1 = bu * bv_d;
                let w2 = bu * bv;
                x += w1 * p.x + w2 * dp_dv.x;
                y += w1 * p.y + w2 * dp_dv.y;
                z += w1 * p.z + w2 * dp_dv.z;
            }
        }
        Vec3::new(x, y, z)
    }
}
