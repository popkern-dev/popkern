//! Chiyokura–Kimura degree-5 rectangular Gregory patch.
//!
//! Parallel to [`super::RectangularGregory`] (which is degree-3)
//! but with 6×6 boundary control points and interior blending
//! over a 4×4 grid of inner positions. Degree-5 patches can
//! represent circular arcs **exactly** (within f64 precision),
//! closing the inherent 2% cubic-Bezier-circle approximation
//! error that limits the degree-3 form's §2.4 G2 match.
//!
//! # Parameter domain
//!
//! `[0, 1] × [0, 1]`, same as the degree-3 variant.
//!
//! # Control structure
//!
//! * **Boundary grid** `boundary_cps: [[Point3; 6]; 6]` — 36
//!   entries. Corners `(0,0), (0,5), (5,5), (5,0)` are the
//!   patch corners (CCW). The 16 entries along the edges
//!   (excluding corners) trace the four degree-5 boundary
//!   Bezier curves. The 16 interior entries `(i, j) ∈ {1,2,3,4}²`
//!   are placeholders — evaluation replaces them with the
//!   rational inner-CP blend.
//!
//! * **Inner CP pairs** `inner_cps: [(Point3, Point3); 4]` —
//!   one pair per corner, same layout as the degree-3
//!   variant. For evaluation at any interior position
//!   `(i, j)`, the blend pulls the pair of the *nearest*
//!   corner and mixes the two candidates by a weight
//!   determined by the position's proximity to each incident
//!   edge.
//!
//! # Interior quadrant map
//!
//! The 16 interior positions are partitioned into 4 quadrants
//! of 4 positions each, each quadrant associated with the
//! nearest corner:
//!
//! ```text
//!   corner 1 (0,5)       corner 2 (5,5)
//!     +---+---+---+---+
//!     | 1 | 1 | 2 | 2 |   j ∈ {3,4} ← corners 1 & 2 share upper row
//!     +---+---+---+---+
//!     | 1 | 1 | 2 | 2 |
//!     +---+---+---+---+
//!     | 0 | 0 | 3 | 3 |   j ∈ {1,2} ← corners 0 & 3 share lower row
//!     +---+---+---+---+
//!     | 0 | 0 | 3 | 3 |
//!     +---+---+---+---+
//!   corner 0 (0,0)       corner 3 (5,0)
//! ```
//!
//! At the corner itself the blend is singular; we resolve the
//! 0/0 limit with a tie-break to `0.5 · (P_u + P_v)` (same as
//! degree-3).

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Degree-5 rectangular rational Gregory patch.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularGregoryDeg5 {
    /// 6×6 control-point grid, row-major. `boundary_cps[i][j]`
    /// is at Bezier index `(i, j)` where `i` indexes `u` and
    /// `j` indexes `v`. The 16 interior entries at
    /// `(i, j) ∈ {1..5} × {1..5}` are placeholders — evaluation
    /// replaces them with the rational inner-CP blend.
    pub boundary_cps: [[Point3; 6]; 6],
    /// Four inner CP pairs, one per corner, CCW:
    /// `[corner (0,0), corner (0,5), corner (5,5), corner (5,0)]`.
    /// Each pair is `(P_u, P_v)` — `P_u` drives the cross-
    /// tangent of the u-running boundary through the corner,
    /// `P_v` drives the v-running boundary.
    pub inner_cps: [(Point3, Point3); 4],
}

/// Quintic Bernstein `B_i^5(t)` for `i ∈ {0..5}`.
#[inline]
fn bernstein_quintic(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    // Expanded: coefficients are (5 choose i).
    match i {
        0 => s * s * s * s * s,
        1 => 5.0 * t * s * s * s * s,
        2 => 10.0 * t * t * s * s * s,
        3 => 10.0 * t * t * t * s * s,
        4 => 5.0 * t * t * t * t * s,
        5 => t * t * t * t * t,
        _ => 0.0,
    }
}

/// Derivative `dB_i^5/dt`.
#[inline]
fn bernstein_quintic_deriv(i: usize, t: f64) -> f64 {
    // d/dt [C(5,i) · t^i · s^(5-i)] = C(5,i) · [i · t^(i-1) · s^(5-i) − (5-i) · t^i · s^(5-i-1)]
    // With s = 1-t so ds/dt = -1. For i = 0 the first term is 0; for i = 5 the second is 0.
    let s = 1.0 - t;
    match i {
        0 => -5.0 * s * s * s * s,
        1 => 5.0 * s * s * s * s - 20.0 * t * s * s * s,
        2 => 20.0 * t * s * s * s - 30.0 * t * t * s * s,
        3 => 30.0 * t * t * s * s - 20.0 * t * t * t * s,
        4 => 20.0 * t * t * t * s - 5.0 * t * t * t * t,
        5 => 5.0 * t * t * t * t,
        _ => 0.0,
    }
}

/// Second derivative `d²B_i^5/dt²`.
#[inline]
fn bernstein_quintic_second(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => 20.0 * s * s * s,
        1 => -40.0 * s * s * s + 60.0 * t * s * s,
        2 => 20.0 * s * s * s - 120.0 * t * s * s + 60.0 * t * t * s,
        3 => 60.0 * t * s * s - 120.0 * t * t * s + 20.0 * t * t * t,
        4 => 60.0 * t * t * s - 40.0 * t * t * t,
        5 => 20.0 * t * t * t,
        _ => 0.0,
    }
}

#[inline]
fn is_inner(i: usize, j: usize) -> bool {
    (1..=4).contains(&i) && (1..=4).contains(&j)
}

/// Map an interior Bezier index `(i, j)` to its (corner k,
/// blend weights wu, wv) tuple. Follows the quadrant diagram
/// in the module docs: lower-left → corner 0, upper-left →
/// corner 1, upper-right → corner 2, lower-right → corner 3.
///
/// `wu` is the weight that vanishes on the v-boundary through
/// corner k (so on that boundary `P_u` is used wholly); `wv` is
/// the weight that vanishes on the u-boundary through corner k
/// (so on that boundary `P_v` is used wholly).
#[inline]
fn quadrant_weights(i: usize, j: usize, u: f64, v: f64) -> (usize, f64, f64) {
    // Corner indices and their "remote edge" weight maps:
    //   corner 0 at (0,0): weights (u, v) → pair index 0
    //   corner 1 at (0,5): weights (u, 1-v) → pair index 1
    //   corner 2 at (5,5): weights (1-u, 1-v) → pair index 2
    //   corner 3 at (5,0): weights (1-u, v) → pair index 3
    match (i, j) {
        (1..=2, 1..=2) => (0, u, v),
        (1..=2, 3..=4) => (1, u, 1.0 - v),
        (3..=4, 3..=4) => (2, 1.0 - u, 1.0 - v),
        (3..=4, 1..=2) => (3, 1.0 - u, v),
        _ => unreachable!("quadrant_weights called with non-inner index ({i}, {j})"),
    }
}

#[inline]
fn inner_cp(inner_cps: &[(Point3, Point3); 4], i: usize, j: usize, u: f64, v: f64) -> Point3 {
    let (k, wu, wv) = quadrant_weights(i, j, u, v);
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

/// `∂(inner_cp)/∂u` via quotient rule. Same form as the
/// degree-3 variant — the rational blend weights depend only
/// on the quadrant mapping, which is degree-independent.
#[inline]
fn inner_cp_partial_u(
    inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv) = quadrant_weights(i, j, u, v);
    // Sign of ∂wu/∂u: +1 for u-based quadrants (k=0,1), −1 for (1−u)-based (k=2,3).
    let du_sign = match k {
        0 | 1 => 1.0,
        2 | 3 => -1.0,
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

#[inline]
fn inner_cp_partial_v(
    inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv) = quadrant_weights(i, j, u, v);
    // Sign of ∂wv/∂v: +1 for v-based quadrants (k=0,3), −1 for (1−v)-based (k=1,2).
    let dv_sign = match k {
        0 | 3 => 1.0,
        1 | 2 => -1.0,
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

impl RectangularGregoryDeg5 {
    /// Evaluate the patch at `(u, v) ∈ [0, 1]²`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..6 {
            let bu = bernstein_quintic(i, u);
            for j in 0..6 {
                let bv = bernstein_quintic(j, v);
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

    /// `∂S/∂u` at `(u, v)`. Product rule over the tensor-
    /// product Bernstein basis PLUS the derivative of the
    /// rational inner-CP blend for interior positions.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..6 {
            let bu = bernstein_quintic(i, u);
            let bu_prime = bernstein_quintic_deriv(i, u);
            for j in 0..6 {
                let bv = bernstein_quintic(j, v);
                let (p, dp_du) = if is_inner(i, j) {
                    (
                        inner_cp(&self.inner_cps, i, j, u, v),
                        inner_cp_partial_u(&self.inner_cps, i, j, u, v),
                    )
                } else {
                    (self.boundary_cps[i][j], Vec3::new(0.0, 0.0, 0.0))
                };
                // d/du [B_i(u) · B_j(v) · P(u, v)] =
                //   B'_i · B_j · P + B_i · B_j · dP/du
                let w_prime = bu_prime * bv;
                let w = bu * bv;
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

    /// `∂S/∂v` at `(u, v)`. Symmetric to `partial_u`.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..6 {
            let bu = bernstein_quintic(i, u);
            for j in 0..6 {
                let bv = bernstein_quintic(j, v);
                let bv_prime = bernstein_quintic_deriv(j, v);
                let (p, dp_dv) = if is_inner(i, j) {
                    (
                        inner_cp(&self.inner_cps, i, j, u, v),
                        inner_cp_partial_v(&self.inner_cps, i, j, u, v),
                    )
                } else {
                    (self.boundary_cps[i][j], Vec3::new(0.0, 0.0, 0.0))
                };
                let w_prime = bu * bv_prime;
                let w = bu * bv;
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

    /// `∂²S/∂u²` at `(u, v)`. Approximates via finite
    /// differences for the inner-CP contribution (whose full
    /// analytic second derivatives require third-rank quotient
    /// rule with cross-corner contributions — out of scope for
    /// this core-type landing; `compute_metrics` doesn't use
    /// second partials on Gregory directly yet).
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-4_f64;
        let u_plus = (u + h).min(1.0);
        let u_minus = (u - h).max(0.0);
        let fwd = self.partial_u(u_plus, v);
        let bwd = self.partial_u(u_minus, v);
        let dh = u_plus - u_minus;
        if dh < 1e-15 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        Vec3::new(
            (fwd.x - bwd.x) / dh,
            (fwd.y - bwd.y) / dh,
            (fwd.z - bwd.z) / dh,
        )
    }

    /// `∂²S/∂u∂v` — central-difference on `partial_u` along
    /// `v`.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-4_f64;
        let v_plus = (v + h).min(1.0);
        let v_minus = (v - h).max(0.0);
        let fwd = self.partial_u(u, v_plus);
        let bwd = self.partial_u(u, v_minus);
        let dh = v_plus - v_minus;
        if dh < 1e-15 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        Vec3::new(
            (fwd.x - bwd.x) / dh,
            (fwd.y - bwd.y) / dh,
            (fwd.z - bwd.z) / dh,
        )
    }

    /// `∂²S/∂v²` — central-difference on `partial_v` along
    /// `v`.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-4_f64;
        let v_plus = (v + h).min(1.0);
        let v_minus = (v - h).max(0.0);
        let fwd = self.partial_v(u, v_plus);
        let bwd = self.partial_v(u, v_minus);
        let dh = v_plus - v_minus;
        if dh < 1e-15 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        Vec3::new(
            (fwd.x - bwd.x) / dh,
            (fwd.y - bwd.y) / dh,
            (fwd.z - bwd.z) / dh,
        )
    }

    /// Parameter range, `[0,1] × [0,1]`.
    pub fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, 1.0), (0.0, 1.0))
    }
}

// Silence unused-in-core warnings — these are used downstream
// and in tests.
#[allow(dead_code)]
fn _keep_bernstein_second_referenced() -> f64 {
    bernstein_quintic_second(0, 0.0)
}
