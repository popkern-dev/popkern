//! §5.3 item #3 extended (the degree-7 spoke-solver work) — Karciauskas-Peters
//! **biseptic** rectangular Gregory patch (degree 7) with
//! free-curvature interior CPs. Extends
//! [`super::RectangularGregoryBiquinticFull`] from degree 5 to
//! degree 7 to provide **4 free DOF per shared spoke** for the
//! spoke-G1 joint solver (up from 2 at degree 5), closing the
//! structural `B_1^5 + B_4^5` cancellation gap documented in
//! the Phase-8 scorecard.
//!
//! # Control structure
//!
//! * `grid: [[Point3; 8]; 8]` — 64 entries. Boundary positions
//!   (28 CPs at `i ∈ {0, 7}` or `j ∈ {0, 7}`, corners counted
//!   once) are fixed by the 4 boundary Bezier curves.
//!   Corner-adjacent interior positions `(1,1), (1,6), (6,1),
//!   (6,6)` are **placeholders** — evaluation replaces them
//!   with the rational inner-CP blend from `corner_inner_cps`.
//!   Non-corner-adjacent interior positions (the remaining 32)
//!   are **free** plain-Bezier CPs, used directly.
//! * `corner_inner_cps: [(Point3, Point3); 4]` — one pair per
//!   corner CCW (`corner 0 = (0,0)`, `corner 1 = (0,7)`,
//!   `corner 2 = (7,7)`, `corner 3 = (7,0)`), same semantics
//!   as biquintic-full.
//!
//! # DOF analysis (G1 and future G2)
//!
//! At boundary `u = 0`, `B_i^7'(0)` is nonzero only for `i = 0`
//! (coef `-7`) and `i = 1` (coef `+7`). Cross-tangent:
//! `∂S/∂u(0, v) = 7 · Σ_j B_j^7(v) · (P_{1,j} − P_{0,j})`.
//! Row `i = 1` has free interior CPs at `j ∈ {2, 3, 4, 5}` — 4
//! per spoke, spanning the 4-D `{B_2^7, B_3^7, B_4^7, B_5^7}`
//! quintic subspace. The mismatch from corner-adjacent
//! rational-blend contributions lives in the `{B_1^7, B_6^7}`
//! complement — but in biseptic the boundary has ENOUGH Bernstein
//! basis coverage for row 1's free CPs to reach the full
//! `{B_2^7..B_5^7}` span, which covers the raised quintic image
//! of the biquintic `{B_1^5 + B_4^5}` gap.
//!
//! For **G2** (phase 10), `B_i^7''(0)` is nonzero for
//! `i = 0, 1, 2`. Row `i = 2` has 6 free interior CPs
//! (`j ∈ {1, 2, 3, 4, 5, 6}` — none are corner-adjacent, since
//! corner-adjacent is only at `(1,1), (1,6), (6,1), (6,6)`).
//! That's 6 more DOF per spoke for the G2 continuity constraint.
//! the polynomial-identity solver work will add `∂²S/∂u² match` linear constraints to the
//! existing solver over these additional unknowns.

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Biseptic (degree 7) rectangular Gregory patch with 32 free-
/// curvature interior CPs in addition to 4 corner inner-CP pairs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RectangularGregoryBiseptic {
    /// 8×8 control-point grid, row-major. `grid[i][j]` is at
    /// Bezier index `(i, j)`. Boundary positions are fixed.
    /// Corner-adjacent interior positions `(1,1), (1,6), (6,1),
    /// (6,6)` are placeholders — overwritten during evaluation
    /// by the rational inner-CP blend. The other 32 interior
    /// positions are free plain-Bezier CPs.
    pub grid: [[Point3; 8]; 8],
    /// Four inner CP pairs, one per corner CCW. Each pair
    /// `(P_u, P_v)` drives the cross-tangent at its corner along
    /// the u-boundary (P_u) and v-boundary (P_v).
    pub corner_inner_cps: [(Point3, Point3); 4],
}

/// Is `(i, j)` one of the 4 corner-adjacent interior positions
/// whose value is determined by the rational inner-CP blend?
#[inline]
pub(crate) fn is_corner_adjacent(i: usize, j: usize) -> bool {
    matches!((i, j), (1, 1) | (1, 6) | (6, 1) | (6, 6))
}

/// Is `(i, j)` a free plain-Bezier interior position?
///
/// Used by the spoke-G1 solver to iterate the 32 free CPs per
/// sub-patch as solve unknowns.
#[inline]
#[allow(dead_code)]
pub(crate) fn is_free_interior(i: usize, j: usize) -> bool {
    (1..=6).contains(&i) && (1..=6).contains(&j) && !is_corner_adjacent(i, j)
}

/// Septic Bernstein `B_i^7(t)` for `i ∈ {0..=7}`.
#[inline]
pub(crate) fn bernstein_septic(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => s.powi(7),
        1 => 7.0 * t * s.powi(6),
        2 => 21.0 * t * t * s.powi(5),
        3 => 35.0 * t.powi(3) * s.powi(4),
        4 => 35.0 * t.powi(4) * s.powi(3),
        5 => 21.0 * t.powi(5) * s * s,
        6 => 7.0 * t.powi(6) * s,
        7 => t.powi(7),
        _ => 0.0,
    }
}

#[inline]
fn bernstein_sextic(i: i32, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => s.powi(6),
        1 => 6.0 * t * s.powi(5),
        2 => 15.0 * t * t * s.powi(4),
        3 => 20.0 * t.powi(3) * s.powi(3),
        4 => 15.0 * t.powi(4) * s * s,
        5 => 6.0 * t.powi(5) * s,
        6 => t.powi(6),
        _ => 0.0,
    }
}

/// `B_i^7'(t) = 7 · (B_{i-1}^6(t) − B_i^6(t))` with `B_{-1}^6 = B_7^6 = 0`.
#[inline]
pub(crate) fn bernstein_septic_deriv(i: usize, t: f64) -> f64 {
    7.0 * (bernstein_sextic(i as i32 - 1, t) - bernstein_sextic(i as i32, t))
}

/// For a corner-adjacent interior index `(i, j)` at degree 7,
/// return the corresponding corner index `k` and the rational-
/// blend weights `(wu, wv)` used in `corner_adjacent_blend`.
#[inline]
fn corner_adjacent_weights(i: usize, j: usize, u: f64, v: f64) -> (usize, f64, f64) {
    match (i, j) {
        (1, 1) => (0, u, v),
        (1, 6) => (1, u, 1.0 - v),
        (6, 6) => (2, 1.0 - u, 1.0 - v),
        (6, 1) => (3, 1.0 - u, v),
        _ => {
            unreachable!("corner_adjacent_weights called with non-corner-adjacent index ({i}, {j})")
        }
    }
}

#[inline]
fn corner_adjacent_blend(
    corner_inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Point3 {
    let (k, wu, wv) = corner_adjacent_weights(i, j, u, v);
    let (p_u, p_v) = corner_inner_cps[k];
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

#[inline]
fn corner_adjacent_blend_partial_u(
    corner_inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv) = corner_adjacent_weights(i, j, u, v);
    let du_sign = match k {
        0 | 1 => 1.0,
        2 | 3 => -1.0,
        _ => unreachable!(),
    };
    let (p_u, p_v) = corner_inner_cps[k];
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
fn corner_adjacent_blend_partial_v(
    corner_inner_cps: &[(Point3, Point3); 4],
    i: usize,
    j: usize,
    u: f64,
    v: f64,
) -> Vec3 {
    let (k, wu, wv) = corner_adjacent_weights(i, j, u, v);
    let dv_sign = match k {
        0 | 3 => 1.0,
        1 | 2 => -1.0,
        _ => unreachable!(),
    };
    let (p_u, p_v) = corner_inner_cps[k];
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

impl RectangularGregoryBiseptic {
    /// Evaluate the patch at `(u, v) ∈ [0, 1]²`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..8 {
            let bu = bernstein_septic(i, u);
            for j in 0..8 {
                let bv = bernstein_septic(j, v);
                let w = bu * bv;
                let p = if is_corner_adjacent(i, j) {
                    corner_adjacent_blend(&self.corner_inner_cps, i, j, u, v)
                } else {
                    self.grid[i][j]
                };
                x += w * p.x;
                y += w * p.y;
                z += w * p.z;
            }
        }
        Point3::new(x, y, z)
    }

    /// `∂S/∂u` at `(u, v)`.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..8 {
            let bu = bernstein_septic(i, u);
            let bu_prime = bernstein_septic_deriv(i, u);
            for j in 0..8 {
                let bv = bernstein_septic(j, v);
                let (p, dp_du) = if is_corner_adjacent(i, j) {
                    (
                        corner_adjacent_blend(&self.corner_inner_cps, i, j, u, v),
                        corner_adjacent_blend_partial_u(&self.corner_inner_cps, i, j, u, v),
                    )
                } else {
                    (self.grid[i][j], Vec3::new(0.0, 0.0, 0.0))
                };
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

    /// `∂S/∂v` at `(u, v)`.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..8 {
            let bu = bernstein_septic(i, u);
            for j in 0..8 {
                let bv = bernstein_septic(j, v);
                let bv_prime = bernstein_septic_deriv(j, v);
                let (p, dp_dv) = if is_corner_adjacent(i, j) {
                    (
                        corner_adjacent_blend(&self.corner_inner_cps, i, j, u, v),
                        corner_adjacent_blend_partial_v(&self.corner_inner_cps, i, j, u, v),
                    )
                } else {
                    (self.grid[i][j], Vec3::new(0.0, 0.0, 0.0))
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

    /// `∂²S/∂u²` — central-difference on `partial_u`. Used by
    /// the G2 solver (the polynomial-identity solver work) and curvature measures.
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

    /// `∂²S/∂u∂v` — central-difference on `partial_u` along `v`.
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

    /// `∂²S/∂v²` — central-difference on `partial_v`.
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
}
