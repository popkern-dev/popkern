//! §5.3 item #3 — Karciauskas-Peters biquintic rectangular
//! Gregory patch with **free-curvature interior CPs**.
//!
//! Extends [`super::RectangularGregoryDeg5`] (which uses 4
//! corner inner-CP pairs over a 6×6 boundary grid, blending ALL
//! 16 interior positions rationally from the 4 corners) by
//! replacing 12 of the 16 interior positions with FREE plain-
//! Bezier control points. Only the 4 **corner-adjacent**
//! interior positions — `(1,1)`, `(1,4)`, `(4,1)`, `(4,4)` —
//! retain the rational blend that gives Gregory its twist-
//! compatibility property at valence-N vertices.
//!
//! The 12 free CPs provide the additional degrees of freedom
//! needed to enforce **geometric G1 along the entire shared
//! spoke** between adjacent sub-patches — not just at the
//! central vertex C (which the simplified Chiyokura-Kimura +
//! Shirman-Séquin cubic Gregory achieves). With 4-inner-CP-
//! pair Gregory (cubic or biquintic), every interior Bezier
//! position is rationally blended from the 4 corner pairs; the
//! per-corner Chiyokura-Kimura twist consumes all the DOF and
//! the spoke-interior cross-tangent is fixed by interpolation,
//! often with a significant normal-mismatch vs. the neighbour.
//! The K-P free-curvature form decouples these, giving a joint
//! linear system over the 12 per-patch free CPs that can
//! enforce spoke G1 at multiple parameter samples along each
//! shared boundary.
//!
//! # Parameter domain
//!
//! `[0, 1] × [0, 1]`, same as the other rectangular variants.
//!
//! # Control structure
//!
//! * `grid: [[Point3; 6]; 6]` — 36 entries. Boundary positions
//!   (20 CPs at `i ∈ {0, 5}` or `j ∈ {0, 5}`) are fixed by the
//!   4 boundary Bezier curves. Corner-adjacent interior
//!   positions `(1,1), (1,4), (4,1), (4,4)` are **placeholders**
//!   — evaluation replaces them with the rational inner-CP
//!   blend from `corner_inner_cps`. Non-corner-adjacent
//!   interior positions (the remaining 12) are **free** plain-
//!   Bezier CPs, used directly.
//! * `corner_inner_cps: [(Point3, Point3); 4]` — one pair per
//!   corner CCW (`corner 0 = (0,0)`, `corner 1 = (0,5)`,
//!   `corner 2 = (5,5)`, `corner 3 = (5,0)`), same meaning as
//!   [`super::RectangularGregoryDeg5`].
//!
//! The 4 corner-adjacent blend positions use the same quadrant-
//! mapped rational-blend formula that corners drive in the 4-
//! CP-pair form; the 12 free positions break the rational-
//! blend inheritance and are user-supplied.
//!
//! # References
//!
//! * Karciauskas & Peters (2003), "Biquintic C1 surface splines
//!   over arbitrary quadrilateral meshes."
//! * Peters, J. (2002), "Geometric continuity." In *Handbook of
//!   Computer Aided Geometric Design*, Elsevier.
//! * Chiyokura-Kimura (1983) — original rational Gregory
//!   (4-inner-CP-pair form).

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Biquintic rectangular Gregory patch with 12 free-curvature
/// interior CPs in addition to 4 corner inner-CP pairs.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RectangularGregoryBiquinticFull {
    /// 6×6 control-point grid, row-major. `grid[i][j]` is at
    /// Bezier index `(i, j)`. Boundary positions are fixed.
    /// Corner-adjacent interior positions `(1,1), (1,4), (4,1),
    /// (4,4)` are placeholders — overwritten during evaluation
    /// by the rational inner-CP blend. The other 12 interior
    /// positions are free plain-Bezier CPs.
    pub grid: [[Point3; 6]; 6],
    /// Four inner CP pairs, one per corner CCW. Each pair
    /// `(P_u, P_v)` drives the cross-tangent at its corner along
    /// the u-boundary (P_u) and v-boundary (P_v).
    pub corner_inner_cps: [(Point3, Point3); 4],
}

/// Is `(i, j)` one of the 4 corner-adjacent interior positions
/// whose value is determined by the rational inner-CP blend?
#[inline]
pub(crate) fn is_corner_adjacent(i: usize, j: usize) -> bool {
    matches!((i, j), (1, 1) | (1, 4) | (4, 1) | (4, 4))
}

/// Is `(i, j)` a free plain-Bezier interior position?
///
/// Used by the spoke-G1 solver (Phase 6) to iterate the 12 free
/// CPs per sub-patch as solve unknowns.
#[inline]
#[allow(dead_code)]
pub(crate) fn is_free_interior(i: usize, j: usize) -> bool {
    (1..=4).contains(&i) && (1..=4).contains(&j) && !is_corner_adjacent(i, j)
}

/// Quintic Bernstein `B_i^5(t)` for `i ∈ {0..5}`.
#[inline]
fn bernstein_quintic(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
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

#[inline]
fn bernstein_quintic_deriv(i: usize, t: f64) -> f64 {
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

/// For a corner-adjacent interior index `(i, j)`, return the
/// corresponding corner index `k` and the rational-blend
/// weights `(wu, wv)` used in `corner_adjacent_blend`.
///
/// The weights vanish on one of the two incident boundary edges
/// (the one passing through the corner in the u or v direction)
/// so that the inner-CP pair's `P_u` / `P_v` selection picks the
/// correct cross-tangent direction on each boundary.
#[inline]
fn corner_adjacent_weights(i: usize, j: usize, u: f64, v: f64) -> (usize, f64, f64) {
    match (i, j) {
        // Corner 0 at (0,0): weights (u, v) — both vanish at the corner.
        (1, 1) => (0, u, v),
        // Corner 1 at (0,5): weights (u, 1-v).
        (1, 4) => (1, u, 1.0 - v),
        // Corner 2 at (5,5): weights (1-u, 1-v).
        (4, 4) => (2, 1.0 - u, 1.0 - v),
        // Corner 3 at (5,0): weights (1-u, v).
        (4, 1) => (3, 1.0 - u, v),
        _ => {
            unreachable!("corner_adjacent_weights called with non-corner-adjacent index ({i}, {j})")
        }
    }
}

/// Rational blend at a corner-adjacent interior position:
/// `P = (wv · P_u + wu · P_v) / (wu + wv)` with a singularity
/// tie-break at `wu + wv → 0` that resolves to `(P_u + P_v) / 2`.
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

impl RectangularGregoryBiquinticFull {
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
        for i in 0..6 {
            let bu = bernstein_quintic(i, u);
            let bu_prime = bernstein_quintic_deriv(i, u);
            for j in 0..6 {
                let bv = bernstein_quintic(j, v);
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
        for i in 0..6 {
            let bu = bernstein_quintic(i, u);
            for j in 0..6 {
                let bv = bernstein_quintic(j, v);
                let bv_prime = bernstein_quintic_deriv(j, v);
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

    /// `∂²S/∂u²` — central-difference on `partial_u`.
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
