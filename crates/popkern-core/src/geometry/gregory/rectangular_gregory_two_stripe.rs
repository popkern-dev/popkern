//! §5.3 item #3 Phase-10 — Peters two-stripe bidegree-7 patch.
//!
//! The structural fix for the `{B_1^n, B_{n-1}^n}` Bernstein
//! modes that the single-stripe biseptic can't reach:
//! **remove the rational blend** at the 4 corner-adjacent
//! positions and make them plain-Bezier free CPs. This gives
//! the spoke-G1 solver access to the full row-1 Bernstein
//! basis (`{B_1^7..B_6^7}` — 6 DOF per spoke side, up from 4
//! in biseptic / 2 in biquintic-full).
//!
//! The name "two-stripe" comes from the Peters (2000)
//! construction: instead of using ONE stripe (row 1) of
//! interior CPs where `(1,1), (1,n-1)` are rational-blend
//! placeholders and the rest are free, we use the FIRST TWO
//! STRIPES (rows 1 & 2) with NO rational-blend placeholders.
//! The corner twist compatibility that was previously
//! enforced by the rational blend is now enforced by the
//! joint solver via explicit corner-G1 constraints between
//! the `(1,1)` plain-Bezier CP and the adjacent row-0 / row-2
//! positions.
//!
//! # Control structure
//!
//! * `grid: [[Point3; 8]; 8]` — 64 plain-Bezier CPs. No
//!   placeholders. No rational-blend overrides. The patch is
//!   a pure bidegree-7 Bezier tensor-product surface, which
//!   greatly simplifies evaluation, partial derivatives, and
//!   the exact NURBS conversion (no denominator clearing —
//!   the patch IS a rational Bezier with unit weights, i.e.
//!   an ordinary polynomial Bezier).
//!
//! # Spoke-G1 DOF per side
//!
//! At boundary `u = 0`, `∂S/∂u(0, v)` depends on rows `i = 0`
//! (boundary, fixed) and `i = 1` (all plain-Bezier free in
//! the two-stripe form). Free row-1 indices: `j ∈ {1..=6}`
//! (j=0, j=7 are corners, shared with neighbours). 6 DOF per
//! spoke side — the solver's reachable span is
//! `{B_1^7..B_6^7}` which covers every possible boundary-
//! cross-tangent mismatch pattern. Full-rank G1.
//!
//! # G2 DOF (the polynomial-identity solver workb)
//!
//! `∂²S/∂u²(0, v)` depends on rows `i = 0, 1, 2`. Row 2 has
//! free interior CPs at `j ∈ {1..=6}` (6 per spoke) plus the
//! non-shared corners at `(2, 0)` and `(2, 7)` which ARE
//! boundary. Combined with row 1's 6 free CPs, spoke G2
//! has 12 DOF per side over 2 constraint dimensions (n · xtan
//! and n · xtan2) — comfortably full-rank.

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Peters two-stripe bidegree-7 Bezier patch. 64 plain-Bezier
/// control points; no rational-blend corner-adjacent
/// placeholders. Corner tangent-plane compatibility is
/// enforced by the joint solver's corner-G1 constraints rather
/// than by an evaluator-level rational blend.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RectangularGregoryTwoStripe {
    /// 8×8 control-point grid, row-major. All 64 entries are
    /// plain-Bezier CPs with no special semantics — unlike
    /// biseptic's rational-blend placeholders at `(1,1),
    /// (1,6), (6,1), (6,6)`.
    pub grid: [[Point3; 8]; 8],
}

/// Septic Bernstein `B_i^7(t)` for `i ∈ {0..=7}`. Shared
/// implementation with [`super::rectangular_biseptic`].
#[inline]
pub(crate) fn bernstein_septic(i: usize, t: f64) -> f64 {
    super::rectangular_biseptic::bernstein_septic(i, t)
}

/// `B_i^7'(t)` shared with [`super::rectangular_biseptic`].
#[inline]
pub(crate) fn bernstein_septic_deriv(i: usize, t: f64) -> f64 {
    super::rectangular_biseptic::bernstein_septic_deriv(i, t)
}

impl RectangularGregoryTwoStripe {
    /// Evaluate the patch at `(u, v) ∈ [0, 1]²`. Pure
    /// bidegree-7 tensor-product Bezier — no rational blend.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..8 {
            let bu = bernstein_septic(i, u);
            for j in 0..8 {
                let bv = bernstein_septic(j, v);
                let w = bu * bv;
                let p = self.grid[i][j];
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
            let bu_prime = bernstein_septic_deriv(i, u);
            for j in 0..8 {
                let bv = bernstein_septic(j, v);
                let w = bu_prime * bv;
                let p = self.grid[i][j];
                out = out + Vec3::new(w * p.x, w * p.y, w * p.z);
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
                let bv_prime = bernstein_septic_deriv(j, v);
                let w = bu * bv_prime;
                let p = self.grid[i][j];
                out = out + Vec3::new(w * p.x, w * p.y, w * p.z);
            }
        }
        out
    }

    /// `∂²S/∂u²` at `(u, v)` — analytic (no finite differencing
    /// needed now that evaluation is plain Bezier).
    pub fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..8 {
            let bu_pp = bernstein_septic_second_deriv(i, u);
            for j in 0..8 {
                let bv = bernstein_septic(j, v);
                let w = bu_pp * bv;
                let p = self.grid[i][j];
                out = out + Vec3::new(w * p.x, w * p.y, w * p.z);
            }
        }
        out
    }

    /// `∂²S/∂u∂v` at `(u, v)`.
    pub fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..8 {
            let bu_p = bernstein_septic_deriv(i, u);
            for j in 0..8 {
                let bv_p = bernstein_septic_deriv(j, v);
                let w = bu_p * bv_p;
                let p = self.grid[i][j];
                out = out + Vec3::new(w * p.x, w * p.y, w * p.z);
            }
        }
        out
    }

    /// `∂²S/∂v²` at `(u, v)`.
    pub fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..8 {
            let bu = bernstein_septic(i, u);
            for j in 0..8 {
                let bv_pp = bernstein_septic_second_deriv(j, v);
                let w = bu * bv_pp;
                let p = self.grid[i][j];
                out = out + Vec3::new(w * p.x, w * p.y, w * p.z);
            }
        }
        out
    }
}

/// `B_i^7''(t) = 7 · 6 · (B_{i-2}^5(t) − 2 B_{i-1}^5(t) + B_i^5(t))`
/// with out-of-range terms = 0.
#[inline]
pub(crate) fn bernstein_septic_second_deriv(i: usize, t: f64) -> f64 {
    let b5 = |i: i32| -> f64 {
        let s = 1.0 - t;
        match i {
            0 => s.powi(5),
            1 => 5.0 * t * s.powi(4),
            2 => 10.0 * t * t * s.powi(3),
            3 => 10.0 * t.powi(3) * s * s,
            4 => 5.0 * t.powi(4) * s,
            5 => t.powi(5),
            _ => 0.0,
        }
    };
    let ii = i as i32;
    42.0 * (b5(ii - 2) - 2.0 * b5(ii - 1) + b5(ii))
}
