//! Longhi triangular Gregory patch evaluator.
//!
//! See the `gregory` module's top-level docs for background on
//! the twist-compatibility problem and the rational inner-CP
//! blending scheme.

use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// A triangular (N=3) rational Gregory patch (Longhi form) on the
/// barycentric simplex `{(u, v, w) : u + v + w = 1, u, v, w ≥ 0}`.
///
/// Stores 9 **boundary** control points (3 corners + 6 edge-thirds
/// CPs) and 3 **inner** control-point candidates — one candidate
/// associated with each corner. The single interior Bezier control
/// point `P_111` of an ordinary cubic triangular Bezier is replaced
/// by a barycentric-rational blend of these 3 candidates:
/// `P_111(u, v, w) = (v·w·P_A + u·w·P_B + u·v·P_C) /
///                    (v·w + u·w + u·v)`, where `P_A`, `P_B`, `P_C`
/// are the candidates for corners 0, 1, 2 respectively. The weight
/// `v·w` vanishes on the two edges through corner 0 — so the
/// patch's cross-derivative across edge 01 (v varies, w = 0) and
/// edge 20 (w varies, v = 0) is determined by `P_A` alone at corner
/// 0, independent of the other candidates. This delivers G1 across
/// each edge with arbitrary neighbours (Longhi 1987).
///
/// The 9 boundary CPs are stored in triangular-Bezier index order
/// by `(i, j, k)` with `i + j + k = 3`:
/// `[300, 030, 003, 210, 120, 021, 012, 102, 201]` —
/// corner0, corner1, corner2, then two CPs along each edge
/// (edge 01: 210 and 120; edge 12: 021 and 012; edge 20: 102 and 201).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriangularGregory {
    /// The 9 boundary control points in triangular-Bezier index
    /// order: three corners followed by six edge-thirds CPs.
    pub boundary_cps: [Point3; 9],
    /// Three interior control-point candidates — one per corner,
    /// in corner order `(0, 1, 2)`. Blended barycentrically at
    /// evaluation time per the Longhi rational scheme.
    pub inner_cps: [Point3; 3],
}

/// Cubic triangular-Bezier Bernstein value
/// `B_{ijk}^3(u, v, w) = 3! / (i! j! k!) · u^i · v^j · w^k`.
#[inline]
fn bernstein(i: usize, j: usize, k: usize, u: f64, v: f64, w: f64) -> f64 {
    let multinomial = match (i, j, k) {
        (3, 0, 0) | (0, 3, 0) | (0, 0, 3) => 1.0,
        (2, 1, 0) | (2, 0, 1) | (1, 2, 0) | (0, 2, 1) | (0, 1, 2) | (1, 0, 2) => 3.0,
        (1, 1, 1) => 6.0,
        _ => 0.0,
    };
    multinomial * u.powi(i as i32) * v.powi(j as i32) * w.powi(k as i32)
}

/// Partial of `B_{ijk}^3(u, v, w)` with respect to one of the
/// barycentric coordinates (`which` is 0 for u, 1 for v, 2 for w).
/// Each barycentric is treated as independent here; the caller
/// handles the `w = 1 − u − v` reduction via the
/// [`TriangularGregory::partial_u_uv`] / `_v_uv` wrappers.
fn bernstein_partial(i: usize, j: usize, k: usize, u: f64, v: f64, w: f64, which: usize) -> f64 {
    let multinomial = match (i, j, k) {
        (3, 0, 0) | (0, 3, 0) | (0, 0, 3) => 1.0,
        (2, 1, 0) | (2, 0, 1) | (1, 2, 0) | (0, 2, 1) | (0, 1, 2) | (1, 0, 2) => 3.0,
        (1, 1, 1) => 6.0,
        _ => 0.0,
    };
    match which {
        0 => {
            if i == 0 {
                return 0.0;
            }
            multinomial * (i as f64) * u.powi((i - 1) as i32) * v.powi(j as i32) * w.powi(k as i32)
        }
        1 => {
            if j == 0 {
                return 0.0;
            }
            multinomial * (j as f64) * u.powi(i as i32) * v.powi((j - 1) as i32) * w.powi(k as i32)
        }
        2 => {
            if k == 0 {
                return 0.0;
            }
            multinomial * (k as f64) * u.powi(i as i32) * v.powi(j as i32) * w.powi((k - 1) as i32)
        }
        _ => 0.0,
    }
}

/// Barycentric-rational blend of the 3 inner candidates.
#[inline]
fn inner_cp(inner: &[Point3; 3], u: f64, v: f64, w: f64) -> Point3 {
    let wa = v * w;
    let wb = u * w;
    let wc = u * v;
    let denom = wa + wb + wc;
    if denom < 1e-15 {
        return Point3::new(
            (inner[0].x + inner[1].x + inner[2].x) / 3.0,
            (inner[0].y + inner[1].y + inner[2].y) / 3.0,
            (inner[0].z + inner[1].z + inner[2].z) / 3.0,
        );
    }
    let ca = wa / denom;
    let cb = wb / denom;
    let cc = wc / denom;
    Point3::new(
        ca * inner[0].x + cb * inner[1].x + cc * inner[2].x,
        ca * inner[0].y + cb * inner[1].y + cc * inner[2].y,
        ca * inner[0].z + cb * inner[1].z + cc * inner[2].z,
    )
}

/// Barycentric partial of the Longhi inner-CP blend.
///
/// Quotient rule: if `P = (wa·A + wb·B + wc·C) / D` with
/// `D = wa + wb + wc`, then
/// `∂P/∂ξ = [(∂wa/∂ξ)·A + (∂wb/∂ξ)·B + (∂wc/∂ξ)·C] / D − P · (∂D/∂ξ) / D`.
fn inner_cp_partial(inner: &[Point3; 3], u: f64, v: f64, w: f64, which: usize) -> Vec3 {
    let wa = v * w;
    let wb = u * w;
    let wc = u * v;
    let denom = wa + wb + wc;
    if denom < 1e-15 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let (dwa, dwb, dwc) = match which {
        0 => (0.0, w, v), // d/du
        1 => (w, 0.0, u), // d/dv
        2 => (v, u, 0.0), // d/dw
        _ => (0.0, 0.0, 0.0),
    };
    let d_denom = dwa + dwb + dwc;
    let blended = inner_cp(inner, u, v, w);
    let numer_x = dwa * inner[0].x + dwb * inner[1].x + dwc * inner[2].x;
    let numer_y = dwa * inner[0].y + dwb * inner[1].y + dwc * inner[2].y;
    let numer_z = dwa * inner[0].z + dwb * inner[1].z + dwc * inner[2].z;
    Vec3::new(
        (numer_x - blended.x * d_denom) / denom,
        (numer_y - blended.y * d_denom) / denom,
        (numer_z - blended.z * d_denom) / denom,
    )
}

/// Maps triangular-Bezier index triples `(i, j, k)` (with
/// `i + j + k = 3`) to storage positions in
/// `TriangularGregory::boundary_cps`. Returns `None` for the
/// interior `(1, 1, 1)` — that CP comes from the rational blend.
#[inline]
fn boundary_slot(i: usize, j: usize, k: usize) -> Option<usize> {
    match (i, j, k) {
        (3, 0, 0) => Some(0),
        (0, 3, 0) => Some(1),
        (0, 0, 3) => Some(2),
        (2, 1, 0) => Some(3),
        (1, 2, 0) => Some(4),
        (0, 2, 1) => Some(5),
        (0, 1, 2) => Some(6),
        (1, 0, 2) => Some(7),
        (2, 0, 1) => Some(8),
        _ => None,
    }
}

impl TriangularGregory {
    /// Cartesian-parameter wrapper over barycentric evaluation.
    /// The mapping is `(u, v, w) = (u, v, 1 − u − v)`.
    pub fn point_at_uv(&self, u: f64, v: f64) -> Point3 {
        let w = 1.0 - u - v;
        self.point_at_bary(u, v, w)
    }

    /// Cartesian-parameter wrapper for `∂S/∂u`.
    ///
    /// With `w = 1 − u − v`, `∂w/∂u = −1`, so
    /// `∂S/∂u_cart = ∂S/∂u_bary − ∂S/∂w_bary`.
    pub fn partial_u_uv(&self, u: f64, v: f64) -> Vec3 {
        let w = 1.0 - u - v;
        self.partial_bary(u, v, w, 0) - self.partial_bary(u, v, w, 2)
    }

    /// Cartesian-parameter wrapper for `∂S/∂v`.
    pub fn partial_v_uv(&self, u: f64, v: f64) -> Vec3 {
        let w = 1.0 - u - v;
        self.partial_bary(u, v, w, 1) - self.partial_bary(u, v, w, 2)
    }

    fn point_at_bary(&self, u: f64, v: f64, w: f64) -> Point3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..=3 {
            for j in 0..=(3 - i) {
                let k = 3 - i - j;
                let b = bernstein(i, j, k, u, v, w);
                let p = match boundary_slot(i, j, k) {
                    Some(slot) => self.boundary_cps[slot],
                    None => inner_cp(&self.inner_cps, u, v, w),
                };
                x += b * p.x;
                y += b * p.y;
                z += b * p.z;
            }
        }
        Point3::new(x, y, z)
    }

    fn partial_bary(&self, u: f64, v: f64, w: f64, which: usize) -> Vec3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..=3 {
            for j in 0..=(3 - i) {
                let k = 3 - i - j;
                let b = bernstein(i, j, k, u, v, w);
                let db = bernstein_partial(i, j, k, u, v, w, which);
                let (p, dp) = match boundary_slot(i, j, k) {
                    Some(slot) => (self.boundary_cps[slot], Vec3::new(0.0, 0.0, 0.0)),
                    None => (
                        inner_cp(&self.inner_cps, u, v, w),
                        inner_cp_partial(&self.inner_cps, u, v, w, which),
                    ),
                };
                x += db * p.x + b * dp.x;
                y += db * p.y + b * dp.y;
                z += db * p.z + b * dp.z;
            }
        }
        Vec3::new(x, y, z)
    }
}
