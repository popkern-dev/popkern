//! Rational Gregory patches — N=4 rectangular and N=3 triangular.
//!
//! Gregory patches (Gregory 1974, Chiyokura-Kimura 1983) resolve
//! the **twist compatibility problem** that arises at corners of
//! ordinary bicubic Bezier patches: at a shared corner, the two
//! incident boundaries each impose a different cross-tangent
//! direction, and an ordinary tensor-product Bezier cannot
//! satisfy both simultaneously. A rational Gregory patch stores
//! *two* candidate inner control points at each corner and blends
//! them rationally by a weight that vanishes on one boundary and
//! grows on the other — giving each boundary its correct cross-
//! tangent plane and thereby G1 continuity with arbitrary
//! neighbours, regardless of valence.
//!
//! # Representation
//!
//! * [`RectangularGregory`] is a bicubic-Bezier-domain patch on
//!   `[0, 1] × [0, 1]` with 4×4 = 16 boundary control points and 4
//!   pairs of inner control points (one pair per corner). Used
//!   for N=4 corner resolution and for three-quad-split N=3.
//! * [`TriangularGregory`] is a cubic-triangular-Bezier-domain
//!   patch on the barycentric simplex `{(u, v, w) : u + v + w = 1,
//!   u, v, w ≥ 0}` with 9 boundary control points and 3 pairs of
//!   inner control points. Used for native N=3 corner resolution.
//!
//! The inner-control-point *pair* at each corner is selected from
//! two candidates (one matching each incident boundary's cross-
//! tangent) and rationally blended at evaluation time. The
//! blending weights are C∞-smooth in the interior and have a
//! single removable singularity at the corner itself; the
//! evaluator handles the corner as a limiting value.
//!
//! # References
//!
//! * Gregory, J. A. (1974). "Smooth Interpolation Without Twist
//!   Constraints." *Computer Aided Geometric Design*.
//! * Chiyokura, H., & Kimura, F. (1983). "Design of Solids with
//!   Free-form Surfaces." *Computer Graphics*, 17(3).
//! * Farin, G. (1993). "Curves and Surfaces for CAGD", Ch. 18.
//! * Piegl, L., & Tiller, W. (1997). "The NURBS Book", §9.

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

/// A triangular (N=3) rational Gregory patch on the barycentric
/// simplex `{(u, v, w) : u + v + w = 1, u, v, w ≥ 0}`.
///
/// Stores 9 **boundary** control points (3 corners + 6 edge
/// midpoints at the thirds of each edge) and 3 pairs of **inner**
/// control points — one pair per corner. Like the rectangular
/// form, each inner pair is rationally blended to give the two
/// boundaries sharing that corner their correct cross-tangents
/// independently.
///
/// The 9 boundary CPs are stored in the standard triangular
/// Bezier order (indices `(i, j, k)` with `i + j + k = 3`):
/// `[300, 030, 003, 210, 120, 021, 012, 102, 201]` —
/// corner0, corner1, corner2, then the two CPs along each edge
/// (edge 01, edge 12, edge 20).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriangularGregory {
    /// The 9 boundary control points in triangular-Bezier index
    /// order: three corners followed by six edge-thirds CPs.
    pub boundary_cps: [Point3; 9],
    /// Three pairs of inner control points — one pair per corner
    /// in corner order `(0, 1, 2)`. Each pair is
    /// `(P_clockwise, P_counterclockwise)` giving the cross-
    /// tangents for the two boundaries meeting at that corner.
    pub inner_cps: [(Point3, Point3); 3],
}

/// A Gregory patch of either rectangular (N=4) or triangular
/// (N=3) topology.
///
/// Boxed because the fixed-size control-point arrays on both
/// variants are large (~400–600 bytes) and including them inline
/// would inflate every `SurfaceGeometry` enum instance to match;
/// Gregory patches are rare (only at corners) so the one extra
/// pointer indirection per evaluation is free compared to the
/// per-NURBS-face cache cost of a bloated enum.
///
/// Use [`SurfaceGeometry::Gregory`](crate::geometry::SurfaceGeometry::Gregory)
/// to attach to a face.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GregoryPatch {
    /// N=4 rectangular patch.
    Rectangular(Box<RectangularGregory>),
    /// N=3 triangular patch.
    Triangular(Box<TriangularGregory>),
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

/// Rationally-blended inner control points at parameter `(u, v)`
/// for a rectangular Gregory patch.
///
/// The four corners are ordered (0,0) → (0,3) → (3,3) → (3,0)
/// counter-clockwise, matching [`RectangularGregory::inner_cps`].
/// The blending weight at each corner vanishes on the boundary
/// opposite its cross-tangent direction — so the corner's
/// contribution is wholly `P_u` on the v-boundary through it and
/// wholly `P_v` on the u-boundary through it. The combined
/// denominator is a smooth function of `(u, v)` strictly positive
/// on the open domain `(0, 1) × (0, 1)`; at an exact corner the
/// limit is taken (arbitrary tie-break — the Bernstein factor at
/// that corner is zero, so the inner CP does not contribute).
#[inline]
fn rect_inner_cp(inner_cps: &[(Point3, Point3); 4], i: usize, j: usize, u: f64, v: f64) -> Point3 {
    // Map (i, j) ∈ {(1,1), (1,2), (2,1), (2,2)} → corner index k.
    // Corner ordering (k): (0,0)=0, (0,3)=1, (3,3)=2, (3,0)=3.
    // The inner CP at Bezier index (i, j) sits one step inside
    // the corresponding corner.
    let (k, wu, wv) = match (i, j) {
        (1, 1) => (0, u, v),       // corner (0,0): wu→0 on u=0 edge, wv→0 on v=0 edge.
        (1, 2) => (1, u, 1.0 - v), // corner (0,3): wu→0 on u=0, wv→0 on v=1.
        (2, 2) => (2, 1.0 - u, 1.0 - v), // corner (3,3): wu→0 on u=1, wv→0 on v=1.
        (2, 1) => (3, 1.0 - u, v), // corner (3,0): wu→0 on u=1, wv→0 on v=0.
        _ => unreachable!("rect_inner_cp called with non-inner index"),
    };
    let (p_u, p_v) = inner_cps[k];
    let denom = wu + wv;
    if denom < 1e-15 {
        // Exact corner: Bernstein factor vanishes; tie-break.
        return Point3::new(
            0.5 * (p_u.x + p_v.x),
            0.5 * (p_u.y + p_v.y),
            0.5 * (p_u.z + p_v.z),
        );
    }
    // Blend: P(u, v) = (wu·P_v + wv·P_u) / (wu + wv).
    // P_u's weight is `wv` because `P_u` is the candidate that
    // matches the u-direction boundary's cross-tangent, which
    // dominates as we move along `v` toward the corner (i.e. as
    // `wu` shrinks). The Chiyokura-Kimura convention.
    let c_u = wv / denom;
    let c_v = wu / denom;
    Point3::new(
        c_u * p_u.x + c_v * p_v.x,
        c_u * p_u.y + c_v * p_v.y,
        c_u * p_u.z + c_v * p_v.z,
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
                let p =
                    if (i, j) == (1, 1) || (i, j) == (1, 2) || (i, j) == (2, 1) || (i, j) == (2, 2)
                    {
                        rect_inner_cp(&self.inner_cps, i, j, u, v)
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
    /// The inner control points depend on `(u, v)` via the
    /// rational blend, so the full derivative is
    /// `Σ_ij ∂B_i^3/∂u · B_j^3 · P_ij + B_i^3 · B_j^3 · ∂P_ij/∂u`.
    /// The second term vanishes on all four boundaries (where one
    /// of the Bernstein factors is zero at the inner indices) —
    /// this is exactly why the rational construction delivers
    /// G1 there without further adjustment. Inside the domain
    /// we compute `∂P_ij/∂u` from the quotient rule.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in 0..4 {
            let bu_d = bernstein_cubic_deriv(i, u);
            let bu = bernstein_cubic(i, u);
            for j in 0..4 {
                let bv = bernstein_cubic(j, v);
                let (p, dp_du) =
                    if (i, j) == (1, 1) || (i, j) == (1, 2) || (i, j) == (2, 1) || (i, j) == (2, 2)
                    {
                        (
                            rect_inner_cp(&self.inner_cps, i, j, u, v),
                            rect_inner_cp_partial_u(&self.inner_cps, i, j, u, v),
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
                let (p, dp_dv) =
                    if (i, j) == (1, 1) || (i, j) == (1, 2) || (i, j) == (2, 1) || (i, j) == (2, 2)
                    {
                        (
                            rect_inner_cp(&self.inner_cps, i, j, u, v),
                            rect_inner_cp_partial_v(&self.inner_cps, i, j, u, v),
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

/// ∂/∂u of the rational inner-CP blend. See [`rect_inner_cp`].
///
/// Let `D = wu + wv`. Then
/// `P = (wv·P_u + wu·P_v) / D`.
/// With `(wu, wv) = (α(u), β(v))` (product of at most a `u`-
/// only factor and a `v`-only factor per the corner table), we
/// get `∂wu/∂u = α'(u)`, `∂wv/∂u = 0`, so
/// `∂P/∂u = [α'·P_v · D − (wv·P_u + wu·P_v)·α'] / D²`
///       `= α' · (P_v − P)/D`.
fn rect_inner_cp_partial_u(
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
    // ∂P/∂u = du_sign · (wv · (P_v - P_u)) / (wu + wv)² ... wait
    // Re-derive: P = (wv·P_u + wu·P_v) / (wu + wv).
    // dP/d wu = (P_v·D - (wv·P_u + wu·P_v)) / D²
    //        = (P_v · (wu+wv) - wv·P_u - wu·P_v) / D²
    //        = wv · (P_v - P_u) / D².
    // dwu/du = du_sign.
    let factor = du_sign * wv / d2;
    Vec3::new(
        factor * (p_v.x - p_u.x),
        factor * (p_v.y - p_u.y),
        factor * (p_v.z - p_u.z),
    )
}

/// ∂/∂v of the rational inner-CP blend. Symmetric to
/// [`rect_inner_cp_partial_u`].
fn rect_inner_cp_partial_v(
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
    // dP/d wv = wu · (P_u - P_v) / D²
    // dwv/dv = dv_sign.
    let factor = dv_sign * wu / d2;
    Vec3::new(
        factor * (p_u.x - p_v.x),
        factor * (p_u.y - p_v.y),
        factor * (p_u.z - p_v.z),
    )
}

impl GregoryPatch {
    /// Evaluate at the given parameters. `(u, v)` is the Cartesian
    /// parameter for [`Rectangular`](GregoryPatch::Rectangular); for
    /// [`Triangular`](GregoryPatch::Triangular), `(u, v)` map into
    /// barycentric coords via `(u, v, 1 − u − v)` and must satisfy
    /// `u ≥ 0`, `v ≥ 0`, `u + v ≤ 1`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        match self {
            GregoryPatch::Rectangular(r) => r.point_at(u, v),
            GregoryPatch::Triangular(t) => t.point_at_uv(u, v),
        }
    }

    /// Partial derivative with respect to `u`. See
    /// [`point_at`](Self::point_at) for the parameter convention.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        match self {
            GregoryPatch::Rectangular(r) => r.partial_u(u, v),
            GregoryPatch::Triangular(t) => t.partial_u_uv(u, v),
        }
    }

    /// Partial derivative with respect to `v`. See
    /// [`point_at`](Self::point_at) for the parameter convention.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        match self {
            GregoryPatch::Rectangular(r) => r.partial_v(u, v),
            GregoryPatch::Triangular(t) => t.partial_v_uv(u, v),
        }
    }
}

// Triangular Gregory stubs — evaluator filled in a follow-up
// commit in Phase 1. For now return boundary-corner placeholders
// so callers that construct `TriangularGregory` don't fault; the
// construction path is not yet wired in.

impl TriangularGregory {
    /// Cartesian-parameter wrapper over barycentric evaluation.
    /// The mapping is `(u, v, w) = (u, v, 1 − u − v)`.
    pub fn point_at_uv(&self, u: f64, v: f64) -> Point3 {
        let w = 1.0 - u - v;
        self.point_at_bary(u, v, w)
    }

    /// Cartesian-parameter wrapper for `∂S/∂u` on the triangular
    /// patch. Equivalent to `∂S/∂u − ∂S/∂w` when the patch is
    /// expressed in barycentrics.
    pub fn partial_u_uv(&self, u: f64, v: f64) -> Vec3 {
        let w = 1.0 - u - v;
        self.partial_bary(u, v, w, 0) - self.partial_bary(u, v, w, 2)
    }

    /// Cartesian-parameter wrapper for `∂S/∂v`.
    pub fn partial_v_uv(&self, u: f64, v: f64) -> Vec3 {
        let w = 1.0 - u - v;
        self.partial_bary(u, v, w, 1) - self.partial_bary(u, v, w, 2)
    }

    /// Barycentric evaluator — filled in Phase 1 C2.
    fn point_at_bary(&self, _u: f64, _v: f64, _w: f64) -> Point3 {
        self.boundary_cps[0]
    }

    /// Barycentric partial placeholder — filled in Phase 1 C2.
    fn partial_bary(&self, _u: f64, _v: f64, _w: f64, _which: usize) -> Vec3 {
        Vec3::new(0.0, 0.0, 0.0)
    }
}

#[cfg(test)]
#[path = "gregory_tests.rs"]
mod tests;
