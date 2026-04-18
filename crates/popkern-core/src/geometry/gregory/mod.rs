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
//!   **pairs** of inner control points (one pair per corner — the
//!   Chiyokura-Kimura form). Used for N=4 corner resolution and
//!   for three-quad-split N=3.
//! * [`TriangularGregory`] is a cubic-triangular-Bezier-domain
//!   patch on the barycentric simplex `{(u, v, w) : u + v + w = 1,
//!   u, v, w ≥ 0}` with 9 boundary control points and 3 inner
//!   **candidates** (one per corner — the Longhi form). Used for
//!   native N=3 corner resolution.
//!
//! In both forms the inner control point at evaluation is a
//! rational blend — two-way for rectangular (per corner), three-
//! way for triangular (global, weighted by barycentric products).
//! The blending weights are C∞-smooth on the open interior and
//! have removable singularities at the corners where the
//! Bernstein basis of the interior control point also vanishes,
//! so the inner CP does not contribute and the limit is chosen
//! as a tie-break without affecting the patch value.
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

mod build;
mod rectangular;
mod rectangular_deg5;
mod rectangular_deg5_to_nurbs;
mod rectangular_to_nurbs;
mod triangular;

pub use build::GregoryBuildError;
pub use rectangular::RectangularGregory;
pub use rectangular_deg5::RectangularGregoryDeg5;
pub use triangular::TriangularGregory;

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
    /// N=4 rectangular cubic patch (bidegree 3×3).
    Rectangular(Box<RectangularGregory>),
    /// N=4 rectangular quintic patch (bidegree 5×5). Same
    /// Chiyokura-Kimura 4-inner-CP-pair structure as
    /// [`Rectangular`](GregoryPatch::Rectangular); higher
    /// boundary degree closes the 2% cubic-Bezier arc
    /// approximation error inherent at the lower degree.
    RectangularDeg5(Box<RectangularGregoryDeg5>),
    /// N=3 triangular patch.
    Triangular(Box<TriangularGregory>),
}

impl GregoryPatch {
    /// Evaluate at the given parameters. `(u, v)` is the Cartesian
    /// parameter for [`Rectangular`](GregoryPatch::Rectangular) and
    /// [`RectangularDeg5`](GregoryPatch::RectangularDeg5); for
    /// [`Triangular`](GregoryPatch::Triangular), `(u, v)` map into
    /// barycentric coords via `(u, v, 1 − u − v)` and must satisfy
    /// `u ≥ 0`, `v ≥ 0`, `u + v ≤ 1`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        match self {
            GregoryPatch::Rectangular(r) => r.point_at(u, v),
            GregoryPatch::RectangularDeg5(r) => r.point_at(u, v),
            GregoryPatch::Triangular(t) => t.point_at_uv(u, v),
        }
    }

    /// Partial derivative with respect to `u`. See
    /// [`point_at`](Self::point_at) for the parameter convention.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        match self {
            GregoryPatch::Rectangular(r) => r.partial_u(u, v),
            GregoryPatch::RectangularDeg5(r) => r.partial_u(u, v),
            GregoryPatch::Triangular(t) => t.partial_u_uv(u, v),
        }
    }

    /// Partial derivative with respect to `v`. See
    /// [`point_at`](Self::point_at) for the parameter convention.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        match self {
            GregoryPatch::Rectangular(r) => r.partial_v(u, v),
            GregoryPatch::RectangularDeg5(r) => r.partial_v(u, v),
            GregoryPatch::Triangular(t) => t.partial_v_uv(u, v),
        }
    }
}

#[cfg(test)]
#[path = "../gregory_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "build_tests.rs"]
mod build_tests;

#[cfg(test)]
#[path = "rectangular_deg5_tests.rs"]
mod rectangular_deg5_tests;

#[cfg(test)]
#[path = "rectangular_to_nurbs_tests.rs"]
mod rectangular_to_nurbs_tests;
