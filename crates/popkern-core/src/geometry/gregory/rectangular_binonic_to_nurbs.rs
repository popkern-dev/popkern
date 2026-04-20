//! Exact bidegree-(9, 9) NURBS conversion for
//! [`RectangularGregoryBiNonic`].
//!
//! A [`RectangularGregoryBiNonic`] is a **plain** (non-rational)
//! tensor-product Bezier patch of bidegree (9, 9). The conversion
//! to a clamped [`NurbsSurface`] is therefore trivial and exact —
//! no rational-denominator clearing, no ε-regularisation — with
//! the same 10×10 control-point grid, all weights `1.0`, and a
//! standard clamped knot vector `[0;10, 1;10]` in both directions.
//!
//! # Identity
//!
//! A single-interval clamped B-spline of degree `n` on `[0, 1]`
//! with `n + 1` control points is exactly the Bernstein-basis
//! evaluation
//!
//! ```text
//!   S(u, v) = Σ_{i, j} P_{ij} · B_i^9(u) · B_j^9(v)
//! ```
//!
//! which is the definition of [`RectangularGregoryBiNonic::point_at`].
//! Hence [`RectangularGregoryBiNonic::to_nurbs_exact`] returns a
//! NURBS that equals the patch pointwise to floating-point
//! precision — no approximation involved.
//!
//! # References
//!
//! * Piegl, L., & Tiller, W. (1997). "The NURBS Book", §1.4
//!   (Bernstein-Bezier ↔ single-interval clamped B-spline).

use super::rectangular_binonic::RectangularGregoryBiNonic;
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

impl RectangularGregoryBiNonic {
    /// Convert to an exact bidegree-(9, 9) clamped [`NurbsSurface`]
    /// on `[0, 1]²`. All weights are `1.0` (non-rational), and
    /// the returned NURBS matches [`Self::point_at`] to floating-
    /// point precision at every `(u, v)`.
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        let mut control_points: Vec<Point3> = Vec::with_capacity(100);
        for row in &self.grid {
            for cp in row {
                control_points.push(*cp);
            }
        }
        let weights = vec![1.0; 100];

        let mut knots = Vec::with_capacity(20);
        knots.extend(std::iter::repeat_n(0.0, 10));
        knots.extend(std::iter::repeat_n(1.0, 10));

        NurbsSurface {
            degree_u: 9,
            degree_v: 9,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            num_u: 10,
            num_v: 10,
            weights,
        }
    }
}

#[cfg(test)]
#[path = "rectangular_binonic_to_nurbs_tests.rs"]
mod tests;
