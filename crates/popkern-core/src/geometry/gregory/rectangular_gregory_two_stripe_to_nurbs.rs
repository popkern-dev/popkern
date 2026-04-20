//! Exact bidegree-(7, 7) Bezier → NURBS conversion for
//! [`RectangularGregoryTwoStripe`]. Trivial compared to the
//! Gregory variants — no rational blend to clear, no
//! denominator to elevate. The 8×8 plain-Bezier control grid
//! IS a clamped bidegree-(7, 7) NURBS with unit weights.

use super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe;
use crate::geometry::nurbs_surface::NurbsSurface;

impl RectangularGregoryTwoStripe {
    /// Convert to a clamped bidegree-(7, 7) NURBS on `[0, 1]²`.
    /// Exact — plain Bezier is a unit-weight NURBS.
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        let mut control_points = Vec::with_capacity(64);
        for row in self.grid.iter() {
            for p in row.iter() {
                control_points.push(*p);
            }
        }
        let weights = vec![1.0; 64];
        let mut knots = Vec::with_capacity(16);
        knots.extend(std::iter::repeat_n(0.0, 8));
        knots.extend(std::iter::repeat_n(1.0, 8));
        NurbsSurface {
            degree_u: 7,
            degree_v: 7,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            num_u: 8,
            num_v: 8,
            weights,
        }
    }
}
