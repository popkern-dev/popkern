//! Rigorous bounds of a NURBS sub-patch.
//!
//! A face trimmed from a large NURBS band (a thread runout corner cut
//! from a multi-wrap helical flank) occupies a small window of its
//! chart, but a box built from the whole surface covers the band: the
//! boolean's face-pair prefilter then pairs the corner with every face
//! near the thread (PK-493: two runout corners boxed 3x their trim).
//! Splitting the surface at the window's bounds gives a patch whose
//! control points — with positive weights — contain it.

use super::NurbsSurface;
use crate::geometry::Point3;

impl NurbsSurface {
    /// Control points of the sub-patch over `u_range` x `v_range`
    /// (clamped to the domain): the patch lies in their convex hull
    /// when every weight is positive. `None` when a weight is not
    /// positive or a split fails.
    pub fn sub_patch_hull_points(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Option<Vec<Point3>> {
        if !self.weights.iter().all(|w| *w > 0.0) {
            return None;
        }
        let (u_lo, u_hi) = (self.knots_u[self.degree_u], self.knots_u[self.num_u]);
        let (v_lo, v_hi) = (self.knots_v[self.degree_v], self.knots_v[self.num_v]);
        let clamp = |x: f64, lo: f64, hi: f64| x.max(lo).min(hi);
        let (u0, u1) = (clamp(u_range.0, u_lo, u_hi), clamp(u_range.1, u_lo, u_hi));
        let (v0, v1) = (clamp(v_range.0, v_lo, v_hi), clamp(v_range.1, v_lo, v_hi));
        if !(u0 < u1 && v0 < v1) {
            return None;
        }
        let interior =
            |x: f64, lo: f64, hi: f64| x > lo + 1.0e-12 * (hi - lo) && x < hi - 1.0e-12 * (hi - lo);
        let mut patch = self.clone();
        if interior(u0, u_lo, u_hi) {
            patch = patch.split_at_u(u0).ok()?.1;
        }
        if interior(u1, u0, u_hi) {
            patch = patch.split_at_u(u1).ok()?.0;
        }
        if interior(v0, v_lo, v_hi) {
            patch = patch.split_at_v(v0).ok()?.1;
        }
        if interior(v1, v0, v_hi) {
            patch = patch.split_at_v(v1).ok()?.0;
        }
        Some(patch.control_points)
    }
}

#[cfg(test)]
#[path = "sub_patch_hull_tests.rs"]
mod tests;
