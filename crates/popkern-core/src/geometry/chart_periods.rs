//! Geometric chart periods: analytic periodicity PLUS the geometric
//! closure of a NURBS chart.
//!
//! A lofted or swept closed profile's wall is one G0-closed NURBS
//! surface: its first and last control rows coincide, so the chart's
//! two u-boundaries are the SAME curve in space, even though
//! `SurfaceGeometry::periodicity` stays `None` (evaluation is not
//! periodic beyond the chart). Consumers that work in chart space —
//! point-in-face classification, the fragment walk's seam-branch
//! identification, and surface-surface MARCHING (which wraps the
//! parameter back into the chart rather than evaluating beyond it) —
//! must treat that closure as a period: without it a chart-wrapping
//! intersection loop is traced as an open arc cut at the seam, its
//! far side never marched (the seed there lies on the consumed
//! polyline), and the boolean strands the loop (PK-440/PK-443).

use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Chart periods for chart-space consumers: analytic periodicity, plus
/// the geometric closure of a NURBS chart (period = the closed
/// direction's knot span). Never use this to evaluate a surface beyond
/// its chart — wrap the parameter into the chart first.
pub fn geometric_chart_periods(geometry: &SurfaceGeometry) -> (Option<f64>, Option<f64>) {
    if let SurfaceGeometry::Nurbs(ns) = geometry {
        let ((u0, u1), (v0, v1)) = ns.parameter_range();
        let (nu, nv) = (ns.num_u, ns.num_v);
        if nu < 2 || nv < 2 || ns.control_points.len() != nu * nv {
            return (None, None);
        }
        let cp = |i: usize, j: usize| ns.control_points[i * nv + j];
        let w = |i: usize, j: usize| ns.weights[i * nv + j];
        let diag = {
            let mut lo = cp(0, 0);
            let mut hi = cp(0, 0);
            for p in &ns.control_points {
                lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
                hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
            }
            lo.distance_to(&hi)
        };
        let eps = 1.0e-9 * (1.0 + diag);
        let u_closed = (0..nv).all(|j| {
            cp(0, j).distance_to(&cp(nu - 1, j)) < eps && (w(0, j) - w(nu - 1, j)).abs() < 1.0e-9
        });
        let v_closed = (0..nu).all(|i| {
            cp(i, 0).distance_to(&cp(i, nv - 1)) < eps && (w(i, 0) - w(i, nv - 1)).abs() < 1.0e-9
        });
        return (u_closed.then_some(u1 - u0), v_closed.then_some(v1 - v0));
    }
    geometry.periodicity()
}

/// Per-direction flags for the marching wrap: true where the surface
/// is analytically periodic OR its NURBS chart closes geometrically.
pub fn marching_periodicity(geometry: &SurfaceGeometry) -> (bool, bool) {
    let (analytic_u, analytic_v) = geometry.is_periodic();
    let (chart_u, chart_v) = geometric_chart_periods(geometry);
    (
        analytic_u || chart_u.is_some(),
        analytic_v || chart_v.is_some(),
    )
}
