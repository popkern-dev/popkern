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

/// Interior knot values where a NURBS chart GEOMETRICALLY creases,
/// per direction, in ascending order. Candidates are knots with
/// multiplicity at or above the degree (parametric C0); each is then
/// verified by the one-sided tangent angle across the iso line at
/// several stations — a rational circle profile carries
/// full-multiplicity knots at its arc junctions while staying G1
/// (splitting there diced every circular-section sweep record and
/// stranded the spring grind), whereas a faceted wall turns by the
/// facet angle. Consumers that fit smooth curves across chart runs
/// must split at these (PK-442: a record fitted across a prism crease
/// rounded the true corner by 1e-2). Analytic surfaces and smooth
/// charts return empty.
pub fn c0_interior_knots(geometry: &SurfaceGeometry) -> (Vec<f64>, Vec<f64>) {
    let SurfaceGeometry::Nurbs(ns) = geometry else {
        return (Vec::new(), Vec::new());
    };
    use crate::geometry::surface::SurfaceEval as _;
    let ((u0, u1), (v0, v1)) = ns.parameter_range();
    let scan = |knots: &[f64], degree: usize| -> Vec<f64> {
        if knots.len() < 2 || degree == 0 {
            return Vec::new();
        }
        let lo = knots[0];
        let hi = knots[knots.len() - 1];
        let mut out: Vec<f64> = Vec::new();
        let mut i = 0;
        while i < knots.len() {
            let v = knots[i];
            let mut m = 1;
            while i + m < knots.len() && (knots[i + m] - v).abs() < 1.0e-12 {
                m += 1;
            }
            if v > lo + 1.0e-12 && v < hi - 1.0e-12 && m >= degree {
                out.push(v);
            }
            i += m;
        }
        out
    };
    // One-sided tangent angle across the iso line, worst of several
    // stations; below the gate the junction is geometrically smooth.
    const CREASE_COS: f64 = 0.999_999; // ~1.4e-3 rad
    let geometric = |candidates: Vec<f64>, is_u: bool| -> Vec<f64> {
        candidates
            .into_iter()
            .filter(|&k| {
                let span = if is_u { u1 - u0 } else { v1 - v0 };
                let h = span * 1.0e-5;
                for station in [0.17_f64, 0.5, 0.83] {
                    let w = if is_u {
                        v0 + (v1 - v0) * station
                    } else {
                        u0 + (u1 - u0) * station
                    };
                    let at = |t: f64| {
                        if is_u {
                            geometry.point_at(t, w)
                        } else {
                            geometry.point_at(w, t)
                        }
                    };
                    let before = at(k - h) - at(k - 2.0 * h);
                    let after = at(k + 2.0 * h) - at(k + h);
                    let (lb, la) = (before.length(), after.length());
                    if lb <= 0.0 || la <= 0.0 {
                        continue;
                    }
                    if before.dot(&after) / (lb * la) < CREASE_COS {
                        return true;
                    }
                }
                false
            })
            .collect()
    };
    let cu = geometric(scan(&ns.knots_u, ns.degree_u), true);
    let cv = geometric(scan(&ns.knots_v, ns.degree_v), false);
    (cu, cv)
}
