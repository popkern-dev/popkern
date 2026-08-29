//! Side-curve construction for exact-boundary sector solves.
//!
//! Nested inside `solve` (shares its conventions); extracted to hold
//! the parent under the 500-line module budget.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};

/// Exactly re-represent a prepared spoke curve in `template`'s basis
/// (degree elevation + knot insertion, both exact operations).
/// Returns `None` when the template basis cannot host the spoke —
/// the caller falls back to a straight-line side.
pub(super) fn spoke_in_basis(spoke: &NurbsCurve, template: &NurbsCurve) -> Option<NurbsCurve> {
    if template.degree < spoke.degree {
        return None;
    }
    let elevated = spoke.elevate_degree_to(template.degree).ok()?;
    let (lo, hi) = template.parameter_range();
    let (spoke_lo, spoke_hi) = elevated.parameter_range();
    if (lo - spoke_lo).abs() > 1.0e-12 || (hi - spoke_hi).abs() > 1.0e-12 {
        return None;
    }
    let mut missing = Vec::new();
    let mut elevated_interior: Vec<f64> =
        elevated.knots[elevated.degree + 1..elevated.knots.len() - elevated.degree - 1].to_vec();
    for &knot in &template.knots[template.degree + 1..template.knots.len() - template.degree - 1] {
        if let Some(position) = elevated_interior
            .iter()
            .position(|&existing| (existing - knot).abs() < 1.0e-12)
        {
            elevated_interior.remove(position);
        } else {
            missing.push(knot);
        }
    }
    if !elevated_interior.is_empty() {
        return None;
    }
    let refined = if missing.is_empty() {
        elevated
    } else {
        elevated.insert_knots(&missing).ok()?
    };
    (refined.control_points.len() == template.control_points.len()).then_some(refined)
}

pub(super) fn line_curve_like(template: &NurbsCurve, start: Point3, end: Point3) -> NurbsCurve {
    let (lo, hi) = template.parameter_range();
    let span = hi - lo;
    let degree = template.degree;
    let cps = (0..template.control_points.len())
        .map(|i| {
            let t = if degree == 0 {
                i as f64 / (template.control_points.len() - 1) as f64
            } else {
                let greville =
                    template.knots[i + 1..=i + degree].iter().sum::<f64>() / degree as f64;
                (greville - lo) / span
            };
            super::lerp(start, end, t.clamp(0.0, 1.0))
        })
        .collect();
    NurbsCurve {
        degree,
        knots: template.knots.clone(),
        control_points: cps,
        weights: template.weights.clone(),
    }
}
