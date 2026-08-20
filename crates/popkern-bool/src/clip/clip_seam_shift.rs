//! Seam-shift utilities for UV clipping (split out of `clip_uv.rs`
//! for the 500-line ratchet).

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use super::clip_uv::TWO_PI;

pub(super) fn is_near_seam(curve: &NurbsCurve) -> bool {
    let quarter = std::f64::consts::FRAC_PI_4;
    let seven_quarter = 7.0 * quarter;
    curve
        .control_points
        .iter()
        .any(|p| p.x < quarter || p.x > seven_quarter)
}

pub(super) fn shift_curve_u(curve: &NurbsCurve, shift: f64) -> NurbsCurve {
    let shifted_points = curve
        .control_points
        .iter()
        .map(|p| {
            let mut u = p.x + shift;
            u = u.rem_euclid(TWO_PI);
            Point3::new(u, p.y, p.z)
        })
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: shifted_points,
        weights: curve.weights.clone(),
    }
}
