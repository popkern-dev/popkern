//! Where the spokes meet the rim (PK-618).
//!
//! The fit models every sector on the regular characteristic kite, whose
//! spoke meets the rim at a right angle at the side's midpoint. On an
//! irregular opening (setback arcs of very different lengths) a spoke
//! from the fitted centre to a short side's parameter midpoint arrives
//! obliquely (measured sine 0.53 on a steep pentagon), the sector chart
//! shears just inside the rim and the cap bends there (the 6x-disparity
//! cube corner's curvature-weighted fit: 30 / r_min). Splitting each side
//! at its point nearest the fitted centre restores near-perpendicular
//! arrivals: the cube corner's fit reads 2.07 / r_min. The off-centre
//! split is exact under the spoke conditions through the quadratic half
//! map ([`super::rim::HalfMap`]); the fraction is clamped to
//! [`FOOT_CLAMP`], 1 − [`FOOT_CLAMP`] (the half map stays monotone for
//! fractions in (0.25, 0.75); measured on the cube corner, 0.2 folds
//! at the setback corner, 0.3 reads 3.29, 0.35 reads 2.07, 0.4 reads
//! 2.44).

use super::{VariationalCap, VariationalCapSide};
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::healing::curvature::principal_curvatures;

/// Smallest split fraction measured from either end of a side.
pub(super) const FOOT_CLAMP: f64 = 0.35;

/// Interior samples per face parameter for [`interior_curvature`].
const SAMPLES: usize = 30;

/// Where the fit split each side.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapLayout {
    /// Every side at its parameter midpoint.
    Midpoints,
    /// Every side at its point nearest the fitted centre (clamped).
    CentreFeet,
}

/// Every side's parameter midpoint.
pub(super) fn midpoints(sides: &[VariationalCapSide]) -> Vec<f64> {
    sides
        .iter()
        .map(|side| {
            let (t0, t1) = side.curve.parameter_range();
            0.5 * (t0 + t1)
        })
        .collect()
}

/// Each side's parameter nearest `cap`'s centre (sampled at 1/400 of
/// the side), its fraction clamped to [`FOOT_CLAMP`], 1 − [`FOOT_CLAMP`].
pub(super) fn centre_feet(sides: &[VariationalCapSide], cap: &VariationalCap) -> Vec<f64> {
    const STEPS: usize = 400;
    sides
        .iter()
        .map(|side| {
            let (t0, t1) = side.curve.parameter_range();
            let mut best = (f64::INFINITY, 0.5);
            for k in 0..=STEPS {
                let f = k as f64 / STEPS as f64;
                let d = side
                    .curve
                    .point_at(t0 + (t1 - t0) * f)
                    .distance_to(&cap.centre);
                if d < best.0 {
                    best = (d, f);
                }
            }
            t0 + (t1 - t0) * best.1.clamp(FOOT_CLAMP, 1.0 - FOOT_CLAMP)
        })
        .collect()
}

/// The largest sampled principal curvature over the faces' interiors
/// (model units), the reading the blend's fairness bound applies to.
pub(super) fn interior_curvature(cap: &VariationalCap) -> f64 {
    let mut worst = 0.0_f64;
    for face in &cap.faces {
        let surface: &dyn SurfaceEval = face;
        let ((u0, u1), (v0, v1)) = surface.parameter_range();
        for i in 1..SAMPLES {
            for j in 1..SAMPLES {
                let u = u0 + (u1 - u0) * i as f64 / SAMPLES as f64;
                let v = v0 + (v1 - v0) * j as f64 / SAMPLES as f64;
                if let Some((k1, k2)) = principal_curvatures(surface, u, v) {
                    worst = worst.max(k1.abs().max(k2.abs()));
                }
            }
        }
    }
    worst
}
