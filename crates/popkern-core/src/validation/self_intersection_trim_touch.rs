//! An intersection curve touching a face's trim from outside (PK-1153).
//!
//! Two faces' surfaces can meet in a curve that reaches one face's
//! trim at a single point and turns back: a ball resting on a box's
//! top edge meets the top plane in a circle tangent to that edge from
//! beyond it. The sample at the touch reads `Boundary` on the top face
//! and `Interior` on the ball, which the overlap test took for a
//! crossing, so a body whose faces only touch there (two shells two
//! tolerances apart, or one general body glued at the point) failed
//! `validate_all`. A crossing enters the face beside the sample; a
//! touch from outside leaves it on both sides.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{BRep, FaceHandle};

use super::self_intersection_uv_containment::TrimContainment;

/// The steps either side of the sample, as fractions of the curve's
/// parameter span, tried in turn: the finest whose sides both leave
/// the trim's boundary band decides (a touch at millimetre scale
/// stays inside the band a sixty-fourth of a small circle away).
const SIDE_STEPS: [f64; 4] = [1.0 / 64.0, 1.0 / 32.0, 1.0 / 16.0, 1.0 / 8.0];

impl BRep {
    /// Whether `curve` at fraction `alpha` of its span lies outside
    /// `face`'s trim either side of the sample, at the finest step
    /// clear of the boundary band and inside it at none: a touch of the
    /// trim from outside, not a crossing. An open curve's end has one side and is
    /// never a touch; a closed curve's sides wrap across its ends.
    pub(super) fn curve_touches_trim_from_outside(
        &self,
        face: FaceHandle,
        surface: &SurfaceGeometry,
        curve: &NurbsCurve,
        alpha: f64,
        linear_tol: f64,
    ) -> bool {
        let (t0, t1) = curve.parameter_range();
        let span = t1 - t0;
        if !span.is_finite() || span <= 0.0 {
            return false;
        }
        let closed = curve.point_at(t0).distance_to(&curve.point_at(t1)) <= linear_tol;
        let side = |a: f64| -> Option<f64> {
            if (0.0..=1.0).contains(&a) {
                Some(a)
            } else if closed {
                Some(a.rem_euclid(1.0))
            } else {
                None
            }
        };
        for step in SIDE_STEPS {
            let mut outside = 0;
            for a in [alpha - step, alpha + step] {
                let Some(a) = side(a) else {
                    return false;
                };
                let p = curve.point_at(t0 + span * a);
                match self.face_trim_classify_point(face, surface, &p, linear_tol) {
                    Some(TrimContainment::Outside) => outside += 1,
                    Some(TrimContainment::Boundary) => {}
                    // The curve enters the face beside the sample.
                    _ => return false,
                }
            }
            if outside == 2 {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
#[path = "self_intersection_trim_touch_tests.rs"]
mod tests;
