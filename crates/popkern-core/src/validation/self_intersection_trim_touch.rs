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
//!
//! The trim test's boundary band is ten linear tolerances wide, so a
//! sample on both faces' boundary bands can lie up to the band off
//! either face: a plate filleted ten tolerances short of its thickness
//! keeps a front strip that wide, and its blend and bottom face,
//! sharing no vertex, read Boundary on both along the strip (PK-1154).
//! Faces touch only where the sample is within the tolerance of each
//! face's boundary curves.

use crate::geometry::Point3;
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

/// Stations along a boundary edge for the coarse distance search.
const EDGE_STATIONS: usize = 32;
/// Ternary-search rounds refining the nearest station.
const EDGE_ROUNDS: usize = 40;

impl BRep {
    /// Whether `point` lies within `reach` of one of `face`'s boundary
    /// curves (its outer loop and holes).
    pub(super) fn point_near_face_boundary(
        &self,
        face: FaceHandle,
        point: &Point3,
        reach: f64,
    ) -> bool {
        let Some(face_data) = self.faces.get(face) else {
            return true;
        };
        for l in std::iter::once(face_data.outer_loop).chain(face_data.inner_loops.iter().copied())
        {
            let Some(start) = self.loops.get(l).map(|l| l.half_edge) else {
                continue;
            };
            let mut cursor = start;
            for _ in 0..=self.half_edges.len() {
                let Some(he) = self.half_edges.get(cursor) else {
                    break;
                };
                let curve = self.edges.get(he.edge).and_then(|e| e.curve.as_ref());
                // An edge with no curve is the segment between its
                // vertices.
                let ends = (
                    self.vertices.get(he.origin).map(|v| v.position),
                    self.half_edges
                        .get(he.next)
                        .and_then(|n| self.vertices.get(n.origin))
                        .map(|v| v.position),
                );
                let at = |s: f64| -> Option<Point3> {
                    match (curve, ends) {
                        (Some(c), _) => {
                            Some(c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * s))
                        }
                        (None, (Some(a), Some(b))) => Some(a + (b - a) * s),
                        _ => None,
                    }
                };
                if at(0.0).is_some() {
                    let dist = |s: f64| at(s).map_or(f64::INFINITY, |p| p.distance_to(point));
                    let step = 1.0 / EDGE_STATIONS as f64;
                    let best = (0..=EDGE_STATIONS)
                        .map(|i| i as f64 * step)
                        .min_by(|a, b| dist(*a).total_cmp(&dist(*b)))
                        .unwrap_or(0.0);
                    let (mut lo, mut hi) = ((best - step).max(0.0), (best + step).min(1.0));
                    for _ in 0..EDGE_ROUNDS {
                        let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
                        if dist(m1) <= dist(m2) {
                            hi = m2;
                        } else {
                            lo = m1;
                        }
                    }
                    if dist(0.5 * (lo + hi)) <= reach {
                        return true;
                    }
                }
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
        }
        false
    }

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
