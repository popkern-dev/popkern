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
/// The parametric band of a sharp trim reading.
const SHARP_BAND: f64 = 1e-12;

impl BRep {
    /// Whether `point` lies within `reach` of one of `face`'s boundary
    /// curves (its outer loop and holes), measured in space: what a
    /// uv in the trim test's boundary band (ten tolerances wide) needs
    /// before it counts as on the face to the tolerance.
    pub fn point_near_face_boundary(&self, face: FaceHandle, point: &Point3, reach: f64) -> bool {
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

    /// `class`, `face`'s trim reading of `point`, confirmed in space: a
    /// `Boundary` (a uv in the trim test's band, ten tolerances wide,
    /// more where the chart is stretched) farther than the linear reach
    /// from every boundary curve is no place on the face's boundary (one
    /// place within the tolerance, PK-1131: the glue reads the same
    /// reach), and reads `Outside` where the trim reads it sharply
    /// outside. A cone's
    /// face touching a frustum's side 2.4e-5 from its cap's rim met
    /// the cap's plane 1.7e-6 outside the cap's rim, in the band, and the
    /// glued body failed as self-intersecting (CF-E 150 at scale 1e-3;
    /// faces touch only within the tolerance of each boundary, PK-1154).
    pub(super) fn trim_reading_in_space(
        &self,
        (face, surface): (FaceHandle, &SurfaceGeometry),
        point: &Point3,
        class: Option<TrimContainment>,
        linear_tol: f64,
    ) -> Option<TrimContainment> {
        if class != Some(TrimContainment::Boundary)
            || self.point_near_face_boundary(
                face,
                point,
                linear_tol * (1.0 + crate::tolerance::LINEAR_ROUNDING),
            )
        {
            return class;
        }
        let (u, v) =
            crate::geometry::uv_projection::project_point_to_uv_multisheet(surface, point)?;
        let sharp = self.face_trim_classify_uv_band(
            face,
            surface,
            super::self_intersection_uv_containment::UvPoint { u, v },
            SHARP_BAND,
        );
        // Only off the face: a point sharply inside but past the reach
        // from the boundary stays on the band, where a curve grazing the
        // boundary from inside (a dimple's rim touching a plate's edge at
        // milli scale, 1.2e-6 in) is the touch, not a crossing.
        match sharp {
            Some(TrimContainment::Outside) => Some(TrimContainment::Outside),
            _ => class,
        }
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

impl BRep {
    /// Whether `point`, where `other` meets `face` inside its trim's
    /// boundary band, is a point of `other` touching a face across
    /// `face`'s boundary there: within the linear reach of that face's
    /// surface and in the band of a definite touch of the two
    /// (`points_in_touch_band`, the one rule for a touch decided by
    /// distance). An edge resting on a face inside a touch's patch
    /// (CC-E 17 and 9 at scale 1e-3: B's rim 5e-7 off A's cone, B's
    /// cone touching A's cone 1e-5 away) has the surfaces of three
    /// faces meeting within the tolerance across the patch, and the
    /// curve of A's cone and B's base plane runs along B's rim there,
    /// reading `Interior` on A's cone and `Boundary` on B's base plane
    /// up to 2.6e-6 inside its rim: it is the touch of the cones, not a
    /// crossing.
    pub(super) fn meets_across_boundary_in_touch_band(
        &self,
        face: FaceHandle,
        (other, other_surface): (FaceHandle, &SurfaceGeometry),
        point: &Point3,
        linear_tol: f64,
    ) -> bool {
        use crate::geometry::intersection::touch_ring::{points_in_touch_band, window_extent};
        let reach = self.tolerance.linear_reach().max(linear_tol);
        let band = 10.0 * linear_tol;
        let Some(face_data) = self.faces.get(face) else {
            return false;
        };
        let other_window =
            crate::validation::face_uv_window::face_trim_uv_window(self, other, other_surface);
        let mut across = Vec::new();
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
                let near = self
                    .edges
                    .get(he.edge)
                    .and_then(|e| e.curve.as_ref())
                    .is_some_and(|c| {
                        let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                        (&c.geometry as &dyn CurveEval)
                            .closest_point(point)
                            .is_some_and(|t| {
                                c.geometry.point_at(t.clamp(lo, hi)).distance_to(point) <= band
                            })
                    });
                if near {
                    for fin in self.radial_fins(he.edge) {
                        let h = self
                            .half_edges
                            .get(fin)
                            .and_then(|f| self.loops.get(f.loop_));
                        if let Some(h) = h.map(|l| l.face)
                            && h != face
                            && h != other
                            && !across.contains(&h)
                        {
                            across.push(h);
                        }
                    }
                }
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
        }
        across.into_iter().any(|h| {
            let Some(surface) = self.faces.get(h).and_then(|f| f.surface.as_ref()) else {
                return false;
            };
            let s = &surface.geometry;
            let window = crate::validation::face_uv_window::face_trim_uv_window(self, h, s);
            let extent = window_extent(s, window).max(window_extent(other_surface, other_window));
            points_in_touch_band(
                std::slice::from_ref(point),
                (other_surface, other_window),
                (s, window),
                (linear_tol, reach, extent),
            )
        })
    }
}

/// The samples of a coincidence region that lie within `reach` of
/// `other`. The pair scan finds coincidence at twice the tolerance, so
/// a region it reports may be two faces up to that far apart: faces
/// between one and two tolerances apart are distinct (PK-1131), and
/// two blocks set 1.5e-6 apart, turned, failed validation as one face
/// lying on the other.
pub(super) fn region_samples_within(
    other: &SurfaceGeometry,
    points: &[Point3],
    reach: f64,
) -> Vec<Point3> {
    use crate::geometry::surface::SurfaceEval;
    points
        .iter()
        .filter(|p| {
            crate::geometry::uv_projection::project_point_to_uv(other, p)
                .is_some_and(|(u, v)| other.point_at(u, v).distance_to(p) <= reach)
        })
        .copied()
        .collect()
}

#[cfg(test)]
#[path = "self_intersection_trim_touch_tests.rs"]
mod tests;
