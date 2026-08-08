//! [`BRep::interior_overlap_samples`] — the validator's Interior
//! overlap classification exposed as localized samples for repair and
//! diagnostic consumers (PK-251).

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{BRep, FaceHandle};

use super::self_intersection_uv_containment::{TrimContainment, UvPoint};
use super::{intersection_sample_points, surface_uses_iterative_projection};

/// One Interior-classified overlap sample between two faces (see
/// [`BRep::interior_overlap_samples`]). UVs are in the caller's
/// argument order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InteriorOverlapSample {
    /// 3D location of the sample on the intersection.
    pub point: Point3,
    /// UV on the first face passed by the caller.
    pub uv_a: (f64, f64),
    /// UV on the second face passed by the caller.
    pub uv_b: (f64, f64),
}

/// A sample counts as Interior overlap when neither trim classifies
/// it `Outside`, it is classifiable on both sides, and it is not
/// exempt Boundary∧Boundary contact — the `_ => return true` arm of
/// the validator's classification.
fn interior_overlap_counts(
    trim_a: Option<TrimContainment>,
    trim_b: Option<TrimContainment>,
) -> bool {
    matches!(
        (trim_a, trim_b),
        (
            Some(TrimContainment::Interior),
            Some(TrimContainment::Interior | TrimContainment::Boundary),
        ) | (
            Some(TrimContainment::Boundary),
            Some(TrimContainment::Interior),
        )
    )
}

impl BRep {
    /// Interior-classified overlap samples between two faces — the
    /// same surface intersection and trim classification
    /// [`BRep::validate_self_intersection`] applies, returned with
    /// each sample's 3D point and its UV on each face (in the
    /// caller's argument order) so repair and diagnostic consumers
    /// can localize the overlap without re-projecting (iterative
    /// 3D→UV projection on free-form surfaces is exactly what the
    /// validator itself avoids). A sample counts when neither trim
    /// classifies it `Outside` and it is not exempt
    /// boundary/shared-vertex contact. Empty when the pair does not
    /// overlap (or one face is stale).
    pub fn interior_overlap_samples(
        &self,
        face_a: FaceHandle,
        face_b: FaceHandle,
    ) -> Vec<InteriorOverlapSample> {
        // Same canonical pair orientation as the self-intersection
        // check: the solver is order-sensitive for NURBS×analytic
        // pairs, and analytic-first is the direction that finds the
        // curve (PK-251).
        let is_free_form = |face: FaceHandle| {
            matches!(
                self.faces
                    .get(face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry),
                Some(SurfaceGeometry::Nurbs(_) | SurfaceGeometry::Gregory(_))
            )
        };
        let swapped = is_free_form(face_a) && !is_free_form(face_b);
        let (face_a, face_b) = if swapped {
            (face_b, face_a)
        } else {
            (face_a, face_b)
        };
        let Some(surface_a) = self
            .faces
            .get(face_a)
            .and_then(|face| face.surface.as_ref())
            .map(|sof| sof.geometry.clone())
        else {
            return Vec::new();
        };
        let Some(surface_b) = self
            .faces
            .get(face_b)
            .and_then(|face| face.surface.as_ref())
            .map(|sof| sof.geometry.clone())
        else {
            return Vec::new();
        };
        let linear_tol = self.tolerance.linear();
        let coarse_tol = crate::tolerance::Tolerance::default()
            .with_linear(linear_tol * 2.0)
            .with_angular(self.tolerance.angular())
            .with_parametric(self.tolerance.parametric_base());
        let positions_a = self.face_outer_loop_positions(face_a);
        let ((u0a, u1a), (v0a, v1a)) =
            crate::geometry::uv_projection::exact_bounded_parameter_range(&positions_a, &surface_a);
        let positions_b = self.face_outer_loop_positions(face_b);
        let ((u0b, u1b), (v0b, v1b)) =
            crate::geometry::uv_projection::exact_bounded_parameter_range(&positions_b, &surface_b);
        // Cylinder × free-form fast path (PK-251): implicit-distance
        // samples with exact free-form UVs — same reason and same
        // machinery as the validator's fast path (the marching solver
        // collapses on grazing fillet-seam contact).
        if let Some(samples) = crate::validation::self_intersection::analytic_free_form_near_samples(
            &surface_a,
            &surface_b,
            (u0b, u1b),
            (v0b, v1b),
            linear_tol,
        ) {
            let mut out = Vec::new();
            for sample in samples {
                let trim_a =
                    self.face_trim_classify_point(face_a, &surface_a, &sample.point, linear_tol);
                let trim_b = self.face_trim_classify_uv(
                    face_b,
                    &surface_b,
                    UvPoint {
                        u: sample.free_form_uv.0,
                        v: sample.free_form_uv.1,
                    },
                );
                if interior_overlap_counts(trim_a, trim_b) {
                    let uv_a = crate::geometry::uv_projection::project_point_to_uv(
                        &surface_a,
                        &sample.point,
                    )
                    .unwrap_or((f64::NAN, f64::NAN));
                    out.push(InteriorOverlapSample {
                        point: sample.point,
                        uv_a,
                        uv_b: sample.free_form_uv,
                    });
                }
            }
            if swapped {
                for sample in &mut out {
                    std::mem::swap(&mut sample.uv_a, &mut sample.uv_b);
                }
            }
            return out;
        }
        let result = crate::geometry::intersection::intersect_surfaces_with_handler(
            &surface_a,
            (u0a, u1a),
            (v0a, v1a),
            &surface_b,
            (u0b, u1b),
            (v0b, v1b),
            &coarse_tol,
            &crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler,
        );
        let Ok(result) = result else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for intersection in &result.intersections {
            // Mirror `intersection_overlaps_face_trims`: free-form
            // sides classify through the intersection's UV curves
            // (iterative 3D→UV projection is unreliable on NURBS —
            // the same reason the validator branches here).
            let use_uv_curves = matches!(intersection, SurfaceSurfaceIntersection::Curve(_))
                && (surface_uses_iterative_projection(&surface_a)
                    || surface_uses_iterative_projection(&surface_b));
            if use_uv_curves {
                let SurfaceSurfaceIntersection::Curve(curve) = intersection else {
                    continue;
                };
                let (t0, t1) = curve.curve_3d.parameter_range();
                let (uv1_t0, uv1_t1) = curve.curve_uv1.parameter_range();
                let (uv2_t0, uv2_t1) = curve.curve_uv2.parameter_range();
                if !t0.is_finite()
                    || !t1.is_finite()
                    || !uv1_t0.is_finite()
                    || !uv1_t1.is_finite()
                    || !uv2_t0.is_finite()
                    || !uv2_t1.is_finite()
                    || (t1 - t0).abs() < f64::EPSILON
                {
                    continue;
                }
                const SAMPLE_COUNT: usize = 33;
                for i in 0..SAMPLE_COUNT {
                    let alpha = i as f64 / (SAMPLE_COUNT - 1) as f64;
                    let uv_a = curve.curve_uv1.point_at(uv1_t0 + (uv1_t1 - uv1_t0) * alpha);
                    let uv_b = curve.curve_uv2.point_at(uv2_t0 + (uv2_t1 - uv2_t0) * alpha);
                    let trim_a = self.face_trim_classify_uv(
                        face_a,
                        &surface_a,
                        UvPoint {
                            u: uv_a.x,
                            v: uv_a.y,
                        },
                    );
                    let trim_b = self.face_trim_classify_uv(
                        face_b,
                        &surface_b,
                        UvPoint {
                            u: uv_b.x,
                            v: uv_b.y,
                        },
                    );
                    if interior_overlap_counts(trim_a, trim_b) {
                        out.push(InteriorOverlapSample {
                            point: curve.curve_3d.point_at(t0 + (t1 - t0) * alpha),
                            uv_a: (uv_a.x, uv_a.y),
                            uv_b: (uv_b.x, uv_b.y),
                        });
                    }
                }
                continue;
            }
            for point in intersection_sample_points(&surface_a, intersection) {
                let trim_a = self.face_trim_classify_point(face_a, &surface_a, &point, linear_tol);
                let trim_b = self.face_trim_classify_point(face_b, &surface_b, &point, linear_tol);
                if interior_overlap_counts(trim_a, trim_b) {
                    let uv_a =
                        crate::geometry::uv_projection::project_point_to_uv(&surface_a, &point)
                            .unwrap_or((f64::NAN, f64::NAN));
                    let uv_b =
                        crate::geometry::uv_projection::project_point_to_uv(&surface_b, &point)
                            .unwrap_or((f64::NAN, f64::NAN));
                    out.push(InteriorOverlapSample { point, uv_a, uv_b });
                }
            }
        }
        if swapped {
            for sample in &mut out {
                std::mem::swap(&mut sample.uv_a, &mut sample.uv_b);
            }
        }
        out
    }
}
