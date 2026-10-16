//! Pair-class fast paths for the self-intersection check: bounded
//! implicit-distance sampling verdicts for cylinder×cylinder and
//! analytic×free-form pairs, where the generic marching solver is
//! disproportionate or pathological (PK-251).

use super::self_intersection_helpers::{
    analytic_free_form_near_samples, cylinder_cylinder_near_samples,
};
use crate::topology::BRep;

impl BRep {
    /// PK-251 fast path: cylinder × free-form (NURBS/Gregory) pair
    /// verdict via implicit-distance samples of the crossing curve,
    /// classified with the same two-tier trim semantics as
    /// `sampled_points_overlap_faces` — but the free-form side
    /// classifies by its EXACT sample UV rather than iterative
    /// projection. `None` when the pair is not this class.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn bounded_analytic_free_form_trim_overlap(
        &self,
        cylinder_face: crate::topology::FaceHandle,
        cylinder_surface: &crate::geometry::surface::SurfaceGeometry,
        free_form_face: crate::topology::FaceHandle,
        free_form_surface: &crate::geometry::surface::SurfaceGeometry,
        free_form_u_range: (f64, f64),
        free_form_v_range: (f64, f64),
        linear_tol: f64,
    ) -> Option<bool> {
        let samples = analytic_free_form_near_samples(
            cylinder_surface,
            free_form_surface,
            free_form_u_range,
            free_form_v_range,
            linear_tol,
        )?;
        Some(self.analytic_free_form_samples_overlap(
            cylinder_face,
            cylinder_surface,
            free_form_face,
            free_form_surface,
            &samples,
            linear_tol,
        ))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn bounded_cylinder_cylinder_trim_overlap(
        &self,
        face_a: crate::topology::FaceHandle,
        surface_a: &crate::geometry::surface::SurfaceGeometry,
        u_range_a: (f64, f64),
        v_range_a: (f64, f64),
        face_b: crate::topology::FaceHandle,
        surface_b: &crate::geometry::surface::SurfaceGeometry,
        u_range_b: (f64, f64),
        v_range_b: (f64, f64),
        linear_tol: f64,
    ) -> Option<bool> {
        let samples = cylinder_cylinder_near_samples(
            surface_a, u_range_a, v_range_a, surface_b, u_range_b, v_range_b, linear_tol,
        )?;
        // The grid keeps points within 10 tolerances of the other
        // cylinder: near the crossing, not on it. Classified raw, a
        // point inside one strip and in the other's boundary band read
        // as an overlap though it sits up to 10 tolerances off the
        // other surface — a full radius on a corner whose r is 10
        // tolerances, where two strips meet only at the ball's vertex
        // (PK-866). Each sample moves onto both surfaces first; one
        // with no crossing within the tolerance is a near miss.
        let crossings: Vec<_> = samples
            .iter()
            .filter_map(|p| cylinder_pair_crossing_near(surface_a, surface_b, *p, linear_tol))
            .collect();
        if crossings.is_empty() {
            return Some(false);
        }
        Some(self.sampled_points_overlap_faces(
            face_a, surface_a, face_b, surface_b, &crossings, linear_tol,
        ))
    }
}

/// The point of both cylinders nearest `start`, by alternating exact
/// projections, or `None` when the sheets do not meet within
/// `linear_tol` there (a near miss; tangent sheets converge slowly but
/// within the iteration cap at tolerance scale).
pub(super) fn cylinder_pair_crossing_near(
    surface_a: &crate::geometry::surface::SurfaceGeometry,
    surface_b: &crate::geometry::surface::SurfaceGeometry,
    start: crate::geometry::Point3,
    linear_tol: f64,
) -> Option<crate::geometry::Point3> {
    use crate::geometry::surface::SurfaceGeometry;
    let project = |surface: &SurfaceGeometry, p: crate::geometry::Point3| match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let offset = p - *origin;
            let along = *origin + *axis * offset.dot(axis);
            let radial = p - along;
            let length = radial.length();
            (length > 0.0).then(|| along + radial * (*radius / length))
        }
        _ => None,
    };
    let mut p = start;
    for _ in 0..64 {
        let on_a = project(surface_a, p)?;
        let on_b = project(surface_b, on_a)?;
        let gap = on_a.distance_to(&on_b);
        p = on_b;
        if gap <= 0.1 * linear_tol {
            return Some(p);
        }
    }
    let gap = project(surface_a, p)?.distance_to(&p);
    (gap <= linear_tol).then_some(p)
}
