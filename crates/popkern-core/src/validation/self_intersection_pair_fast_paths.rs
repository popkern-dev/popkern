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
        if samples.is_empty() {
            return Some(false);
        }
        Some(self.sampled_points_overlap_faces(
            face_a, surface_a, face_b, surface_b, &samples, linear_tol,
        ))
    }
}
