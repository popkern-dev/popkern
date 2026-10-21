//! Adaptive precision for geometric predicates.
//!
//! The foundational types ([`DoublePrecision`] and [`adaptive_sign`]) live in
//! `popkern_core::adaptive_precision` and are re-exported here for
//! backwards compatibility. This module adds `robust_point_in_solid`,
//! which depends on bool-specific types, and `body_face_aabbs`, the
//! per-face boxes every point-in-solid query takes.

pub use popkern_core::adaptive_precision::{DoublePrecision, adaptive_sign};

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::classify::point_in_solid;
use crate::classify_winding::point_in_solid_winding;
use crate::error::PointClassification;

/// The per-face boxes the point-in-solid queries take
/// ([`robust_point_in_solid`], [`crate::point_in_solid`] and their
/// `_final` forms): every face with surface geometry, boxed as the
/// boolean's intersection stage boxes it (trim-scale where the surface
/// allows). Callers outside the boolean had no way to build them.
pub fn body_face_aabbs(brep: &BRep) -> Vec<(FaceHandle, Aabb)> {
    crate::intersect::compute_face_aabbs(brep)
}

/// Robust point-in-solid classification with adaptive fallback.
///
/// First tries multi-ray classification (fast). If the result is
/// ambiguous (OnBoundary), falls back to winding number (more robust
/// but slower). This provides the best of both methods.
pub fn robust_point_in_solid(
    point: &Point3,
    solid: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> PointClassification {
    // Fast path: multi-ray classification
    match point_in_solid(point, solid, face_aabbs, tolerance) {
        Ok(PointClassification::Inside) => PointClassification::Inside,
        Ok(PointClassification::Outside) => PointClassification::Outside,
        Ok(PointClassification::OnBoundary) | Err(_) => {
            // Ambiguous or failed — fall back to winding number
            point_in_solid_winding(point, solid)
        }
    }
}

/// [`robust_point_in_solid_final`] against a prepared solid (face
/// windows kept across classifications).
pub(crate) fn robust_point_in_solid_final_on(
    point: &Point3,
    rays: &crate::classify::RaySolid<'_>,
    tolerance: &Tolerance,
) -> PointClassification {
    match crate::classify::point_in_solid_final_on(point, rays, tolerance) {
        Ok(PointClassification::Inside) => PointClassification::Inside,
        Ok(PointClassification::Outside) => PointClassification::Outside,
        Ok(PointClassification::OnBoundary) | Err(_) => point_in_solid_winding(point, rays.solid),
    }
}

/// [`robust_point_in_solid`] against a FINAL body: hits through the
/// holes of holed faces are genuine voids and don't count (see
/// `point_in_solid_final`).
pub fn robust_point_in_solid_final(
    point: &Point3,
    solid: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> PointClassification {
    match crate::classify::point_in_solid_final(point, solid, face_aabbs, tolerance) {
        Ok(PointClassification::Inside) => PointClassification::Inside,
        Ok(PointClassification::Outside) => PointClassification::Outside,
        Ok(PointClassification::OnBoundary) | Err(_) => {
            // Ambiguous or failed — fall back to winding number
            point_in_solid_winding(point, solid)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn robust_classification_inside_cube() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let aabbs = crate::intersect::compute_face_aabbs(&cube);
        let pt = Point3::new(0.5, 0.5, 0.5);
        let result = robust_point_in_solid(&pt, &cube, &aabbs, &Tolerance::default());
        assert_eq!(result, PointClassification::Inside);
    }

    #[test]
    fn robust_classification_outside_cube() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let aabbs = crate::intersect::compute_face_aabbs(&cube);
        let pt = Point3::new(5.0, 5.0, 5.0);
        let result = robust_point_in_solid(&pt, &cube, &aabbs, &Tolerance::default());
        assert_eq!(result, PointClassification::Outside);
    }

    #[test]
    fn adaptive_sign_clear() {
        assert_eq!(adaptive_sign(1.0, 1e-10), Some(std::cmp::Ordering::Greater));
        assert_eq!(adaptive_sign(-1.0, 1e-10), Some(std::cmp::Ordering::Less));
    }

    #[test]
    fn adaptive_sign_ambiguous() {
        assert_eq!(adaptive_sign(1e-15, 1e-10), None);
    }

    #[test]
    fn double_precision_add() {
        let a = DoublePrecision::from_f64(1.0);
        let b = DoublePrecision::from_f64(1e-16);
        let sum = a.dd_add(b);
        // Should preserve the small value that f64 would lose
        assert!((sum.to_f64() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn double_precision_mul() {
        let a = DoublePrecision::from_f64(3.0);
        let b = DoublePrecision::from_f64(7.0);
        let product = a.dd_mul(b);
        assert!((product.to_f64() - 21.0).abs() < 1e-14);
    }
}
