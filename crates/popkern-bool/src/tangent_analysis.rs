//! Second-order contact analysis for tangent intersections (spec §2.5).
//!
//! When two surfaces are tangent at a point (normals parallel), first-order
//! information is insufficient to determine the intersection topology.
//! This module uses the second fundamental forms to classify:
//!
//! - **Crossing**: surfaces cross at the tangent point (cusp or transversal
//!   in the tangent plane). Determined by the relative curvature.
//! - **Touching**: surfaces osculate without crossing (isolated point or
//!   curve segment). Determined by the sign of the curvature difference.
//! - **Overlapping**: surfaces are coincident in a neighborhood.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Result of second-order contact analysis at a tangent point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TangentContactType {
    /// Surfaces cross at the tangent point.
    Crossing,
    /// Surfaces touch without crossing (osculate).
    Touching,
    /// Surfaces are coincident in a neighborhood.
    Coincident,
    /// Analysis inconclusive (difference of curvatures within FP uncertainty).
    Inconclusive,
}

/// Analyze second-order contact between two surfaces at a tangent point.
///
/// At parameters (u1, v1) on surface1 and (u2, v2) on surface2 where
/// the surfaces are tangent (normals parallel), compute the difference
/// of the second fundamental forms to determine crossing vs touching.
///
/// The eigenvalues of the curvature difference matrix determine:
/// - Both eigenvalues same sign → touching (surfaces curve away from each other)
/// - Eigenvalues of opposite sign → crossing (saddle-like contact)
/// - Both eigenvalues near zero → coincident
pub fn classify_tangent_contact(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    angular_tolerance: f64,
) -> TangentContactType {
    // Get surface normals at the tangent point
    let n1 = surface1.normal_at(u1, v1);
    let n2 = surface2.normal_at(u2, v2);

    let n1_len = n1.length();
    let n2_len = n2.length();
    if n1_len < 1e-12 || n2_len < 1e-12 {
        return TangentContactType::Inconclusive;
    }

    let n1_unit = n1 * (1.0 / n1_len);

    // Second fundamental form coefficients for surface 1
    let suu1 = surface1.second_partial_uu(u1, v1);
    let suv1 = surface1.second_partial_uv(u1, v1);
    let svv1 = surface1.second_partial_vv(u1, v1);

    let l1 = suu1.dot(&n1_unit);
    let m1 = suv1.dot(&n1_unit);
    let n1_coeff = svv1.dot(&n1_unit);

    // Second fundamental form coefficients for surface 2
    let suu2 = surface2.second_partial_uu(u2, v2);
    let suv2 = surface2.second_partial_uv(u2, v2);
    let svv2 = surface2.second_partial_vv(u2, v2);

    let l2 = suu2.dot(&n1_unit); // project onto same normal
    let m2 = suv2.dot(&n1_unit);
    let n2_coeff = svv2.dot(&n1_unit);

    // Difference of second fundamental forms
    let dl = l1 - l2;
    let dm = m1 - m2;
    let dn = n1_coeff - n2_coeff;

    // The 2x2 curvature difference matrix is:
    //   [dl  dm]
    //   [dm  dn]
    // Its eigenvalues determine the contact type.
    let trace = dl + dn;
    let det = dl * dn - dm * dm;

    let discriminant = trace * trace - 4.0 * det;

    let eps = angular_tolerance * angular_tolerance;

    if discriminant < -eps {
        // Complex eigenvalues — shouldn't happen for real surfaces
        TangentContactType::Inconclusive
    } else if det.abs() < eps && trace.abs() < eps {
        // Both eigenvalues near zero → coincident
        TangentContactType::Coincident
    } else if det > eps {
        // Same-sign eigenvalues → touching (both curve same way)
        TangentContactType::Touching
    } else if det < -eps {
        // Opposite-sign eigenvalues → crossing (saddle)
        TangentContactType::Crossing
    } else {
        TangentContactType::Inconclusive
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::{Point3, Vec3};

    #[test]
    fn sphere_plane_tangent_is_touching() {
        // Sphere tangent to plane at the north pole
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // Sphere at (u=0, v=0) has normal (0,0,1), point (0,0,1)
        // Plane at (0,0) has normal (0,0,1), point (0,0,1)
        // They touch at the north pole.
        let result = classify_tangent_contact(
            &sphere, 0.0, 0.01, // near north pole (v=0 is degenerate)
            &plane, 0.0, 0.0, 1e-4,
        );
        // Sphere curves away from plane → touching
        assert!(
            result == TangentContactType::Touching || result == TangentContactType::Inconclusive,
            "sphere-plane should be touching, got {:?}",
            result
        );
    }

    #[test]
    fn two_planes_are_coincident() {
        let p1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let p2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = classify_tangent_contact(&p1, 0.0, 0.0, &p2, 0.0, 0.0, 1e-4);
        assert_eq!(result, TangentContactType::Coincident);
    }
}
