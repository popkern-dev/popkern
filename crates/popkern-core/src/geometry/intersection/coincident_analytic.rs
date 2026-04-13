//! Analytic coincidence detection for same-type surface pairs.
//!
//! Provides an O(1) check that two analytic surfaces of the same type represent
//! the same geometric locus (within tolerance). Returns `false` for NURBS and
//! mixed-type pairs — callers should fall back to sample-based detection for
//! those cases.

use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Returns `true` if two analytic surfaces of the same type represent the same
/// geometric locus within the given tolerance.
///
/// This is a fast O(1) comparison that handles Plane, Cylinder, Cone, Sphere,
/// and Torus pairs. For NURBS surfaces or mixed-type pairs it always returns
/// `false` — the caller should fall back to sample-based coincidence detection.
pub fn are_surfaces_coincident_analytic(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> bool {
    match (surface1, surface2) {
        // Plane-Plane: normals parallel (or anti-parallel) AND origins on the same plane.
        (
            SurfaceGeometry::Plane {
                origin: o1,
                normal: n1,
                ..
            },
            SurfaceGeometry::Plane {
                origin: o2,
                normal: n2,
                ..
            },
        ) => {
            let normals_parallel = n1.dot(n2).abs() > 1.0 - tolerance.angular();
            let d = *o1 - *o2;
            let origin_on_plane = d.dot(n1).abs() < tolerance.linear();
            normals_parallel && origin_on_plane
        }

        // Cylinder-Cylinder: axes parallel, same radius, colinear axes.
        (
            SurfaceGeometry::Cylinder {
                origin: o1,
                axis: a1,
                radius: r1,
            },
            SurfaceGeometry::Cylinder {
                origin: o2,
                axis: a2,
                radius: r2,
            },
        ) => {
            let axes_parallel = a1.dot(a2).abs() > 1.0 - tolerance.angular();
            let same_radius = (r1 - r2).abs() < tolerance.linear();
            // Colinear axes: perpendicular component of origin separation < tolerance
            let d = *o1 - *o2;
            let along = d.dot(a1);
            let perp = d - *a1 * along;
            // `perp.length() < linear` rewritten as adaptive squared-
            // distance check from the origin to the perp tip.
            let axes_colinear = crate::adaptive_precision::adaptive_dist_sq_lt(
                (0.0, 0.0, 0.0),
                (perp.x, perp.y, perp.z),
                tolerance.linear() * tolerance.linear(),
            );
            axes_parallel && same_radius && axes_colinear
        }

        // Cone-Cone: same apex, parallel axes, same half_angle.
        (
            SurfaceGeometry::Cone {
                apex: apex1,
                axis: a1,
                half_angle: ha1,
            },
            SurfaceGeometry::Cone {
                apex: apex2,
                axis: a2,
                half_angle: ha2,
            },
        ) => {
            let same_apex = apex1.distance_to(apex2) < tolerance.linear();
            let axes_parallel = a1.dot(a2).abs() > 1.0 - tolerance.angular();
            let same_half_angle = (ha1 - ha2).abs() < tolerance.angular();
            same_apex && axes_parallel && same_half_angle
        }

        // Sphere-Sphere: same center, same radius.
        (
            SurfaceGeometry::Sphere {
                center: c1,
                radius: r1,
            },
            SurfaceGeometry::Sphere {
                center: c2,
                radius: r2,
            },
        ) => {
            let same_center = c1.distance_to(c2) < tolerance.linear();
            let same_radius = (r1 - r2).abs() < tolerance.linear();
            same_center && same_radius
        }

        // Torus-Torus: same center, parallel axes, same major_radius, same minor_radius.
        (
            SurfaceGeometry::Torus {
                center: c1,
                axis: a1,
                major_radius: maj1,
                minor_radius: min1,
            },
            SurfaceGeometry::Torus {
                center: c2,
                axis: a2,
                major_radius: maj2,
                minor_radius: min2,
            },
        ) => {
            let same_center = c1.distance_to(c2) < tolerance.linear();
            let axes_parallel = a1.dot(a2).abs() > 1.0 - tolerance.angular();
            let same_major = (maj1 - maj2).abs() < tolerance.linear();
            let same_minor = (min1 - min2).abs() < tolerance.linear();
            same_center && axes_parallel && same_major && same_minor
        }

        // NURBS or mixed types: not handled analytically.
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Point3, Vec3};

    fn default_tolerance() -> Tolerance {
        Tolerance::default()
    }

    // -----------------------------------------------------------------------
    // Plane tests
    // -----------------------------------------------------------------------

    #[test]
    fn identical_planes_are_coincident() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        assert!(are_surfaces_coincident_analytic(
            &plane,
            &plane,
            &default_tolerance()
        ));
    }

    #[test]
    fn same_plane_different_origin_on_plane() {
        let p1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let p2 = SurfaceGeometry::Plane {
            origin: Point3::new(99.0, -42.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        assert!(are_surfaces_coincident_analytic(
            &p1,
            &p2,
            &default_tolerance()
        ));
    }

    #[test]
    fn same_plane_flipped_normal() {
        let p1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let p2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, -1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, -1.0, 0.0),
        };
        assert!(are_surfaces_coincident_analytic(
            &p1,
            &p2,
            &default_tolerance()
        ));
    }

    #[test]
    fn parallel_offset_planes_not_coincident() {
        let p1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let p2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        assert!(!are_surfaces_coincident_analytic(
            &p1,
            &p2,
            &default_tolerance()
        ));
    }

    #[test]
    fn perpendicular_planes_not_coincident() {
        let p1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let p2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        assert!(!are_surfaces_coincident_analytic(
            &p1,
            &p2,
            &default_tolerance()
        ));
    }

    // -----------------------------------------------------------------------
    // Cylinder tests
    // -----------------------------------------------------------------------

    #[test]
    fn identical_cylinders_are_coincident() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        assert!(are_surfaces_coincident_analytic(
            &cyl,
            &cyl,
            &default_tolerance()
        ));
    }

    #[test]
    fn same_cylinder_different_origin_on_axis() {
        let c1 = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        let c2 = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 5.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        assert!(are_surfaces_coincident_analytic(
            &c1,
            &c2,
            &default_tolerance()
        ));
    }

    #[test]
    fn different_radius_cylinders_not_coincident() {
        let c1 = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        let c2 = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
        };
        assert!(!are_surfaces_coincident_analytic(
            &c1,
            &c2,
            &default_tolerance()
        ));
    }

    #[test]
    fn offset_axis_cylinders_not_coincident() {
        let c1 = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        let c2 = SurfaceGeometry::Cylinder {
            origin: Point3::new(1.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 5.0,
        };
        assert!(!are_surfaces_coincident_analytic(
            &c1,
            &c2,
            &default_tolerance()
        ));
    }

    // -----------------------------------------------------------------------
    // Sphere tests
    // -----------------------------------------------------------------------

    #[test]
    fn identical_spheres_are_coincident() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(1.0, 2.0, 3.0),
            radius: 4.0,
        };
        assert!(are_surfaces_coincident_analytic(
            &sphere,
            &sphere,
            &default_tolerance()
        ));
    }

    #[test]
    fn different_center_spheres_not_coincident() {
        let s1 = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let s2 = SurfaceGeometry::Sphere {
            center: Point3::new(1.0, 0.0, 0.0),
            radius: 1.0,
        };
        assert!(!are_surfaces_coincident_analytic(
            &s1,
            &s2,
            &default_tolerance()
        ));
    }

    // -----------------------------------------------------------------------
    // Cone tests
    // -----------------------------------------------------------------------

    #[test]
    fn identical_cones_are_coincident() {
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 1.0),
            axis: Vec3::new(0.0, 0.0, -1.0),
            half_angle: std::f64::consts::FRAC_PI_4,
        };
        assert!(are_surfaces_coincident_analytic(
            &cone,
            &cone,
            &default_tolerance()
        ));
    }

    // -----------------------------------------------------------------------
    // Torus tests
    // -----------------------------------------------------------------------

    #[test]
    fn identical_tori_are_coincident() {
        let torus = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 2.0,
            minor_radius: 0.5,
        };
        assert!(are_surfaces_coincident_analytic(
            &torus,
            &torus,
            &default_tolerance()
        ));
    }

    // -----------------------------------------------------------------------
    // Mixed type / NURBS tests
    // -----------------------------------------------------------------------

    #[test]
    fn plane_vs_cylinder_not_coincident() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        assert!(!are_surfaces_coincident_analytic(
            &plane,
            &cyl,
            &default_tolerance()
        ));
    }

    #[test]
    fn nurbs_returns_false() {
        let nurbs = SurfaceGeometry::Nurbs(crate::geometry::nurbs_surface::NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0, 1.0, 1.0, 1.0],
        });
        assert!(!are_surfaces_coincident_analytic(
            &nurbs,
            &nurbs,
            &default_tolerance()
        ));
    }
}
