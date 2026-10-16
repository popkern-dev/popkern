//! Surface-type to display-color mapping for the debug viewer.

use popkern_core::geometry::surface::SurfaceGeometry;

/// Map a surface geometry type to an RGBA display color.
///
/// Each analytic surface kind gets a distinct color so the user can
/// visually distinguish planes, cylinders, spheres, etc.
pub fn surface_color(surface: &SurfaceGeometry) -> [f32; 4] {
    match surface {
        SurfaceGeometry::Plane { .. } => [0.7, 0.7, 0.7, 1.0], // gray
        SurfaceGeometry::Cylinder { .. } => [0.3, 0.5, 0.9, 1.0], // blue
        SurfaceGeometry::Sphere { .. } => [0.3, 0.8, 0.4, 1.0], // green
        SurfaceGeometry::Cone { .. } => [0.9, 0.6, 0.2, 1.0],  // orange
        SurfaceGeometry::Torus { .. } => [0.7, 0.3, 0.8, 1.0], // purple
        SurfaceGeometry::Nurbs(_) => [0.9, 0.9, 0.9, 1.0],     // white
        SurfaceGeometry::Gregory(_) => [0.9, 0.4, 0.6, 1.0], // magenta — distinguishes §5.3 corners
        SurfaceGeometry::PolyhedralSpline(_) => [0.6, 0.4, 0.9, 1.0], // violet — PnS cap sectors
    }
}

/// Color for selected entity highlight.
pub fn selection_color() -> [f32; 4] {
    [1.0, 0.9, 0.0, 1.0] // yellow
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn surface_color_plane_is_gray() {
        use popkern_core::geometry::{Point3, Vec3};
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let c = surface_color(&plane);
        assert!((c[0] - 0.7).abs() < f32::EPSILON);
        assert!((c[3] - 1.0).abs() < f32::EPSILON);
    }

    #[test]
    fn surface_color_cylinder_is_blue() {
        use popkern_core::geometry::{Point3, Vec3};
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let c = surface_color(&cyl);
        assert!((c[0] - 0.3).abs() < f32::EPSILON);
        assert!((c[1] - 0.5).abs() < f32::EPSILON);
        assert!((c[2] - 0.9).abs() < f32::EPSILON);
    }

    #[test]
    fn surface_color_sphere_is_green() {
        use popkern_core::geometry::Point3;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: popkern_core::geometry::Vec3::Z,
        };
        let c = surface_color(&sphere);
        assert!((c[1] - 0.8).abs() < f32::EPSILON);
    }

    #[test]
    fn surface_color_cone_is_orange() {
        use popkern_core::geometry::{Point3, Vec3};
        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 1.0),
            axis: Vec3::new(0.0, 0.0, -1.0),
            half_angle: 0.5,
        };
        let c = surface_color(&cone);
        assert!((c[0] - 0.9).abs() < f32::EPSILON);
        assert!((c[1] - 0.6).abs() < f32::EPSILON);
    }

    #[test]
    fn surface_color_torus_is_purple() {
        use popkern_core::geometry::{Point3, Vec3};
        let torus = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 2.0,
            minor_radius: 0.5,
        };
        let c = surface_color(&torus);
        assert!((c[0] - 0.7).abs() < f32::EPSILON);
        assert!((c[1] - 0.3).abs() < f32::EPSILON);
    }

    #[test]
    fn surface_color_nurbs_is_white() {
        use popkern_core::geometry::Point3;
        let nurbs = popkern_core::geometry::nurbs_surface::NurbsSurface {
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
        };
        let c = surface_color(&SurfaceGeometry::Nurbs(nurbs));
        assert!((c[0] - 0.9).abs() < f32::EPSILON);
    }

    #[test]
    fn selection_color_is_yellow() {
        let c = selection_color();
        assert!((c[0] - 1.0).abs() < f32::EPSILON);
        assert!((c[1] - 0.9).abs() < f32::EPSILON);
        assert!((c[2] - 0.0).abs() < f32::EPSILON);
        assert!((c[3] - 1.0).abs() < f32::EPSILON);
    }
}
