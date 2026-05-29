//! Surface-family equivalence predicates for exact-boundary commit support lookup.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

pub(super) fn equivalent_surface(a: &SurfaceGeometry, b: &SurfaceGeometry, tolerance: f64) -> bool {
    match (a, b) {
        (
            SurfaceGeometry::Plane {
                origin: origin_a,
                normal: normal_a,
                ..
            },
            SurfaceGeometry::Plane {
                origin: origin_b,
                normal: normal_b,
                ..
            },
        ) => equivalent_planes(*origin_a, *normal_a, *origin_b, *normal_b, tolerance),
        (
            SurfaceGeometry::Cylinder {
                origin: origin_a,
                axis: axis_a,
                radius: radius_a,
            },
            SurfaceGeometry::Cylinder {
                origin: origin_b,
                axis: axis_b,
                radius: radius_b,
            },
        ) => equivalent_cylinders(
            *origin_a, *axis_a, *radius_a, *origin_b, *axis_b, *radius_b, tolerance,
        ),
        _ => false,
    }
}

fn equivalent_planes(
    origin_a: Point3,
    normal_a: Vec3,
    origin_b: Point3,
    normal_b: Vec3,
    tolerance: f64,
) -> bool {
    let len_a = normal_a.length();
    let len_b = normal_b.length();
    if len_a <= 1.0e-12 || len_b <= 1.0e-12 {
        return false;
    }
    let parallel = (normal_a.dot(&normal_b) / (len_a * len_b)).abs() >= 1.0 - 1.0e-10;
    let offset = ((origin_b - origin_a).dot(&normal_a) / len_a).abs();
    parallel && offset <= tolerance.max(1.0e-9)
}

fn equivalent_cylinders(
    origin_a: Point3,
    axis_a: Vec3,
    radius_a: f64,
    origin_b: Point3,
    axis_b: Vec3,
    radius_b: f64,
    tolerance: f64,
) -> bool {
    let len_a = axis_a.length();
    let len_b = axis_b.length();
    if len_a <= 1.0e-12 || len_b <= 1.0e-12 {
        return false;
    }
    let parallel = (axis_a.dot(&axis_b) / (len_a * len_b)).abs() >= 1.0 - 1.0e-10;
    let axis_gap = (origin_b - origin_a).cross(&axis_a).length() / len_a;
    parallel && axis_gap <= tolerance.max(1.0e-9) && (radius_a - radius_b).abs() <= tolerance
}
