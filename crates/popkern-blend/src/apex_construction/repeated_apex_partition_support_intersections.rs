use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use super::partition_cylinder_geometry::point_lies_on_cylinder_face;

pub(super) fn support_plane_faces_containing_point(
    brep: &BRep,
    cylinder_face: FaceHandle,
    point: Point3,
) -> Vec<FaceHandle> {
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            ((point - *origin).dot(normal).abs() <= tol).then_some(face)
        })
        .collect()
}

pub(super) fn line_segment_lies_on_cylinder_and_plane(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: Point3,
    end: Point3,
    plane_origin: Point3,
    plane_normal: Vec3,
) -> bool {
    let direction = end - start;
    if direction.length() <= brep.tolerance.linear().max(1.0e-12) {
        return false;
    }
    let tol = brep.tolerance.linear().max(1.0e-9);
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let point = start + direction * t;
        if (point - plane_origin).dot(&plane_normal).abs() > tol {
            return false;
        }
        if !point_lies_on_cylinder_face(brep, cylinder_face, point) {
            return false;
        }
    }
    true
}

pub(super) fn plane_plane_cylinder_joint_points(
    brep: &BRep,
    cylinder_face: FaceHandle,
    plane_a: FaceHandle,
    plane_b: FaceHandle,
) -> Vec<Point3> {
    let Some(SurfaceGeometry::Cylinder {
        origin: cylinder_origin,
        axis,
        radius,
    }) = brep
        .faces
        .get(cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let Some(SurfaceGeometry::Plane {
        origin: plane_a_origin,
        normal: plane_a_normal,
        ..
    }) = brep
        .faces
        .get(plane_a)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let Some(SurfaceGeometry::Plane {
        origin: plane_b_origin,
        normal: plane_b_normal,
        ..
    }) = brep
        .faces
        .get(plane_b)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let n1 = *plane_a_normal;
    let n2 = *plane_b_normal;
    let direction = n1.cross(&n2);
    let denom = direction.dot(&direction);
    if denom <= 1.0e-24 {
        return Vec::new();
    }
    let c1 = n1.dot(&(*plane_a_origin - Point3::new(0.0, 0.0, 0.0)));
    let c2 = n2.dot(&(*plane_b_origin - Point3::new(0.0, 0.0, 0.0)));
    let point_vec = (n2.cross(&direction) * c1 + direction.cross(&n1) * c2) * (1.0 / denom);
    let line_point = Point3::new(point_vec.x, point_vec.y, point_vec.z);
    line_cylinder_intersections(line_point, direction, *cylinder_origin, *axis, *radius)
        .into_iter()
        .filter(|point| point_lies_on_cylinder_face(brep, cylinder_face, *point))
        .filter(|point| {
            let tol = brep.tolerance.linear().max(1.0e-9);
            ((*point - *plane_a_origin).dot(plane_a_normal)).abs() <= tol
                && ((*point - *plane_b_origin).dot(plane_b_normal)).abs() <= tol
        })
        .collect()
}

fn line_cylinder_intersections(
    line_point: Point3,
    line_direction: Vec3,
    cylinder_origin: Point3,
    cylinder_axis: Vec3,
    radius: f64,
) -> Vec<Point3> {
    let Some(axis_unit) = cylinder_axis.normalized() else {
        return Vec::new();
    };
    let rel = line_point - cylinder_origin;
    let direction_perp = line_direction - axis_unit * line_direction.dot(&axis_unit);
    let rel_perp = rel - axis_unit * rel.dot(&axis_unit);
    let a = direction_perp.dot(&direction_perp);
    let b = 2.0 * rel_perp.dot(&direction_perp);
    let c = rel_perp.dot(&rel_perp) - radius * radius;
    if a.abs() <= 1.0e-24 {
        return Vec::new();
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < -1.0e-12 {
        return Vec::new();
    }
    if discriminant.abs() <= 1.0e-12 {
        let t = -b / (2.0 * a);
        return vec![line_point + line_direction * t];
    }
    let root = discriminant.sqrt();
    let t0 = (-b - root) / (2.0 * a);
    let t1 = (-b + root) / (2.0 * a);
    vec![
        line_point + line_direction * t0,
        line_point + line_direction * t1,
    ]
}
