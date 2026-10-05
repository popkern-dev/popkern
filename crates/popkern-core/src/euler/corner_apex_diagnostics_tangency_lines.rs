//! Line intersections used as exact seeds for three-face tangency diagnostics.

use super::solve_3x3;
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

#[derive(Clone, Copy)]
pub(super) struct PlaneEquation {
    pub(super) normal: Vec3,
    pub(super) offset: f64,
}

#[derive(Clone, Copy)]
pub(super) struct LineEquation {
    point: Point3,
    direction: Vec3,
}

#[derive(Clone, Copy)]
pub(super) struct CylinderEquation {
    pub(super) origin: Point3,
    pub(super) axis: Vec3,
    pub(super) radius: f64,
}

#[derive(Clone, Copy)]
pub(super) struct ConeEquation {
    pub(super) apex: Point3,
    pub(super) axis: Vec3,
    pub(super) half_angle: f64,
}

pub(super) fn plane_plane_offset_line(
    first: PlaneEquation,
    first_offset_delta: f64,
    second: PlaneEquation,
    second_offset_delta: f64,
    near: Point3,
) -> Option<LineEquation> {
    let direction = first.normal.cross(&second.normal).normalized()?;
    let near = Vec3::new(near.x, near.y, near.z);
    let point = solve_3x3(
        [first.normal, second.normal, direction],
        Vec3::new(
            first.offset + first_offset_delta,
            second.offset + second_offset_delta,
            direction.dot(&near),
        ),
    )?;
    Some(LineEquation {
        point: Point3::new(point.x, point.y, point.z),
        direction,
    })
}

pub(super) fn line_sphere_intersections(
    line: LineEquation,
    center: Point3,
    radius: f64,
) -> Vec<Point3> {
    let to_point = line.point - center;
    roots(
        1.0,
        2.0 * to_point.dot(&line.direction),
        to_point.dot(&to_point) - radius * radius,
    )
    .into_iter()
    .map(move |t| line.point + line.direction * t)
    .collect()
}

pub(super) fn line_cylinder_intersections(
    line: LineEquation,
    cylinder: CylinderEquation,
    radius: f64,
) -> Vec<Point3> {
    let from_axis = line.point - cylinder.origin;
    let direction = line.direction - cylinder.axis * line.direction.dot(&cylinder.axis);
    let offset = from_axis - cylinder.axis * from_axis.dot(&cylinder.axis);
    let a = direction.dot(&direction);
    let b = 2.0 * offset.dot(&direction);
    let c = offset.dot(&offset) - radius * radius;
    roots(a, b, c)
        .into_iter()
        .map(move |t| line.point + line.direction * t)
        .collect()
}

pub(super) fn line_cone_intersections(
    line: LineEquation,
    cone: ConeEquation,
    offset_delta: f64,
) -> Vec<Point3> {
    let from_apex = line.point - cone.apex;
    let h0 = from_apex.dot(&cone.axis);
    let hd = line.direction.dot(&cone.axis);
    let radial0 = from_apex - cone.axis * h0;
    let radial_direction = line.direction - cone.axis * hd;
    let sin = cone.half_angle.psin();
    let cos2 = cone.half_angle.pcos().powi(2);
    let base = offset_delta + sin * h0;
    let a = cos2 * radial_direction.dot(&radial_direction) - (sin * hd).powi(2);
    let b = 2.0 * (cos2 * radial0.dot(&radial_direction) - base * sin * hd);
    let c = cos2 * radial0.dot(&radial0) - base * base;
    roots(a, b, c)
        .into_iter()
        .filter(move |t| offset_delta + sin * (h0 + hd * *t) >= -1e-9)
        .map(move |t| line.point + line.direction * t)
        .collect()
}

fn roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() <= 1e-24 {
        return (b.abs() > 1e-24).then(|| -c / b).into_iter().collect();
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < -1e-12 {
        Vec::new()
    } else if discriminant.abs() <= 1e-12 {
        vec![-b / (2.0 * a)]
    } else {
        let root = discriminant.sqrt();
        vec![(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]
    }
}
