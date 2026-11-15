//! The chamfer's foot on a curved face: the point of the face, square
//! to the edge, a chamfer distance from it in a straight line.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

/// The point of `surface` in the plane through `edge_point` square to
/// `tangent`, `distance` from `edge_point` in a straight line: Newton on
/// (|S − p|² − d², (S − p)·t) from the seed `(u, v)`. On a plane, or a
/// face whose section square to the edge is straight (a cylinder or a
/// cone along its ruling), it is the tangent step itself; on a sphere
/// or a torus the step projected back falls short of the chord (a
/// hemisphere's rim chamfered 1 on the ball reached radius 4.9029, not
/// 4.9). `None` where Newton does not converge within the seed's
/// neighbourhood, and the projected seed stands.
pub(in crate::chamfer_curve) fn chord_point_on_section(
    surface: &SurfaceGeometry,
    (mut u, mut v): (f64, f64),
    edge_point: &Point3,
    tangent: &Vec3,
    distance: f64,
) -> Option<Point3> {
    let scale = distance.max(f64::MIN_POSITIVE);
    for _ in 0..30 {
        let s = surface.point_at(u, v);
        let r = s - *edge_point;
        let f0 = r.dot(&r) - distance * distance;
        let f1 = r.dot(tangent);
        if f0.abs() <= 1e-13 * scale * scale && f1.abs() <= 1e-13 * scale {
            return Some(s);
        }
        let su = surface.partial_u(u, v);
        let sv = surface.partial_v(u, v);
        let (a, b) = (2.0 * r.dot(&su), 2.0 * r.dot(&sv));
        let (c, d) = (su.dot(tangent), sv.dot(tangent));
        let det = a * d - b * c;
        if det.abs() < 1e-300 {
            return None;
        }
        u -= (d * f0 - b * f1) / det;
        v -= (-c * f0 + a * f1) / det;
    }
    let s = surface.point_at(u, v);
    let r = s - *edge_point;
    ((r.length() - distance).abs() <= 1e-9 * scale && r.dot(tangent).abs() <= 1e-9 * scale)
        .then_some(s)
}
