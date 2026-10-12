//! Conservative boxes around edge curves, for exact prescreens.
//!
//! `CurveEval::curve_bounding_box` samples the curve, so it can miss an
//! extreme between samples and cannot reject. These boxes contain the
//! whole curve by construction: a line segment's ends, a circle's or
//! ellipse's full extent, a positive-weight NURBS curve's control-point
//! hull. Curves without such a bound get none, and the caller tests them.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::{Point3, Vec3};

fn from_points(points: impl IntoIterator<Item = Point3>) -> Option<Aabb> {
    let mut points = points.into_iter();
    let first = points.next()?;
    let (mut min, mut max) = (first, first);
    for p in points {
        min = Point3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
        max = Point3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
    }
    Some(Aabb { min, max })
}

fn around(center: Point3, half: Vec3) -> Aabb {
    Aabb {
        min: Point3::new(center.x - half.x, center.y - half.y, center.z - half.z),
        max: Point3::new(center.x + half.x, center.y + half.y, center.z + half.z),
    }
}

/// A box containing `curve` over `t_range`, or `None` where no
/// conservative bound is cheap (parabolas, hyperbolas, NURBS with a
/// non-positive weight).
pub(crate) fn conservative_curve_box(curve: &CurveGeometry, t_range: (f64, f64)) -> Option<Aabb> {
    match curve {
        CurveGeometry::Line { origin, direction } => {
            from_points([t_range.0, t_range.1].map(|t| *origin + *direction * t))
        }
        CurveGeometry::Circle { center, radius, .. }
        | CurveGeometry::Arc { center, radius, .. } => {
            Some(around(*center, Vec3::new(*radius, *radius, *radius)))
        }
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            ..
        } => Some(around(
            *center,
            Vec3::new(
                u_basis.x.abs() + v_basis.x.abs(),
                u_basis.y.abs() + v_basis.y.abs(),
                u_basis.z.abs() + v_basis.z.abs(),
            ),
        )),
        CurveGeometry::Nurbs(nurbs) => {
            if nurbs.weights.iter().any(|&w| w <= 0.0) {
                return None;
            }
            from_points(nurbs.control_points.iter().copied())
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::curve::CurveEval;

    fn contains(bx: &Aabb, p: Point3) -> bool {
        let e = 1e-12;
        p.x >= bx.min.x - e
            && p.y >= bx.min.y - e
            && p.z >= bx.min.z - e
            && p.x <= bx.max.x + e
            && p.y <= bx.max.y + e
            && p.z <= bx.max.z + e
    }

    fn assert_contains_curve(curve: &CurveGeometry, t_range: (f64, f64)) {
        let bx = conservative_curve_box(curve, t_range).expect("box");
        for k in 0..=1000 {
            let t = t_range.0 + (t_range.1 - t_range.0) * k as f64 / 1000.0;
            assert!(
                contains(&bx, curve.point_at(t)),
                "{curve:?} leaves its box at t = {t}"
            );
        }
    }

    #[test]
    fn boxes_contain_their_curves() {
        let line = CurveGeometry::Line {
            origin: Point3::new(1.0, -2.0, 0.5),
            direction: Vec3::new(0.3, 0.4, -1.2),
        };
        assert_contains_curve(&line, (-2.0, 3.0));
        // An arc whose x extreme (angle 0) lies between any coarse
        // sampling of its span: the box must still reach x = 2.
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
            start_angle: -0.7,
            end_angle: 0.9,
        };
        assert_contains_curve(&arc, (-0.7, 0.9));
        let ellipse = CurveGeometry::Ellipse {
            center: Point3::new(1.0, 1.0, 1.0),
            u_basis: Vec3::new(3.0, 1.0, 0.0),
            v_basis: Vec3::new(0.0, 0.5, 2.0),
            start_angle: 0.0,
            end_angle: std::f64::consts::TAU,
        };
        assert_contains_curve(&ellipse, (0.0, std::f64::consts::TAU));
        let nurbs = popkern_core::geometry::NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 2.0, 0.0),
                Point3::new(2.0, 0.0, 1.0),
            ],
            weights: vec![1.0, 0.5, 1.0],
        };
        let curve = CurveGeometry::Nurbs(nurbs);
        assert_contains_curve(&curve, (0.0, 1.0));
    }

    #[test]
    fn unbounded_kinds_get_no_box() {
        let nurbs = popkern_core::geometry::NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 0.0],
        };
        assert!(conservative_curve_box(&CurveGeometry::Nurbs(nurbs), (0.0, 1.0)).is_none());
    }
}
