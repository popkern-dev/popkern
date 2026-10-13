use super::NurbsSurface;
use crate::geometry::{Point3, SurfaceEval};

/// The sub-patch hull contains every point of the patch and is smaller
/// than the whole surface's hull.
#[test]
fn sub_patch_hull_contains_its_patch() {
    // A bicubic wavy sheet, 6x6 control points.
    let n = 6;
    let mut control_points = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let (x, y) = (i as f64, j as f64);
            control_points.push(Point3::new(x, y, (x * 1.3).sin() * (y * 0.7).cos()));
        }
    }
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0, 1.0];
    let surface = NurbsSurface {
        degree_u: 3,
        degree_v: 3,
        knots_u: knots.clone(),
        knots_v: knots,
        control_points,
        weights: vec![1.0; n * n],
        num_u: n,
        num_v: n,
    };
    let (u, v) = ((0.2, 0.45), (0.6, 0.9));
    let hull = surface.sub_patch_hull_points(u, v).expect("hull");
    let lo = hull
        .iter()
        .fold(Point3::new(f64::MAX, f64::MAX, f64::MAX), |a, p| {
            Point3::new(a.x.min(p.x), a.y.min(p.y), a.z.min(p.z))
        });
    let hi = hull
        .iter()
        .fold(Point3::new(f64::MIN, f64::MIN, f64::MIN), |a, p| {
            Point3::new(a.x.max(p.x), a.y.max(p.y), a.z.max(p.z))
        });
    for i in 0..=20 {
        for j in 0..=20 {
            let p = surface.point_at(
                u.0 + (u.1 - u.0) * f64::from(i) / 20.0,
                v.0 + (v.1 - v.0) * f64::from(j) / 20.0,
            );
            for (x, a, b) in [(p.x, lo.x, hi.x), (p.y, lo.y, hi.y), (p.z, lo.z, hi.z)] {
                assert!(
                    x >= a - 1e-12 && x <= b + 1e-12,
                    "patch point outside its hull"
                );
            }
        }
    }
    assert!(
        hi.x - lo.x < 0.5 * 5.0 && hi.y - lo.y < 0.5 * 5.0,
        "trim-scale, not the whole sheet"
    );
}
