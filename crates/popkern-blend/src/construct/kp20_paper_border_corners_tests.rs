use super::*;
use crate::center_curve_helpers::fit_nurbs_through_points;

/// Degree-1 polyline through `samples` points of the line
/// `origin + t·dir`, `t ∈ [0, length]` — the shape of the G2 sweep's
/// contact trims (walker samples on a straight contact line).
fn polyline(origin: Point3, dir: Vec3, length: f64, samples: usize) -> NurbsCurve {
    let points: Vec<Point3> = (0..samples)
        .map(|i| origin + dir * (length * i as f64 / (samples - 1) as f64))
        .collect();
    fit_nurbs_through_points(&points).expect("polyline")
}

fn unit(angle: f64) -> Vec3 {
    Vec3::new(angle.cos(), angle.sin(), 0.0)
}

/// The blunt-pentagon class: two coplanar contact lines meeting at a
/// shallow angle, both spanning ~3 units with the crossing 0.2 from
/// one end. The former window scan stalled at 3.4e-3 here; the
/// Gauss–Newton polish must land on the exact crossing.
#[test]
fn shallow_angle_coplanar_lines_close_to_the_exact_crossing() {
    let crossing = Point3::new(0.3, -0.2, 1.7);
    let half = 10.5_f64.to_radians();
    let dir_a = unit(std::f64::consts::FRAC_PI_2 + half);
    let dir_b = unit(std::f64::consts::FRAC_PI_2 - half);
    // Each line runs from 0.2 before the crossing to 2.8 past it.
    let a = polyline(crossing - dir_a * 0.2, dir_a, 3.0, 20);
    let b = polyline(crossing - dir_b * 0.2, dir_b, 3.0, 20);
    let (gap, corner) = closest_between_refined(&a, &b, &|_| true);
    assert!(gap < 1.0e-10, "gap {gap:.3e}");
    assert!(
        corner.distance_to(&crossing) < 1.0e-9,
        "corner {corner:?} vs crossing {crossing:?}"
    );
}

/// Skew pair (no intersection exists — the octagon8_mixed class): the
/// polish reports the true line separation and the midpoint between
/// the closest points.
#[test]
fn skew_lines_report_their_true_separation() {
    let crossing = Point3::new(0.0, 0.0, 0.0);
    let half = 10.5_f64.to_radians();
    let dir_a = unit(std::f64::consts::FRAC_PI_2 + half);
    let dir_b = unit(std::f64::consts::FRAC_PI_2 - half);
    let lift = Vec3::new(0.0, 0.0, 0.01);
    let a = polyline(crossing - dir_a * 0.2, dir_a, 3.0, 20);
    let b = polyline(crossing - dir_b * 0.2 + lift, dir_b, 3.0, 20);
    let (gap, corner) = closest_between_refined(&a, &b, &|_| true);
    assert!((gap - 0.01).abs() < 1.0e-10, "gap {gap:.6e}");
    assert!(
        corner.distance_to(&(crossing + lift * 0.5)) < 1.0e-9,
        "corner {corner:?}"
    );
}

/// The material-side gate still binds: a crossing on the rejected
/// side is never returned, and the accepted answer stays in-domain.
#[test]
fn rejected_side_crossings_are_not_returned() {
    let crossing = Point3::new(0.0, 0.0, 0.0);
    let dir_a = unit(1.2);
    let dir_b = unit(-1.2);
    let a = polyline(crossing - dir_a * 1.0, dir_a, 2.0, 20);
    let b = polyline(crossing - dir_b * 1.0, dir_b, 2.0, 20);
    let accept = |p: Point3| p.x > 0.05;
    let (gap, corner) = closest_between_refined(&a, &b, &accept);
    assert!(accept(corner), "corner {corner:?} on the rejected side");
    assert!(gap > 0.0 && gap.is_finite());
}
