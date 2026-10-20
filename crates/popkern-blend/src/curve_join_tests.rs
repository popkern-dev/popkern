//! G1/G2/G3 joins, reparameterization-aware (PK-507).

use super::*;
use popkern_core::geometry::NurbsCurve;

const TOL: JoinTolerances = JoinTolerances {
    gap: 1e-9,
    angle: 1e-9,
    curvature: 1e-9,
    curvature_derivative: 1e-7,
};

/// The rational quadratic arc of radius `r` about the origin from
/// angle `a0` to `a1` (at most a half turn), in z = 0.
fn arc(r: f64, a0: f64, a1: f64) -> CurveGeometry {
    let half = 0.5 * (a1 - a0);
    let w = half.pcos();
    let mid = 0.5 * (a0 + a1);
    let p = |a: f64, s: f64| Point3::new(s * r * a.pcos(), s * r * a.psin(), 0.0);
    CurveGeometry::Nurbs(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![p(a0, 1.0), p(mid, 1.0 / w), p(a1, 1.0)],
        weights: vec![1.0, w, 1.0],
    })
}

fn line(a: Point3, b: Point3) -> CurveGeometry {
    CurveGeometry::Nurbs(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    })
}

#[test]
fn two_arcs_of_one_circle_join_g3() {
    let (a, b) = (arc(2.0, 0.0, 1.2), arc(2.0, 1.2, 2.0));
    let r = analyze_join(&JoinJet::of(&a, 1.0), &JoinJet::of(&b, 0.0), &TOL);
    assert_eq!(r.continuity, 3, "{r:?}");
    assert!(r.g2_residual < 1e-12 && r.g3_residual < 1e-10, "{r:?}");
}

#[test]
fn equal_radii_bending_opposite_ways_are_g1_not_g2() {
    // An S: the arc about the origin ends at (0, 2) heading -x; the
    // next arc, about (0, 4), starts there heading -x and bends the
    // other way.
    let a = arc(2.0, 0.0, std::f64::consts::FRAC_PI_2);
    let b = {
        let CurveGeometry::Nurbs(n) = arc(2.0, -std::f64::consts::FRAC_PI_2, -std::f64::consts::PI)
        else {
            unreachable!()
        };
        CurveGeometry::Nurbs(NurbsCurve {
            control_points: n
                .control_points
                .iter()
                .map(|p| Point3::new(p.x, p.y + 4.0, p.z))
                .collect(),
            ..n
        })
    };
    let r = analyze_join(&JoinJet::of(&a, 1.0), &JoinJet::of(&b, 0.0), &TOL);
    assert!(r.gap < 1e-12, "{r:?}");
    assert_eq!(r.continuity, 1, "{r:?}");
    // |kappa| is 1/2 on both sides; the vectors oppose: the jump is 1.
    assert!((r.curvature_jump - 1.0).abs() < 1e-12, "{r:?}");
}

#[test]
fn an_arc_into_its_tangent_line_is_g1() {
    let a = arc(2.0, 0.0, std::f64::consts::FRAC_PI_2);
    let b = line(Point3::new(0.0, 2.0, 0.0), Point3::new(-3.0, 2.0, 0.0));
    let r = analyze_join(&JoinJet::of(&a, 1.0), &JoinJet::of(&b, 0.0), &TOL);
    assert_eq!(r.continuity, 1, "{r:?}");
    assert!((r.curvature_jump - 0.5).abs() < 1e-12, "{r:?}");
}

/// A cubic's jet at u, and the jet of the same curve reparameterized
/// by t(s) with rates (beta, gamma, delta) at the join.
fn cubic_jets(u: f64, beta: f64, gamma: f64, delta: f64) -> (JoinJet, JoinJet) {
    let c = |k: usize| -> Vec3 {
        let (a, b, c2, d) = (
            Vec3::new(3.0, -1.0, 2.0),
            Vec3::new(-2.0, 4.0, 1.0),
            Vec3::new(0.5, -3.0, 2.5),
            Vec3::new(1.0, 2.0, 0.5),
        );
        match k {
            1 => a + b * (2.0 * u) + c2 * (3.0 * u * u),
            2 => b * 2.0 + c2 * (6.0 * u),
            3 => c2 * 6.0,
            _ => d + a * u + b * (u * u) + c2 * (u * u * u),
        }
    };
    let p = c(0);
    let end = JoinJet {
        point: Point3::new(p.x, p.y, p.z),
        d: [c(1), c(2), c(3)],
    };
    let start = JoinJet {
        point: end.point,
        d: [
            c(1) * beta,
            c(2) * (beta * beta) + c(1) * gamma,
            c(3) * beta.powi(3) + c(2) * (3.0 * beta * gamma) + c(1) * delta,
        ],
    };
    (end, start)
}

#[test]
fn a_nonlinear_reparameterization_is_still_g3_and_its_rates_recovered() {
    let (end, start) = cubic_jets(0.4, 2.0, 0.7, -0.3);
    let r = analyze_join(&end, &start, &TOL);
    assert_eq!(r.continuity, 3, "{r:?}");
    assert!((r.beta - 2.0).abs() < 1e-12 && (r.gamma - 0.7).abs() < 1e-12);
    assert!((r.delta + 0.3).abs() < 1e-12, "{r:?}");
    // A third derivative off the identity breaks G3, not G2.
    let mut broken = start;
    broken.d[2] = broken.d[2] + Vec3::new(0.0, 0.0, 1e-3);
    let r = analyze_join(&end, &broken, &TOL);
    assert_eq!(r.continuity, 2, "{r:?}");
    assert!(r.g3_residual > 1e-6 && r.g2_residual < 1e-12, "{r:?}");
}

#[test]
fn opposed_tangents_and_gaps_fail_early() {
    let (end, mut start) = cubic_jets(0.4, 1.0, 0.0, 0.0);
    start.d[0] = start.d[0] * -1.0;
    assert_eq!(analyze_join(&end, &start, &TOL).continuity, 0);
    let (end, mut start) = cubic_jets(0.4, 1.0, 0.0, 0.0);
    start.point = Point3::new(start.point.x + 1e-3, start.point.y, start.point.z);
    assert_eq!(analyze_join(&end, &start, &TOL).continuity, -1);
}
