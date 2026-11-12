// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

//! The off-axis cylinder×sphere pierce: the branch count of each regime
//! and every curve, and both its pcurves' images, on both surfaces.

use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::cylinder_sphere::try_cylinder_sphere;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::{PI, TAU};

fn tol() -> Tolerance {
    Tolerance::default().with_linear(1e-6).with_angular(1e-6)
}

/// A cylinder of radius 6 along y through (x, 0, 0) against a ball of
/// radius 20 on the origin (its pole along `pole`).
fn pierce(x: f64, pole: Vec3) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(x, -30.0, 0.0),
        axis: Vec3::new(0.0, 1.0, 0.0),
        radius: 6.0,
    };
    let ball = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 20.0,
        axis: pole,
    };
    try_cylinder_sphere(
        &cyl,
        ((0.0, TAU), (0.0, 60.0)),
        &ball,
        ((0.0, TAU), (0.0, PI)),
        &tol(),
    )
}

/// Every curve lies on both surfaces, and so do its pcurves' images.
fn assert_on_both(curves: &[SurfaceSurfaceIntersection], x: f64, pole: Vec3) {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(x, -30.0, 0.0),
        axis: Vec3::new(0.0, 1.0, 0.0),
        radius: 6.0,
    };
    let ball = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 20.0,
        axis: pole,
    };
    for c in curves {
        let SurfaceSurfaceIntersection::Curve(c) = c else {
            panic!("a point where a curve belongs");
        };
        let (t0, t1) = c.curve_3d.parameter_range();
        let (a0, a1) = c.curve_uv1.parameter_range();
        let (b0, b1) = c.curve_uv2.parameter_range();
        for k in 0..=400 {
            let f = k as f64 / 400.0;
            let p = c.curve_3d.point_at(t0 + (t1 - t0) * f);
            let rho = ((p.x - x).powi(2) + p.z * p.z).sqrt();
            assert!(
                (rho - 6.0).abs() < 1e-6,
                "off the cylinder by {}",
                rho - 6.0
            );
            let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!((r - 20.0).abs() < 1e-6, "off the ball by {}", r - 20.0);
            let q1 = c.curve_uv1.point_at(a0 + (a1 - a0) * f);
            let q2 = c.curve_uv2.point_at(b0 + (b1 - b0) * f);
            assert!(cyl.point_at(q1.x, q1.y).distance_to(&p) < 1e-6);
            assert!(ball.point_at(q2.x, q2.y).distance_to(&p) < 1e-6);
        }
    }
}

#[test]
fn a_bore_through_the_ball_off_its_centre_meets_it_in_two_loops() {
    let pole = Vec3::new(0.0, 0.0, 1.0);
    for x in [0.5, 5.0, 8.0] {
        let curves = pierce(x, pole).expect("answered across the pole axis");
        assert_eq!(curves.len(), 2, "x {x}");
        assert_on_both(&curves, x, pole);
    }
}

/// Along the pole axis a bore holding the pole winds its loops round
/// the chart's poles: the scan leaves those to the marching path. One
/// beside the pole (x 8, the pole 2 outside its wall) stays in the
/// chart and is answered on both surfaces.
#[test]
fn a_bore_round_the_pole_is_left_to_the_marching_path() {
    let pole = Vec3::new(0.0, 1.0, 0.0);
    for x in [0.5, 5.0] {
        assert!(pierce(x, pole).is_none(), "x {x}");
    }
    let curves = pierce(8.0, pole).expect("beside the pole");
    assert_eq!(curves.len(), 2);
    assert_on_both(&curves, 8.0, pole);
}

/// Breaking out of the ball's side (|x| + r > R): one loop round the
/// break-out.
#[test]
fn a_bore_breaking_out_of_the_side_meets_it_in_one_loop() {
    let curves = pierce(16.0, Vec3::new(0.0, 0.0, 1.0)).expect("answered");
    assert_eq!(curves.len(), 1);
    assert_on_both(&curves, 16.0, Vec3::new(0.0, 0.0, 1.0));
}

#[test]
fn a_bore_beside_the_ball_misses_it() {
    assert!(
        pierce(27.0, Vec3::new(0.0, 0.0, 1.0))
            .expect("answered")
            .is_empty()
    );
}

/// Viviani's curve with the ball's pole at its double point: a cylinder
/// of radius 5 on x = 5 inside a ball of radius 10, its pole on x. The
/// branches end at the pole (10, 0, 0), where the chart's longitude is
/// arbitrary; they take its limit along the branch, and the pierce
/// serves the two lobes (each root over the whole turn, from the double
/// point back to it) exactly, every pcurve image on its surface
/// (PK-1139: the arbitrary longitude read as a pole sweep and declined
/// to the marcher, whose records stopped 7.5e-5 short).
#[test]
fn vivianis_curve_through_the_balls_pole_serves() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(5.0, 0.0, -12.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let ball = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 10.0,
        axis: Vec3::new(1.0, 0.0, 0.0),
    };
    let found = try_cylinder_sphere(
        &cyl,
        ((0.0, TAU), (0.0, 24.0)),
        &ball,
        ((0.0, TAU), (0.0, PI)),
        &tol(),
    )
    .expect("the pierce serves");
    assert_eq!(found.len(), 2, "a lobe per root");
    let pole = Point3::new(10.0, 0.0, 0.0);
    for f in &found {
        let SurfaceSurfaceIntersection::Curve(c) = f else {
            panic!("a curve");
        };
        let (t0, t1) = c.curve_3d.parameter_range();
        let ends = [c.curve_3d.point_at(t0), c.curve_3d.point_at(t1)];
        assert!(ends.iter().all(|e| e.distance_to(&pole) < 1e-9), "{ends:?}");
        for k in 0..=200 {
            let t = t0 + (t1 - t0) * f64::from(k) / 200.0;
            let p = c.curve_3d.point_at(t);
            let (q1, q2) = (c.curve_uv1.point_at(t), c.curve_uv2.point_at(t));
            assert!(
                cyl.point_at(q1.x, q1.y).distance_to(&p) < 1e-5,
                "uv1 at {t}"
            );
            assert!(
                ball.point_at(q2.x, q2.y).distance_to(&p) < 1e-5,
                "uv2 at {t}"
            );
        }
    }
}

/// Viviani's operands with the ball's poles on the curve: a cylinder of
/// radius 5 on x = 5 inside a ball of radius 10 whose pole axis is z,
/// in the order `swapped` says.
fn viviani_poles_on_the_curve(swapped: bool) -> Vec<SurfaceSurfaceIntersection> {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(5.0, 0.0, -12.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let ball = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 10.0,
        axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let (c, b) = (((0.0, TAU), (0.0, 24.0)), ((0.0, TAU), (0.0, PI)));
    if swapped {
        try_cylinder_sphere(&ball, b, &cyl, c, &tol())
    } else {
        try_cylinder_sphere(&cyl, c, &ball, b, &tol())
    }
    .expect("the pierce serves")
}

/// Each lobe of Viviani's figure eight passes through a pole of the
/// ball, where the sphere chart's longitude turns half a turn: the scan
/// cuts it there into two arcs from the double point to the pole, four
/// in all, every pcurve image within a tenth of the tolerance of the
/// curve, and each arc's sphere pcurve ends at the pole at the
/// longitude the arc leaves it along (the chart's longitude of the arc
/// a micro-unit away). The branch declined to the marching path as a
/// sweep round the pole, which traced nine records, near-duplicate arcs
/// 2e-3 apart and loops at the double point among them (PK-1169).
#[test]
fn vivianis_curve_with_the_balls_poles_on_it_is_cut_at_them() {
    let ball = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 10.0,
        axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(5.0, 0.0, -12.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let double = Point3::new(10.0, 0.0, 0.0);
    let poles = [Point3::new(0.0, 0.0, 10.0), Point3::new(0.0, 0.0, -10.0)];
    for swapped in [false, true] {
        let found = viviani_poles_on_the_curve(swapped);
        assert_eq!(found.len(), 4, "an arc per lobe and side of its pole");
        let mut reached = [0, 0];
        for f in &found {
            let SurfaceSurfaceIntersection::Curve(c) = f else {
                panic!("a curve");
            };
            let (t0, t1) = c.curve_3d.parameter_range();
            let ends = [c.curve_3d.point_at(t0), c.curve_3d.point_at(t1)];
            let at_double = ends.iter().position(|e| e.distance_to(&double) < 1e-9);
            let at_pole = ends
                .iter()
                .position(|e| poles.iter().any(|p| e.distance_to(p) < 1e-9));
            let (Some(d), Some(k)) = (at_double, at_pole) else {
                panic!("an arc from the double point to a pole: {ends:?}");
            };
            assert_ne!(d, k);
            reached[usize::from(ends[k].z < 0.0)] += 1;
            let (on_cyl, on_ball) = if swapped {
                (&c.curve_uv2, &c.curve_uv1)
            } else {
                (&c.curve_uv1, &c.curve_uv2)
            };
            for i in 0..=400 {
                let t = t0 + (t1 - t0) * f64::from(i) / 400.0;
                let p = c.curve_3d.point_at(t);
                let (q1, q2) = (on_cyl.point_at(t), on_ball.point_at(t));
                assert!(
                    cyl.point_at(q1.x, q1.y).distance_to(&p) < 1e-7,
                    "on the rod at {t}"
                );
                assert!(
                    ball.point_at(q2.x, q2.y).distance_to(&p) < 1e-7,
                    "on the ball at {t}"
                );
            }
            // The pole end's longitude is the arc's limit there.
            let (t_pole, t_in) = if k == 0 {
                (t0, t0 + (t1 - t0) * 1e-7)
            } else {
                (t1, t1 - (t1 - t0) * 1e-7)
            };
            let end_u = on_ball.point_at(t_pole).x;
            let near = c.curve_3d.point_at(t_in);
            let limit = near.y.atan2(near.x);
            let turn = (end_u - limit).rem_euclid(TAU);
            assert!(
                turn.min(TAU - turn) < 1e-5,
                "the pole end's longitude {end_u} for the limit {limit}"
            );
        }
        assert_eq!(reached, [2, 2], "two arcs reach each pole");
    }
}

/// The fit is checked against each pcurve's own surface, the sphere's
/// first when the sphere is the caller's first surface: checked against
/// the cylinder there it never passed, and every such pierce (a boolean
/// passes the ball first when it is the first operand) densified to the
/// sample cap, 16384 a turn. Both orders now fit the same samples.
#[test]
fn both_operand_orders_fit_the_same_samples() {
    let count = |swapped: bool| -> Vec<usize> {
        viviani_poles_on_the_curve(swapped)
            .iter()
            .map(|f| match f {
                SurfaceSurfaceIntersection::Curve(c) => c.curve_3d.control_points.len(),
                _ => 0,
            })
            .collect()
    };
    let (straight, swapped) = (count(false), count(true));
    assert_eq!(straight, swapped);
    assert!(straight.iter().all(|&n| n < 4096), "{straight:?}");
}

/// A branch is cut at a pole it passes: an open branch with the pole
/// between two samples gains it there and is cut in two; a closed one
/// through a sample at the pole is opened there into one arc from the
/// pole round to it; a pole off the branch by more than a step leaves
/// it whole.
#[test]
fn a_branch_is_cut_at_each_pole_it_passes() {
    use super::split_at_poles;
    let line: Vec<Point3> = (0..=10)
        .map(|i| Point3::new(f64::from(i), 0.0, 0.0))
        .collect();
    let pole = Point3::new(4.5, 0.0, 0.0);
    let pieces = split_at_poles(line.clone(), false, &[pole]);
    assert_eq!(pieces.len(), 2);
    assert!(pieces.iter().all(|(_, closed)| !closed));
    assert_eq!(pieces[0].0.last(), Some(&pole));
    assert_eq!(pieces[1].0.first(), Some(&pole));
    assert_eq!(pieces[0].0.len() + pieces[1].0.len(), 13);
    let far = Point3::new(4.5, 3.0, 0.0);
    assert_eq!(split_at_poles(line, false, &[far]).len(), 1);
    // A closed square loop, its last sample its first.
    let ring = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    ];
    let corner = Point3::new(1.0, 1.0, 0.0);
    let pieces = split_at_poles(ring, true, &[corner]);
    assert_eq!(pieces.len(), 1);
    let (arc, closed) = &pieces[0];
    assert!(!closed);
    assert_eq!((arc.first(), arc.last()), (Some(&corner), Some(&corner)));
    assert_eq!(arc.len(), 5);
}
