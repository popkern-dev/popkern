//! The chart service's unit tests (design §9).

use super::descriptor::{Axis, ChartDescriptor, Continuity, PeriodKind, SingularKind};
use super::ops::{canonical, distortion_radius, lift};
use super::regular::{RegularChart, regular_chart};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

fn sphere() -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(1.0, -2.0, 0.5),
        radius: 2.0,
        axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

#[test]
fn every_analytic_kind_reports_its_periods_and_singular_rows() {
    let tol = Tolerance::default();
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let d = ChartDescriptor::of(&plane, &tol);
    assert!(d.u.is_none() && d.v.is_none() && d.singular_rows.is_empty());
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.5,
    };
    let d = ChartDescriptor::of(&cyl, &tol);
    let u = d.u.expect("cylinder closes along u");
    assert_eq!(
        (u.kind, u.continuity),
        (PeriodKind::Analytic, Continuity::Smooth)
    );
    assert!(d.v.is_none() && d.singular_rows.is_empty());
    let d = ChartDescriptor::of(&sphere(), &tol);
    assert!(d.u.is_some());
    assert_eq!(d.singular_rows.len(), 2, "{:?}", d.singular_rows);
    assert!(
        d.singular_rows
            .iter()
            .all(|r| r.kind == SingularKind::Chart && r.axis == Axis::V)
    );
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: 0.4,
    };
    let d = ChartDescriptor::of(&cone, &tol);
    assert!(d.u.is_some());
    assert_eq!(d.singular_rows.len(), 1, "{:?}", d.singular_rows);
    assert_eq!(
        d.singular_rows[0].kind,
        SingularKind::Apex,
        "a cone's apex is a geometric singularity"
    );
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 1.0,
    };
    let d = ChartDescriptor::of(&torus, &tol);
    assert!(d.u.is_some() && d.v.is_some() && d.singular_rows.is_empty());
}

#[test]
fn a_nurbs_ball_has_chart_singular_poles_and_a_closed_chart() {
    let tol = Tolerance::default();
    let s = sphere();
    let nurbs = SurfaceGeometry::Nurbs(
        s.to_nurbs_in_range((0.0, std::f64::consts::TAU), (0.0, std::f64::consts::PI)),
    );
    let d = ChartDescriptor::of(&nurbs, &tol);
    assert!(
        d.u.is_some() || d.v.is_some(),
        "the twin closes round its axis"
    );
    assert_eq!(d.singular_rows.len(), 2, "{:?}", d.singular_rows);
    assert!(
        d.singular_rows
            .iter()
            .all(|r| r.kind == SingularKind::Chart)
    );
}

#[test]
fn the_pole_projection_round_trips_beside_and_at_the_pole() {
    let tol = Tolerance::default();
    let s = sphere();
    let d = ChartDescriptor::of(&s, &tol);
    for &(u, v) in &[
        (0.3, 1e-6),
        (2.0, 1e-3),
        (5.5, 0.05),
        (1.0, std::f64::consts::PI - 1e-4),
    ] {
        let chart = regular_chart(&d, &s, (u, v), 0.2);
        assert!(
            matches!(chart, RegularChart::PoleProjection { .. }),
            "({u}, {v}): {chart:?}"
        );
        let xy = chart.from_surface(&s, (u, v));
        let back = chart.to_surface(&s, xy).expect("inverse converges");
        let (p, q) = (s.point_at(u, v), s.point_at(back.0, back.1));
        assert!(
            p.distance_to(&q) <= 1e-12,
            "({u}, {v}) -> {xy:?} -> {back:?}: {}",
            p.distance_to(&q)
        );
    }
}

#[test]
fn a_shifted_chart_moves_the_cut_half_a_period_away() {
    let tol = Tolerance::default();
    let s = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let d = ChartDescriptor::of(&s, &tol);
    let chart = regular_chart(&d, &s, (1e-4, 0.5), 0.01);
    assert!(matches!(chart, RegularChart::Shifted { .. }), "{chart:?}");
    // Either side of the cut reads as neighbours in the shifted chart.
    let a = chart.from_surface(&s, (std::f64::consts::TAU - 1e-4, 0.5));
    let b = chart.from_surface(&s, (1e-4, 0.5));
    assert!((a.0 - b.0).abs() < 1e-3, "{a:?} {b:?}");
    let back = chart.to_surface(&s, a).expect("shift");
    assert!(
        s.point_at(back.0, back.1)
            .distance_to(&s.point_at(std::f64::consts::TAU - 1e-4, 0.5))
            < 1e-12
    );
}

#[test]
fn a_lift_across_the_cut_is_continuous_and_counts_its_winding() {
    let tol = Tolerance::default();
    let s = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let d = ChartDescriptor::of(&s, &tol);
    let tau = std::f64::consts::TAU;
    let ring: Vec<(f64, f64)> = (0..=16)
        .map(|k| ((k as f64 * tau / 16.0 + 6.0).rem_euclid(tau), 0.0))
        .collect();
    let lifted = lift(&d, &ring);
    assert!(
        lifted
            .points
            .windows(2)
            .all(|w| (w[1].0 - w[0].0).abs() < 1.0)
    );
    assert_eq!(lifted.winding, (1, 0));
    let ((u, _), (k, _)) = canonical(&d, (tau + 0.25, 0.0));
    assert!((u - 0.25).abs() < 1e-12 && k == -1);
}

#[test]
fn the_distortion_radius_of_a_sphere_pole_is_its_closed_form() {
    let tol = Tolerance::default();
    let s = sphere();
    let d = ChartDescriptor::of(&s, &tol);
    let r = distortion_radius(&d, &s, &d.singular_rows[0]);
    // |S_v| / |S_u| = 1 / sin v = 10 at sin v = 0.1: the chord 2R sin(v/2).
    let v = (0.1f64).asin();
    let want = 2.0 * 2.0 * (v / 2.0).sin();
    assert!((r - want).abs() < 1e-6, "{r} vs {want}");
}

/// A charted surface over a pole projection: its points are the
/// surface's, its partials match finite differences of its points, and
/// its orientation (the local partials' cross product) agrees with the
/// surface normal, so a loop running counter-clockwise in the face chart
/// runs counter-clockwise in the cap chart too.
#[test]
fn a_charted_pole_projection_evaluates_the_surface_with_consistent_partials() {
    use super::charted::ChartedSurface;
    let tol = Tolerance::default();
    let s = sphere();
    let d = ChartDescriptor::of(&s, &tol);
    for at in [(0.3, 1e-6), (1.0, std::f64::consts::PI - 1e-4)] {
        let chart = regular_chart(&d, &s, at, 0.5);
        let c = ChartedSurface::new(&s, chart, ((-1.0, 1.0), (-1.0, 1.0)));
        for &(x, y) in &[(0.0, 0.0), (0.3, -0.2), (-0.6, 0.5)] {
            let p = c.point_at(x, y);
            let (u, v) = s.closest_point(&p).expect("on the sphere");
            assert!(s.point_at(u, v).distance_to(&p) <= 1e-12);
            let h = 1e-6;
            let fx = (c.point_at(x + h, y) - c.point_at(x - h, y)) * (0.5 / h);
            let fy = (c.point_at(x, y + h) - c.point_at(x, y - h)) * (0.5 / h);
            let (px, py) = c.first_partials(x, y);
            assert!((fx - px).length() <= 1e-6, "{x},{y}: {fx:?} vs {px:?}");
            assert!((fy - py).length() <= 1e-6, "{x},{y}: {fy:?} vs {py:?}");
            // The face chart's own orientation, read a step off the row
            // where its partials vanish.
            let (fu, fv) = c.face_uv(x, y).expect("inverse");
            let fv = fv.clamp(1e-3, std::f64::consts::PI - 1e-3);
            let (su, sv) = s.first_partials(fu, fv);
            assert!(
                px.cross(&py).dot(&su.cross(&sv)) > 0.0,
                "{at:?} {x},{y}: orientation flips"
            );
            let back = c.closest_point(&p).expect("projects");
            assert!(c.point_at(back.0, back.1).distance_to(&p) <= 1e-12);
        }
    }
}

/// Interpolation in a regular chart: two points either side of a sphere's
/// pole meet near the pole (the face-chart midpoint lands half a period
/// round, on the far side of the sphere), and two points either side of
/// the cut meet on the near side.
#[test]
fn interpolation_runs_through_the_pole_and_across_the_cut() {
    use super::ops::interpolate;
    let tol = Tolerance::default();
    let s = sphere();
    let d = ChartDescriptor::of(&s, &tol);
    for row in &d.singular_rows {
        assert!(row.distortion > 0.19 && row.distortion < 0.21, "{row:?}");
    }
    let at = |uv: (f64, f64)| (uv, s.point_at(uv.0, uv.1));
    // Either side of the north pole (v = 0), half a period apart in u.
    let (a, b) = (at((0.3, 0.02)), at((0.3 + std::f64::consts::PI, 0.02)));
    let (uv, p) = interpolate(&d, &s, a, b, 0.5).expect("converges");
    let pole = s.point_at(0.0, 0.0);
    assert!(
        p.distance_to(&pole) < 1e-9,
        "{uv:?} {}",
        p.distance_to(&pole)
    );
    // Either side of the cut, at the equator: the midpoint is at u = 0.
    let (a, b) = (at((0.05, 1.5)), at((std::f64::consts::TAU - 0.05, 1.5)));
    let (uv, p) = interpolate(&d, &s, a, b, 0.5).expect("face chart");
    assert!(p.distance_to(&s.point_at(0.0, 1.5)) < 1e-12, "{uv:?}");
}

/// A sphere's cap rows sit where the normal has turned 45°: colatitude
/// π/4 about the north pole and 3π/4 about the south.
#[test]
fn a_spheres_cap_rows_sit_at_a_quarter_turn() {
    use super::atlas::cap_row;
    let tol = Tolerance::default();
    let s = sphere();
    let d = ChartDescriptor::of(&s, &tol);
    assert_eq!(d.singular_rows.len(), 2);
    for row in &d.singular_rows {
        let cap = cap_row(&d, &s, row).expect("cap row");
        let want = if row.at < 1.0 {
            std::f64::consts::FRAC_PI_4
        } else {
            3.0 * std::f64::consts::FRAC_PI_4
        };
        assert!(!cap.whole);
        assert!((cap.at - want).abs() < 1e-5, "{row:?}: {}", cap.at);
    }
}
