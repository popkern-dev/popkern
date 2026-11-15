use super::*;
use crate::geometry::Vec3;

/// The biquadratic NURBS sphere of radius 1 scaled by `s` along each
/// axis (the bench's ellipsoid net): nine-point rational circles in u,
/// a five-point rational semicircle from pole to pole in v.
fn ellipsoid(s: [f64; 3]) -> NurbsSurface {
    let w = 0.5_f64.sqrt();
    let ci = [
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 0.0],
        [-1.0, -1.0],
        [0.0, -1.0],
        [1.0, -1.0],
        [1.0, 0.0],
    ];
    let cw = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let mer = [[0.0, -1.0], [1.0, -1.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let mw = [1.0, w, 1.0, w, 1.0];
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for i in 0..9 {
        for j in 0..5 {
            control_points.push(Point3::new(
                s[0] * mer[j][0] * ci[i][0],
                s[1] * mer[j][0] * ci[i][1],
                s[2] * mer[j][1],
            ));
            weights.push(cw[i] * mw[j]);
        }
    }
    NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        knots_v: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
        control_points,
        num_u: 9,
        num_v: 5,
        weights,
    }
}

fn plane(origin: Point3, normal: Vec3) -> SurfaceGeometry {
    let u_axis = if normal.z.abs() < 0.9 {
        Vec3::Z.cross(&normal).normalized().expect("axis")
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let v_axis = normal.cross(&u_axis);
    SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis,
        v_axis,
    }
}

const WIDE: Window = ((-1e3, 1e3), (-1e3, 1e3));
const CHART: Window = ((0.0, 1.0), (0.0, 1.0));

/// One closed curve, on the ellipse `x²/a² + y²/b² = 1 − z²/c²` at
/// height `z`, both pcurves mapping onto it at its own parameter.
fn assert_section(
    found: &[SurfaceSurfaceIntersection],
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    s: [f64; 3],
    z: f64,
) {
    assert_eq!(found.len(), 1, "one section");
    let SurfaceSurfaceIntersection::Curve(c) = &found[0] else {
        panic!("a curve");
    };
    assert!(c.is_closed);
    let (t0, t1) = c.curve_3d.parameter_range();
    let k2 = 1.0 - (z / s[2]) * (z / s[2]);
    for k in 0..=200 {
        let t = t0 + (t1 - t0) * f64::from(k) / 200.0;
        let p = c.curve_3d.point_at(t);
        assert!((p.z - z).abs() < 1e-12, "height {}", p.z);
        let e = (p.x / s[0]).powi(2) + (p.y / s[1]).powi(2) - k2;
        assert!(e.abs() < 1e-12, "off the ellipse by {e}");
        for (uv, surf) in [(&c.curve_uv1, s1), (&c.curve_uv2, s2)] {
            let q = uv.point_at(t);
            assert!(surf.point_at(q.x, q.y).distance_to(&p) < 1e-9);
        }
    }
}

#[test]
fn an_ellipsoid_cut_square_to_its_axis_meets_the_plane_in_one_iso_curve() {
    let tol = Tolerance::default();
    let s = [6.0, 4.0, 3.0];
    let e = SurfaceGeometry::Nurbs(ellipsoid(s));
    for z in [0.0, 0.5, -0.7, -0.3, 2.9] {
        let p = plane(Point3::new(1.0, -2.0, z), Vec3::Z);
        let found = try_iso_level_plane(&e, CHART, &p, WIDE, &tol).expect("level rows");
        assert_section(&found, &e, &p, s, z);
        // Either argument order.
        let found = try_iso_level_plane(&p, WIDE, &e, CHART, &tol).expect("level rows");
        assert_section(&found, &p, &e, s, z);
    }
}

#[test]
fn a_plane_past_the_ellipsoid_gives_no_curve_and_one_at_a_pole_is_left_to_contacts() {
    let tol = Tolerance::default();
    let e = SurfaceGeometry::Nurbs(ellipsoid([6.0, 4.0, 3.0]));
    let past = plane(Point3::new(0.0, 0.0, 3.5), Vec3::Z);
    let found = try_iso_level_plane(&e, CHART, &past, WIDE, &tol).expect("level rows");
    assert!(found.is_empty(), "{}", found.len());
    // Touching at a pole the root lies on the chart's boundary row: the
    // arm declines and contact discovery types the touch.
    for z in [3.0, -3.0] {
        let p = plane(Point3::new(0.0, 0.0, z), Vec3::Z);
        assert!(
            try_iso_level_plane(&e, CHART, &p, WIDE, &tol).is_none(),
            "z {z}"
        );
    }
}

#[test]
fn a_plane_square_to_rows_that_are_not_level_keeps_the_marcher() {
    let tol = Tolerance::default();
    let e = SurfaceGeometry::Nurbs(ellipsoid([6.0, 4.0, 3.0]));
    let tilted = plane(
        Point3::new(0.0, 0.0, 0.5),
        Vec3::new(0.3, 0.0, 1.0).normalized().expect("n"),
    );
    assert!(try_iso_level_plane(&e, CHART, &tilted, WIDE, &tol).is_none());
    let side = plane(Point3::new(1.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
    assert!(try_iso_level_plane(&e, CHART, &side, WIDE, &tol).is_none());
}
