use super::*;

#[test]
fn curve_turn_of_a_circular_arc_is_its_subtended_angle() {
    let r = 0.04;
    for &(a0, a1) in &[(0.0, 0.5 * std::f64::consts::PI), (0.3, 0.5), (1.0, 1.01)] {
        let at = |t: f64| Point3::new(r * f64::cos(t), r * f64::sin(t), 0.0);
        let turn = curve_turn(at(a0), at(0.5 * (a0 + a1)), at(a1));
        assert!((turn - (a1 - a0)).abs() < 1e-9, "{turn} vs {}", a1 - a0);
    }
}

#[test]
fn curve_turn_of_a_line_is_zero() {
    let turn = curve_turn(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(2.0, 2.0, 2.0),
    );
    assert!(turn.abs() < 1e-7);
}

#[test]
fn normal_spread_of_a_cylinder_cell_is_its_angular_width() {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.1,
    };
    let spread = normal_spread(
        (
            &cyl,
            &popkern_core::geometry::chart::ChartDescriptor::of(
                &cyl,
                &popkern_core::tolerance::Tolerance::default(),
            ),
        ),
        (0.0, 0.4),
        (0.0, 5.0),
    );
    assert!((spread - 0.4).abs() < 1e-9, "{spread}");
}

/// The coil_spring defect in closed form (PK-587): a chord along
/// the TOP circle of a thin torus has identical vertex normals
/// (spread 0), yet its midpoint lies over the tube `atan(s / r)`
/// around from the top, `s = R (1 − cos(Δu/2))` the chord's
/// sagitta.
#[test]
fn chord_foot_error_sees_a_torus_top_chord_that_spread_cannot() {
    let (big_r, r) = (1.0, 0.05);
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: big_r,
        minor_radius: r,
    };
    let v = 0.5 * std::f64::consts::PI;
    for du in [0.05_f64, 0.1, 0.2] {
        let (pa, pb) = (torus.point_at(0.0, v), torus.point_at(du, v));
        let n = |u: f64| {
            crate::vertex_normal::limit_normal_at(
                &popkern_core::geometry::chart::ChartDescriptor::of(
                    &torus,
                    &popkern_core::tolerance::Tolerance::default(),
                ),
                &torus,
                (u, v),
            )
            .expect("normal")
        };
        let (na, nb) = (n(0.0), n(du));
        assert!(angle_between(na, nb) < 1e-12);
        let err = chord_foot_error(
            (
                &torus,
                &popkern_core::geometry::chart::ChartDescriptor::of(
                    &torus,
                    &popkern_core::tolerance::Tolerance::default(),
                ),
            ),
            (0.5 * du, v),
            (pa, pb),
            (na, nb),
        )
        .expect("foot");
        let expected = (big_r * (1.0 - (0.5 * du).cos()) / r).atan();
        assert!(
            (err - expected).abs() < 1e-4 * expected,
            "du {du}: {err} vs {expected}"
        );
    }
}

/// Radial sagitta is invisible to the measure: a chord across a
/// cylinder, straight or diagonal, has its foot on the radial
/// projection, whose normal is the shaded one. Measured on an EXACT
/// rational quarter-cylinder (the evaluated path, in a non-angular
/// chart), this is what licenses the analytic short-circuit. The
/// sphere is the same identity in closed form: the shaded direction
/// `(pa − c) + (pb − c) = 2 (m − c)` is the foot normal.
#[test]
fn chord_foot_error_is_zero_on_an_exact_rational_cylinder() {
    let (r, w) = (0.1, std::f64::consts::FRAC_1_SQRT_2);
    let cps = [(r, 0.0), (r, r), (0.0, r)];
    let nurbs = popkern_core::geometry::NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: cps
            .iter()
            .flat_map(|&(x, y)| [Point3::new(x, y, 0.0), Point3::new(x, y, 2.0)])
            .collect(),
        num_u: 3,
        num_v: 2,
        weights: vec![1.0, 1.0, w, w, 1.0, 1.0],
    };
    let s = SurfaceGeometry::Nurbs(nurbs);
    for (a, b) in [((0.1, 0.0), (0.6, 0.0)), ((0.1, 0.05), (0.7, 0.95))] {
        let (pa, pb) = (s.point_at(a.0, a.1), s.point_at(b.0, b.1));
        assert!((pa.x.hypot(pa.y) - r).abs() < 1e-12, "exact cylinder");
        let n = |uv: (f64, f64)| {
            crate::vertex_normal::limit_normal_at(
                &popkern_core::geometry::chart::ChartDescriptor::of(
                    &s,
                    &popkern_core::tolerance::Tolerance::default(),
                ),
                &s,
                (uv.0, uv.1),
            )
            .expect("normal")
        };
        let mid = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
        let err = chord_foot_error(
            (
                &s,
                &popkern_core::geometry::chart::ChartDescriptor::of(
                    &s,
                    &popkern_core::tolerance::Tolerance::default(),
                ),
            ),
            mid,
            (pa, pb),
            (n(a), n(b)),
        )
        .expect("foot");
        // Two Gauss–Newton steps leave <= 1.3e-4 rad of truncation
        // in this non-angular chart (six leave 2e-8): 130x below
        // the smallest foot share in use (1.5° at a 3° bound).
        assert!(err < 2e-4, "{err}");
    }
}

#[test]
fn normal_spread_of_a_plane_is_zero() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(
        normal_spread(
            (
                &plane,
                &popkern_core::geometry::chart::ChartDescriptor::of(
                    &plane,
                    &popkern_core::tolerance::Tolerance::default()
                )
            ),
            (-3.0, 3.0),
            (-3.0, 3.0)
        ) < 1e-12
    );
}
