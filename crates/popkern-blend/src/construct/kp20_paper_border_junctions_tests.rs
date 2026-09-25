use super::*;

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

fn strip(
    origin: Point3,
    axis: Vec3,
    radius: f64,
    rails: [(Point3, Point3); 2],
) -> FilletSurfaceResult {
    FilletSurfaceResult {
        surface: SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        },
        trim_curve_1: line(rails[0].0, rails[0].1),
        trim_curve_2: line(rails[1].0, rails[1].1),
        continuity_downgrade: None,
    }
}

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// The notched-cube saddle at (1, 0.6, 0.6). S0 is the concave strip
/// along the notch edge (r 0.06). S1 and S2 are the convex strips
/// along x = 1 (r 0.08 along y, r 0.10 along z).
fn notch() -> [FilletSurfaceResult; 3] {
    [
        strip(
            p(0.6, 0.66, 0.66),
            Vec3::new(1.0, 0.0, 0.0),
            0.06,
            [
                (p(1.0, 0.66, 0.6), p(0.6, 0.66, 0.6)),
                (p(1.0, 0.6, 0.66), p(0.6, 0.6, 0.66)),
            ],
        ),
        strip(
            p(0.92, 0.6, 0.52),
            Vec3::new(0.0, 1.0, 0.0),
            0.08,
            [
                (p(0.92, 0.6, 0.6), p(0.92, 1.0, 0.6)),
                (p(1.0, 0.6, 0.52), p(1.0, 1.0, 0.52)),
            ],
        ),
        strip(
            p(0.9, 0.5, 0.6),
            Vec3::new(0.0, 0.0, 1.0),
            0.10,
            [
                (p(0.9, 0.6, 0.6), p(0.9, 0.6, 1.0)),
                (p(1.0, 0.5, 0.6), p(1.0, 0.5, 1.0)),
            ],
        ),
    ]
}

const VERTEX: Point3 = Point3 {
    x: 1.0,
    y: 0.6,
    z: 0.6,
};

/// The floor and wall are convex face corners, and x = 1 is reflex.
#[test]
fn the_notch_saddle_has_one_reflex_face() {
    let [s0, s1, s2] = notch();
    assert_eq!(
        classify_junction(&s0, &s1, VERTEX),
        Ok(JunctionKind::Corner)
    );
    assert_eq!(
        classify_junction(&s1, &s2, VERTEX),
        Ok(JunctionKind::Bridge)
    );
    assert_eq!(
        classify_junction(&s2, &s0, VERTEX),
        Ok(JunctionKind::Corner)
    );
}

/// The bridge leaves S1's rail level with S1's floor corner, and
/// arrives on S2's rail level with S2's wall corner. It leaves along
/// rail S1 and arrives along rail S2, with zero curvature at both
/// ends (the strips' rulings are straight), and stays in the x = 1
/// plane.
#[test]
fn the_bridge_continues_both_rails_on_the_reflex_face() {
    let [_, s1, s2] = notch();
    let bridge = bridge_between(&s1, &s2, p(0.92, 0.66, 0.6), p(0.9, 0.6, 0.66), VERTEX)
        .expect("bridge builds");
    assert!(bridge.start.distance_to(&p(1.0, 0.66, 0.52)) < 1e-12);
    assert!(bridge.end.distance_to(&p(1.0, 0.5, 0.66)) < 1e-12);
    for (t, along) in [
        (0.0, Vec3::new(0.0, -1.0, 0.0)),
        (1.0, Vec3::new(0.0, 0.0, 1.0)),
    ] {
        let d1 = bridge.curve.derivative_at(t, 1);
        let d2 = bridge.curve.derivative_at(t, 2);
        let tangent = d1.normalized().expect("tangent");
        assert!(tangent.dot(&along) > 1.0 - 1e-12, "t={t}: {tangent:?}");
        // Curvature |d1 x d2| / |d1|^3 (the ellipse quarter's was 5.5).
        let kappa = d1.cross(&d2).length() / d1.length().powi(3);
        assert!(kappa < 1e-4, "t={t}: end curvature {kappa:.3e}");
    }
    for k in 0..=16 {
        let q = bridge.curve.point_at(k as f64 / 16.0);
        assert!((q.x - 1.0).abs() < 1e-12, "off the face: {q:?}");
        assert!(
            q.y < 0.66 + 1e-12 && q.z < 0.66 + 1e-12,
            "beyond the feet: {q:?}"
        );
    }
}

/// A convex cube corner has no bridge: every face corner is 90°.
#[test]
fn a_convex_cube_corner_has_only_corners() {
    let r = 0.1;
    let o = p(0.0, 0.0, 0.0);
    let along = |axis: Vec3, a: Vec3, b: Vec3| {
        strip(
            o + a * r + b * r,
            axis,
            r,
            [(o + a * r, o + a * r + axis), (o + b * r, o + b * r + axis)],
        )
    };
    let (x, y, z) = (
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );
    // Rail order puts each strip's NON-shared-face rail first: the two
    // non-shared rails of adjacent strips are coplanar too (the plane
    // at the contact offset), and must not be read as the shared pair.
    let strips = [along(x, z, y), along(y, x, z), along(z, y, x)];
    for k in 0..3 {
        assert_eq!(
            classify_junction(&strips[k], &strips[(k + 1) % 3], o),
            Ok(JunctionKind::Corner),
            "junction {k}"
        );
    }
}
