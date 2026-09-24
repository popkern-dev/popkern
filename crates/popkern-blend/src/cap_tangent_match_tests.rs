use super::*;
use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;

/// A convex strip end: supports y = 0 and x = 0 (outward normals −y
/// and −x), ball of radius r centred at (r, r, 1), spine leaving +z.
/// `opening` rotates the second contact, so the section arc spans
/// `opening` radians (90° for the square corner).
fn endpoint(r: f64, opening: f64) -> CapEndpoint {
    let center = Point3::new(r, r, 1.0);
    let a = Vec3::new(0.0, -1.0, 0.0);
    let b = Vec3::new(-opening.sin(), -opening.cos(), 0.0);
    CapEndpoint {
        contact_1: center + a * r,
        contact_2: center + b * r,
        normal_1: a,
        normal_2: b,
        center,
        spine_tangent_outward: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    }
}

/// The v = 0 row is the section arc: every sample is on the circle of
/// radius r about the centre, in the plane normal to the spine, and
/// the arc's ends are the two contacts.
#[test]
fn the_section_row_is_the_end_section_arc() {
    for opening in [0.6, std::f64::consts::FRAC_PI_2, 2.4] {
        let ep = endpoint(0.1, opening);
        let s = tangent_match_surface(&ep).expect("cap builds");
        for k in 0..=16 {
            let p = s.point_at(k as f64 / 16.0, 0.0);
            let d = p - ep.center;
            assert!(d.z.abs() < 1e-12, "off the section plane: {d:?}");
            assert!(
                (d.length() - 0.1).abs() < 1e-12,
                "off the circle: {}",
                d.length()
            );
        }
        let ends = [s.point_at(0.0, 0.0), s.point_at(1.0, 0.0)];
        for contact in [ep.contact_1, ep.contact_2] {
            assert!(
                ends.iter().any(|e| (*e - contact).length() < 1e-12),
                "contact {contact:?} is not an arc end"
            );
        }
    }
}

/// Across the section the cap shares the strip's normal, which points
/// away from the axis, and so is G1. Along the meridian, the lateral
/// offset from the ruling shrinks cubically: zero normal curvature,
/// matching the cylinder's straight rulings (G2).
#[test]
fn the_cap_is_curvature_continuous_with_a_cylinder_strip() {
    let ep = endpoint(0.1, std::f64::consts::FRAC_PI_2);
    let s = tangent_match_surface(&ep).expect("cap builds");
    for k in 0..=8 {
        let u = k as f64 / 8.0;
        let p = s.point_at(u, 0.0);
        let radial = (p - ep.center).normalized().expect("off the axis");
        let n = s.normal_at(u, 0.0).normalized().expect("regular normal");
        assert!(
            n.dot(&radial) > 1.0 - 1e-12,
            "u={u}: normal {n:?} vs {radial:?}"
        );
        let lateral = |h: f64| {
            let q = s.point_at(u, h) - p;
            (q - Vec3::new(0.0, 0.0, q.z)).length()
        };
        let ratio = lateral(2e-2) / lateral(1e-2);
        assert!(ratio > 7.0, "u={u}: lateral ratio {ratio} is not cubic");
    }
}

/// The cap closes at a pole r beyond the station, and each meridian is
/// the surface's own iso-curve.
#[test]
fn the_cap_closes_at_its_pole_along_its_meridians() {
    let ep = endpoint(0.1, 2.4);
    let s = tangent_match_surface(&ep).expect("cap builds");
    assert_eq!(s.num_u, 5, "a 2.4 rad arc takes two segments");
    let pole = ep.center + Vec3::new(0.0, 0.0, 0.1);
    for k in 0..=8 {
        assert!((s.point_at(k as f64 / 8.0, 1.0) - pole).length() < 1e-12);
    }
    for (at_end, u) in [(false, 0.0), (true, 1.0)] {
        let m = meridian(&s, at_end);
        for k in 0..=8 {
            let v = k as f64 / 8.0;
            assert!((m.point_at(v) - s.point_at(u, v)).length() < 1e-12);
        }
    }
}

/// Degenerate endpoints decline: coincident contacts, contacts on the
/// axis, zero radius.
#[test]
fn degenerate_endpoints_decline() {
    let mut ep = endpoint(0.1, std::f64::consts::FRAC_PI_2);
    ep.contact_2 = ep.contact_1;
    assert!(tangent_match_surface(&ep).is_none());
    let mut ep = endpoint(0.1, std::f64::consts::FRAC_PI_2);
    ep.contact_1 = ep.center + Vec3::new(0.0, 0.0, 0.05);
    assert!(tangent_match_surface(&ep).is_none());
    let mut ep = endpoint(0.1, std::f64::consts::FRAC_PI_2);
    ep.radius = 0.0;
    assert!(tangent_match_surface(&ep).is_none());
}
