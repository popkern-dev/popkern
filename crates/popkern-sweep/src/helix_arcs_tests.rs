use super::*;
use crate::types::HelixPitch;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_math::PortableMath;

fn spec(pitch: f64, radius: f64, turns: f64, left_handed: bool) -> HelixSpec {
    HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: radius,
        end_radius: radius,
        pitch: HelixPitch::Constant(pitch),
        turns,
        start_angle: 0.3,
        left_handed,
    }
}

/// The largest radius error and the largest rise error (z against the
/// unwrapped angle's share of the pitch) over a dense sampling.
fn helix_errors(s: &HelixSpec, pitch: f64) -> (f64, f64) {
    let curve = match build_constant_helix(s, pitch) {
        Ok(c) => c,
        Err(e) => panic!("helix build failed: {e}"),
    };
    let winding = if s.left_handed { -1.0 } else { 1.0 };
    let (mut radius_err, mut rise_err) = (0.0f64, 0.0f64);
    let mut unwrapped = 0.0;
    // The axis frame is the builder's, so unwrap from the first sample.
    let first = curve.point_at(0.0);
    let mut previous = first.y.patan2(first.x);
    let samples = 4000;
    for k in 0..=samples {
        let p = curve.point_at(k as f64 / samples as f64);
        let angle = p.y.patan2(p.x);
        let mut delta = angle - previous;
        while delta > std::f64::consts::PI {
            delta -= 2.0 * std::f64::consts::PI;
        }
        while delta < -std::f64::consts::PI {
            delta += 2.0 * std::f64::consts::PI;
        }
        unwrapped += delta;
        previous = angle;
        radius_err = radius_err.max((p.x.phypot(p.y) - s.start_radius).abs());
        let ideal = winding * unwrapped / (2.0 * std::f64::consts::PI) * pitch;
        rise_err = rise_err.max((p.z - ideal).abs());
    }
    (radius_err, rise_err)
}

#[test]
fn rise_tracks_the_angle_on_thread_pitches() {
    // M3, M6 and M10 coarse pitches on their major radii; both hands,
    // whole and fractional turns.
    for &(pitch, radius) in &[(0.5, 1.5), (1.0, 3.0), (1.5, 5.0)] {
        for &turns in &[1.0, 3.25, 12.0] {
            for &left in &[false, true] {
                let (r, z) = helix_errors(&spec(pitch, radius, turns, left), pitch);
                assert!(r < 1e-12 * radius.max(1.0) * 10.0, "radius error {r}");
                assert!(
                    z < 1e-6 * pitch,
                    "rise error {z} at pitch {pitch} turns {turns} left {left}"
                );
            }
        }
    }
}

#[test]
fn arcs_are_degree_five_with_shared_joints() {
    let s = spec(0.5, 1.5, 2.0, false);
    let curve = match build_constant_helix(&s, 0.5) {
        Ok(c) => c,
        Err(e) => panic!("{e}"),
    };
    let arcs = 2 * ARCS_PER_TURN;
    assert_eq!(curve.degree, HELIX_DEGREE);
    assert_eq!(curve.control_points.len(), arcs * HELIX_DEGREE + 1);
    assert_eq!(
        curve.knots.len(),
        curve.control_points.len() + HELIX_DEGREE + 1
    );
    // Endpoints are exact.
    let end = curve.point_at(1.0);
    assert!((end.z - 1.0).abs() < 1e-14);
    let start = curve.point_at(0.0);
    assert!((start.x.phypot(start.y) - 1.5).abs() < 1e-14 && start.z.abs() < 1e-14);
}

#[test]
fn joints_are_tangent_continuous() {
    let s = spec(1.0, 3.0, 1.0, false);
    let curve = match build_constant_helix(&s, 1.0) {
        Ok(c) => c,
        Err(e) => panic!("{e}"),
    };
    for a in 1..ARCS_PER_TURN {
        let u = a as f64 / ARCS_PER_TURN as f64;
        let before = curve.point_at(u) - curve.point_at(u - 1e-7);
        let after = curve.point_at(u + 1e-7) - curve.point_at(u);
        let (Some(b), Some(f)) = (before.normalized(), after.normalized()) else {
            panic!("degenerate tangent at joint {a}");
        };
        assert!(b.dot(&f) > 1.0 - 1e-10, "joint {a} kinks: {}", b.dot(&f));
    }
}
