use super::*;
use popkern_core::geometry::curve::CurveEval;

fn constant_helix_spec() -> HelixSpec {
    HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 5.0,
        end_radius: 5.0,
        pitch: HelixPitch::Constant(2.0),
        turns: 2.0,
        start_angle: 0.0,
        left_handed: false,
    }
}

#[test]
fn constant_helix_builds_valid_curve() {
    let spec = constant_helix_spec();
    let (curve, warnings) = build_helix_path(&spec).unwrap(); // ok in tests
    assert!(warnings.is_empty());
    assert!(curve.control_points.len() >= 9);
    let p0 = curve.point_at(0.0);
    let r0 = (p0.x * p0.x + p0.y * p0.y).sqrt();
    assert!((r0 - 5.0).abs() < 0.5, "expected radius ~5.0, got {}", r0);
    assert!(p0.z.abs() < 0.1, "expected z~0 at t=0, got {}", p0.z);
    let p1 = curve.point_at(1.0);
    assert!((p1.z - 4.0).abs() < 0.5, "expected z~4.0, got {}", p1.z);
}

#[test]
fn constant_helix_quarter_turn_check() {
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 1.0,
        end_radius: 1.0,
        pitch: HelixPitch::Constant(4.0),
        turns: 1.0,
        start_angle: 0.0,
        left_handed: false,
    };
    let (curve, _) = build_helix_path(&spec).unwrap(); // ok in tests
    // At t=0.25 (quarter turn), z should be ~1.0 and radius ~1.0.
    let p = curve.point_at(0.25);
    let r = (p.x * p.x + p.y * p.y).sqrt();
    assert!((r - 1.0).abs() < 0.3, "expected radius ~1.0, got {}", r);
    assert!((p.z - 1.0).abs() < 0.5, "expected z~1.0, got {}", p.z);
}

#[test]
fn variable_helix_emits_warning() {
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 1.0,
        end_radius: 3.0,
        pitch: HelixPitch::Constant(2.0),
        turns: 1.0,
        start_angle: 0.0,
        left_handed: false,
    };
    let (curve, warnings) = build_helix_path(&spec).unwrap(); // ok in tests
    assert!(
        warnings
            .iter()
            .any(|w| matches!(w, SweepWarning::HelixApproximated))
    );
    assert!(!curve.control_points.is_empty());
}

fn circle_nurbs(r: f64) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, 0.0),
            Point3::new(r, r, 0.0),
            Point3::new(0.0, r, 0.0),
            Point3::new(-r, r, 0.0),
            Point3::new(-r, 0.0, 0.0),
            Point3::new(-r, -r, 0.0),
            Point3::new(0.0, -r, 0.0),
            Point3::new(r, -r, 0.0),
            Point3::new(r, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

#[test]
fn detect_torus_for_circle_on_circular_helix() {
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 5.0,
        end_radius: 5.0,
        pitch: HelixPitch::Constant(0.0),
        turns: 1.0,
        start_angle: 0.0,
        left_handed: false,
    };
    let profile = SweepProfile::Wire(circle_nurbs(0.5));
    let result = detect_analytic_helix(&spec, &profile);
    assert!(result.is_some(), "expected torus detection");
    match result {
        Some(SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        }) => {
            assert!((major_radius - 5.0).abs() < 1e-6);
            assert!((minor_radius - 0.5).abs() < 0.1);
        }
        _ => panic!("expected Torus variant"),
    }
}

#[test]
fn detect_torus_returns_none_for_nonzero_pitch() {
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 5.0,
        end_radius: 5.0,
        pitch: HelixPitch::Constant(1.0),
        turns: 1.0,
        start_angle: 0.0,
        left_handed: false,
    };
    let profile = SweepProfile::Wire(circle_nurbs(1.0));
    assert!(detect_analytic_helix(&spec, &profile).is_none());
}

/// Radial component (in the axis-z plane) of a point.
fn radial(p: &Point3) -> (f64, f64, f64) {
    ((p.x * p.x + p.y * p.y).sqrt(), p.x, p.y)
}

#[test]
fn phased_helix_is_the_rotated_unphased_helix() {
    let base = constant_helix_spec();
    let phase = std::f64::consts::FRAC_PI_3;
    let phased = HelixSpec {
        start_angle: phase,
        ..base.clone()
    };
    let (c0, w0) = build_helix_path(&base).unwrap(); // ok in tests
    let (c1, w1) = build_helix_path(&phased).unwrap(); // ok in tests
    assert!(w0.is_empty() && w1.is_empty(), "phase must stay exact");
    for i in 0..=16 {
        let t = i as f64 / 16.0;
        let p0 = c0.point_at(t);
        let p1 = c1.point_at(t);
        // Same height, same radius, rotated by exactly the phase.
        assert!((p0.z - p1.z).abs() < 1e-12, "height differs at t={t}");
        let (r0, x0, y0) = radial(&p0);
        let (r1, x1, y1) = radial(&p1);
        assert!((r0 - r1).abs() < 1e-9, "radius differs at t={t}");
        let delta = (y1.atan2(x1) - y0.atan2(x0)).rem_euclid(2.0 * PI);
        assert!(
            (delta - phase).abs() < 1e-9,
            "phase at t={t}: expected {phase}, got {delta}"
        );
    }
}

#[test]
fn left_hand_helix_winds_clockwise_about_the_axis() {
    let right = constant_helix_spec();
    let left = HelixSpec {
        left_handed: true,
        ..right.clone()
    };
    let (cr, _) = build_helix_path(&right).unwrap(); // ok in tests
    let (cl, _) = build_helix_path(&left).unwrap(); // ok in tests
    // Winding sign = z-component of radial(t0) × radial(t1).
    let wind = |c: &NurbsCurve| {
        let p0 = c.point_at(0.02);
        let p1 = c.point_at(0.06);
        p0.x * p1.y - p0.y * p1.x
    };
    assert!(wind(&cr) > 0.0, "right-hand helix must wind CCW");
    assert!(wind(&cl) < 0.0, "left-hand helix must wind CW");
    // Both still advance along +axis.
    assert!(cl.point_at(1.0).z > cl.point_at(0.0).z);
    // Mirror identity: left(t) == right(t) reflected through the
    // u-z plane of the axis frame (exact arcs on both sides).
    for i in 0..=16 {
        let t = i as f64 / 16.0;
        let pr = cr.point_at(t);
        let pl = cl.point_at(t);
        assert!((pr.z - pl.z).abs() < 1e-12);
        let (rr, ..) = radial(&pr);
        let (rl, ..) = radial(&pl);
        assert!((rr - rl).abs() < 1e-9, "radius differs at t={t}");
    }
}

#[test]
fn phased_helix_radius_stays_exact_everywhere() {
    let spec = HelixSpec {
        start_angle: 1.234,
        left_handed: true,
        ..constant_helix_spec()
    };
    let (curve, _) = build_helix_path(&spec).unwrap(); // ok in tests
    for i in 0..=200 {
        let t = i as f64 / 200.0;
        let (r, ..) = radial(&curve.point_at(t));
        assert!(
            (r - 5.0).abs() < 1e-9,
            "rational-arc radius drifted at t={t}: {r}"
        );
    }
}

#[test]
fn invalid_specs_return_errors() {
    let mut spec = constant_helix_spec();
    spec.turns = 0.0;
    assert!(build_helix_path(&spec).is_err());
    spec.turns = 2.0;
    spec.start_radius = 0.0;
    assert!(build_helix_path(&spec).is_err());
    spec.start_radius = 5.0;
    spec.pitch = HelixPitch::Constant(-1.0);
    assert!(build_helix_path(&spec).is_err());
    spec.pitch = HelixPitch::Constant(1.0);
    spec.start_angle = f64::NAN;
    assert!(build_helix_path(&spec).is_err());
}

/// Canonical recognition must not ALIAS on whole-turn periodic
/// paths: at a turn count divisible by 4 the recognizers' former
/// fixed 5-sample path ladder landed every sample one full turn
/// apart, every sampled row was a perfect translate of row 0, and a
/// helical wall passed the extruded-plane translation test — both
/// flanks of a thread-groove tool were minted as Planes and the
/// boolean intersected phantom planes (PK-460). Span-resolving
/// sampling (two samples per knot span) cannot alias.
/// Build an M6-thread-tool-shaped helical sweep (square wire, radius
/// 3, pitch 1, axis-referenced frames — the fastener engine's
/// convention: outward radial authors as −x) and measure the radial
/// ripple of each NURBS wall: sampled ρ = distance-to-axis, per-face
/// min/max/mean. The outer and inner walls are ideal cylinders, so
/// their (max − min) is exactly the skinning error of the swept wall
/// at this turn count.
fn wall_radial_ripple(turns: f64) -> Vec<(f64, f64)> {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    let result = build_square_tool_sweep(turns);
    let mut rows = Vec::new();
    for (_, f) in result.brep.faces.iter() {
        let Some(s) = f.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Nurbs(_) = &s.geometry else {
            continue;
        };
        let (ur, vr) = s.geometry.parameter_range();
        let nv = ((turns * 64.0) as usize).max(256);
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for i in 0..=8 {
            let u = ur.0 + (ur.1 - ur.0) * i as f64 / 8.0;
            for j in 0..=nv {
                let v = vr.0 + (vr.1 - vr.0) * j as f64 / nv as f64;
                let p = s.geometry.point_at(u, v);
                let rho = (p.x * p.x + p.y * p.y).sqrt();
                lo = lo.min(rho);
                hi = hi.max(rho);
            }
        }
        rows.push((lo, hi));
    }
    rows
}

/// The square-wire helical tool the ripple metrics run on (fastener
/// engine conventions: outward radial authors as −x, axis-referenced
/// frames).
fn build_square_tool_sweep(turns: f64) -> crate::types::SweepResult {
    let profile_wire = popkern_core::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![
            Point3::new(-3.1, -0.2, 0.0),
            Point3::new(-2.9, -0.2, 0.0),
            Point3::new(-2.9, 0.2, 0.0),
            Point3::new(-3.1, 0.2, 0.0),
            Point3::new(-3.1, -0.2, 0.0),
        ],
        weights: vec![1.0; 5],
    };
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 3.0,
        end_radius: 3.0,
        pitch: HelixPitch::Constant(1.0),
        turns,
        start_angle: 0.0,
        left_handed: false,
    };
    let axis = popkern_core::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, turns)],
        weights: vec![1.0, 1.0],
    };
    let options = crate::types::SweepOptions {
        frame: crate::types::FrameMethod::Auxiliary(axis),
        ..crate::types::SweepOptions::default()
    };
    crate::helical_sweep(
        &crate::SweepProfile::Wire(profile_wire),
        &spec,
        &options,
        popkern_core::tolerance::Tolerance::default(),
        None,
    )
    .expect("helical sweep")
}

#[test]
#[ignore = "heavy: PK-495 wall-ripple sweep probe — run by hand with PK495_TURNS; the 100-turn bound runs in the default gate"]
fn probe_wall_radial_ripple() {
    let turns: f64 = std::env::var("PK495_TURNS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(20.0);
    for (lo, hi) in wall_radial_ripple(turns) {
        println!(
            "PK495_RIPPLE turns={turns} rho=[{lo:.6},{hi:.6}] span={:.2e}",
            hi - lo
        );
    }
}

/// PK-495 regression: a 100-turn tool wall's skinning ripple stays at
/// the sub-40-turn baseline. Under the old absolute frame cap
/// (`(6·spans).min(1024)`, bit at ~43 turns) plus the fixed 100-cell
/// arc-length table (51-turn placement beat), the ideal-cylinder
/// walls of this sweep rippled 6.8e-2 radially at 100 turns and the
/// radial level-curve SSI arm declined the thread flank pairs
/// (measured, this probe). Baseline at ≤40 turns is 2.4e-3; the
/// bound leaves ~2× headroom without readmitting the defect (a
/// starved wall is an order of magnitude past it).
#[test]
fn hundred_turn_wall_ripple_stays_at_baseline() {
    let worst = wall_radial_ripple(100.0)
        .iter()
        .map(|(lo, hi)| {
            // Ideal-cylinder walls only: their full span IS the
            // skinning error. Flank-analog walls span the profile.
            if hi - lo < 0.1 { hi - lo } else { 0.0 }
        })
        .fold(0.0_f64, f64::max);
    assert!(
        worst < 5.0e-3,
        "100-turn wall ripple {worst:.2e} regressed past the 2.4e-3 baseline"
    );
}

/// PK-495 regression: a 150-turn helical sweep BUILDS — the
/// self-intersection gate's span-tracking grid, once capped at 400
/// samples, stepped 135° of smooth frame rotation per sample on this
/// coil and declined it as a false fold (120 turns passed at 108° by
/// luck). Suspicious steps now refine before convicting, so smooth
/// rotation at any density cannot false-positive.
#[test]
fn long_coil_sweep_is_not_a_false_fold() {
    // Build only — the ripple sampling is the expensive part and the
    // ripple bound is pinned at 100 turns above.
    let result = build_square_tool_sweep(150.0);
    assert!(!result.brep.faces.is_empty(), "150-turn sweep is empty");
}

#[test]
fn four_turn_helical_walls_are_not_recognized_as_planes() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let profile_wire = popkern_core::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![
            Point3::new(-0.1, -0.05, 0.0),
            Point3::new(0.1, -0.05, 0.0),
            Point3::new(0.1, 0.05, 0.0),
            Point3::new(-0.1, 0.05, 0.0),
            Point3::new(-0.1, -0.05, 0.0),
        ],
        weights: vec![1.0; 5],
    };
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 3.0,
        end_radius: 3.0,
        pitch: HelixPitch::Constant(1.0),
        turns: 4.0,
        start_angle: 0.0,
        left_handed: false,
    };
    let result = crate::helical_sweep(
        &crate::SweepProfile::Wire(profile_wire),
        &spec,
        &crate::SweepOptions::default(),
        popkern_core::tolerance::Tolerance::default(),
        None,
    )
    .expect("4-turn helical sweep");
    let planes = result
        .brep
        .faces
        .iter()
        .filter(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }))
        })
        .count();
    // The two end caps are legitimately planar; the four helical
    // walls must stay NURBS.
    assert_eq!(
        planes, 2,
        "{planes} planar faces on a 4-turn helical sweep (walls aliased into planes)"
    );
}
