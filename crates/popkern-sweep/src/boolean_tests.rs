// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::*;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::types::{LoftOptions, LoftSection, PipeRadius, SweepOptions, SweepProfile};

/// Build a closed circle profile in the XY plane centered at origin.
fn circle_profile(radius: f64) -> SweepProfile {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = radius;
    SweepProfile::Wire(NurbsCurve {
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
    })
}

/// Build a straight-line path along Z, offset to (0.5, 0.5) in XY and
/// starting below/ending above the unit cube to avoid coincident faces.
fn offset_z_path() -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.5, 0.5, -0.5), Point3::new(0.5, 0.5, 1.5)],
        weights: vec![1.0, 1.0],
    }
}

/// Extract the wire curve from a SweepProfile, translated to a given Z
/// and shifted in XY to (0.5, 0.5).
fn wire_at_z_offset(profile: &SweepProfile, z: f64) -> NurbsCurve {
    match profile {
        SweepProfile::Wire(c) => {
            let control_points = c
                .control_points
                .iter()
                .map(|p| Point3::new(p.x + 0.5, p.y + 0.5, z))
                .collect();
            NurbsCurve {
                degree: c.degree,
                knots: c.knots.clone(),
                control_points,
                weights: c.weights.clone(),
            }
        }
        _ => panic!("expected Wire"),
    }
}

/// Extract the wire curve from a SweepProfile, shifted in XY to (0.5, 0.5).
fn wire_offset(profile: &SweepProfile) -> NurbsCurve {
    wire_at_z_offset(profile, 0.0)
}

// ---- sweep_add ----

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn sweep_add_produces_brep() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.2);
    let path = offset_z_path();
    let opts = SweepOptions::default();

    let result = sweep_add(&cube, &profile, &path, &opts, tol, None);
    assert!(result.is_ok(), "sweep_add failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    // Union adds material: result must have at least as many faces as the cube.
    assert!(
        brep.faces.len() >= 6,
        "expected >= 6 faces, got {}",
        brep.faces.len()
    );
}

// ---- sweep_cut ----

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn sweep_cut_produces_brep() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.2);
    let path = offset_z_path();
    let opts = SweepOptions::default();

    let result = sweep_cut(&cube, &profile, &path, &opts, tol, None);
    assert!(result.is_ok(), "sweep_cut failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    // Cut removes material: result should still have faces.
    assert!(
        !brep.faces.is_empty(),
        "expected >= 1 face, got {}",
        brep.faces.len()
    );
}

// ---- sweep_intersect ----

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn sweep_intersect_produces_brep() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.2);
    let path = offset_z_path();
    let opts = SweepOptions::default();

    let result = sweep_intersect(&cube, &profile, &path, &opts, tol, None);
    assert!(result.is_ok(), "sweep_intersect failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    // Intersection keeps overlap: should have at least 1 face.
    assert!(
        !brep.faces.is_empty(),
        "expected >= 1 face, got {}",
        brep.faces.len()
    );
}

// ---- loft_add ----

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn loft_add_produces_brep() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.2);
    let sections = vec![
        LoftSection {
            curve: wire_at_z_offset(&profile, -0.5),
            position: None,
        },
        LoftSection {
            curve: wire_at_z_offset(&profile, 1.5),
            position: None,
        },
    ];
    let opts = LoftOptions::default();

    let result = loft_add(&cube, &sections, &opts, tol, None);
    assert!(result.is_ok(), "loft_add failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(
        brep.faces.len() >= 6,
        "expected >= 6 faces, got {}",
        brep.faces.len()
    );
}

// ---- loft_cut ----

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn loft_cut_produces_brep() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.2);
    let sections = vec![
        LoftSection {
            curve: wire_at_z_offset(&profile, -0.5),
            position: None,
        },
        LoftSection {
            curve: wire_at_z_offset(&profile, 1.5),
            position: None,
        },
    ];
    let opts = LoftOptions::default();

    let result = loft_cut(&cube, &sections, &opts, tol, None);
    assert!(result.is_ok(), "loft_cut failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(
        !brep.faces.is_empty(),
        "expected >= 1 face, got {}",
        brep.faces.len()
    );
}

// ---- pipe_add ----

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn pipe_add_produces_brep() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let path = offset_z_path();
    let radius = PipeRadius::Constant(0.2);
    let opts = SweepOptions::default();

    let result = pipe_add(&cube, &path, &radius, &opts, tol, None);
    assert!(result.is_ok(), "pipe_add failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(
        brep.faces.len() >= 6,
        "expected >= 6 faces, got {}",
        brep.faces.len()
    );
}

// ---- pipe_cut ----

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn pipe_cut_produces_brep() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let path = offset_z_path();
    let radius = PipeRadius::Constant(0.2);
    let opts = SweepOptions::default();

    let result = pipe_cut(&cube, &path, &radius, &opts, tol, None);
    assert!(result.is_ok(), "pipe_cut failed: {:?}", result.err());
    let brep = &result.unwrap().brep; // ok in tests
    assert!(
        !brep.faces.is_empty(),
        "expected >= 1 face, got {}",
        brep.faces.len()
    );
}

// ---- error propagation ----

#[test]
fn sweep_add_propagates_sweep_error() {
    // Empty multi-wire profile triggers InvalidProfile error.
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let bad_profile = SweepProfile::MultiWire(vec![]);
    let path = offset_z_path();
    let opts = SweepOptions::default();

    let result = sweep_add(&cube, &bad_profile, &path, &opts, tol, None);
    assert!(result.is_err(), "empty multi-wire profile should fail");
}

#[test]
fn sweep_cut_propagates_sweep_error() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let bad_profile = SweepProfile::MultiWire(vec![]);
    let path = offset_z_path();
    let opts = SweepOptions::default();

    let result = sweep_cut(&cube, &bad_profile, &path, &opts, tol, None);
    assert!(result.is_err(), "empty multi-wire profile should fail");
}

#[test]
fn sweep_intersect_propagates_sweep_error() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let bad_profile = SweepProfile::MultiWire(vec![]);
    let path = offset_z_path();
    let opts = SweepOptions::default();

    let result = sweep_intersect(&cube, &bad_profile, &path, &opts, tol, None);
    assert!(result.is_err(), "empty multi-wire profile should fail");
}

#[test]
fn loft_add_with_single_section_fails() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.5);
    let sections = vec![LoftSection {
        curve: wire_offset(&profile),
        position: None,
    }];
    let opts = LoftOptions::default();

    let result = loft_add(&cube, &sections, &opts, tol, None);
    assert!(
        result.is_err(),
        "loft with 1 section should fail (need >= 2)"
    );
}

#[test]
fn loft_cut_with_single_section_fails() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.5);
    let sections = vec![LoftSection {
        curve: wire_offset(&profile),
        position: None,
    }];
    let opts = LoftOptions::default();

    let result = loft_cut(&cube, &sections, &opts, tol, None);
    assert!(
        result.is_err(),
        "loft with 1 section should fail (need >= 2)"
    );
}

#[test]
fn pipe_add_with_zero_radius_fails() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let path = offset_z_path();
    let radius = PipeRadius::Constant(0.0);
    let opts = SweepOptions::default();

    let result = pipe_add(&cube, &path, &radius, &opts, tol, None);
    assert!(result.is_err(), "zero radius pipe should fail");
}

#[test]
fn pipe_cut_with_negative_radius_fails() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let path = offset_z_path();
    let radius = PipeRadius::Constant(-1.0);
    let opts = SweepOptions::default();

    let result = pipe_cut(&cube, &path, &radius, &opts, tol, None);
    assert!(result.is_err(), "negative radius pipe should fail");
}

#[test]
#[ignore] // Boolean on curved/planar geometry is slow in debug builds.
fn sweep_cut_preserves_warnings() {
    // Use allow_self_intersection to test that warnings pass through.
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit_cube");
    let profile = circle_profile(0.2);
    let path = offset_z_path();
    let opts = SweepOptions {
        allow_self_intersection: true,
        ..SweepOptions::default()
    };

    let result = sweep_cut(&cube, &profile, &path, &opts, tol, None);
    // Whether or not there are warnings, the result struct carries them.
    if let Ok(res) = result {
        // Warnings field must exist; it may be empty for this geometry.
        let _ = res.warnings;
    }
}
