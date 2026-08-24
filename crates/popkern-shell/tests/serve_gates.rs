//! Serve-gate and decline-path tests for popkern-shell (PK-321/PK-330):
//! collapse declines with binding entities, curved-body validation
//! declines (un-narrowed by PK-324), open-sheet containment, and
//! approximation-honesty warnings.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn shell_past_feature_size_declines_offset_too_large() {
    // A 0.6 wall on a unit cube crosses the opposite walls — a pure
    // plane-plane collision the retired midpoint-normal detector was
    // structurally blind to. The collapse gate declines it with the
    // binding edge and a usable max-safe estimate.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        0.6,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    match result {
        Err(ShellError::OffsetTooLarge {
            requested,
            max_safe,
            binding,
        }) => {
            assert!((requested - 0.6).abs() < 1e-12);
            assert!(max_safe > 0.0 && max_safe < 0.6, "max_safe {max_safe}");
            assert!(binding.contains("edge"), "{binding}");
        }
        other => panic!("expected OffsetTooLarge decline, got {other:?}"),
    }
}

#[test]
fn offset_shrink_past_feature_size_declines() {
    // Shrinking a unit cube by 0.6 per side inverts it — the collapse
    // gate must decline rather than serve an inside-out body.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = popkern_shell::offset::offset_solid(&brep, -0.6, &ShellOptions::default(), tol());
    assert!(
        matches!(result, Err(ShellError::OffsetTooLarge { .. })),
        "expected OffsetTooLarge decline, got {result:?}"
    );
}

#[test]
fn shell_cylinder_declines_until_curved_recompute() {
    // Curved-body shells are geometrically wrong today: edge recompute
    // keeps the un-offset clone curves for non-planar surface pairs, so
    // the assembled body fails structural validation. The serve gate
    // turns that from silent wrong geometry into a diagnosed decline.
    // PK-324 (general SSI edge/vertex recompute) un-narrows this — when
    // it lands, this pin flips to a validate_all-clean serving test.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    match result {
        Err(ShellError::ValidationFailed { errors }) => {
            assert!(!errors.is_empty());
        }
        Ok(_) => panic!(
            "cylinder shell served — if curved recompute landed (PK-324), \
             upgrade this test to assert validate_all-clean serving"
        ),
        Err(other) => panic!("expected ValidationFailed decline, got {other:?}"),
    }
}

#[test]
fn shell_thin_prism_declines_offset_too_large() {
    // A 0.2 wall in a prism only 0.05 thick in Z cannot exist: the skin
    // inverts. The collapse gate declines with the binding edge (the old
    // contract accepted "any valid result or error"; the honest outcome
    // is exactly one thing).
    use popkern_core::geometry::{Point3, Vec3};

    let profile = &[
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let brep = BRep::prism(profile, Vec3::new(0.0, 0.0, 1.0), 0.05, tol()).unwrap();
    let result = shell(
        &brep,
        0.2,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    assert!(
        matches!(result, Err(ShellError::OffsetTooLarge { .. })),
        "expected OffsetTooLarge decline, got {result:?}"
    );
}

#[test]
fn shell_midplane_past_feature_size_declines() {
    // 1.2 total wall on a unit cube: the inward half (0.6) inverts.
    // The decline reports the caller's full-thickness units.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        1.2,
        ShellDirection::Midplane,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    match result {
        Err(ShellError::OffsetTooLarge {
            requested,
            max_safe,
            ..
        }) => {
            assert!((requested - 1.2).abs() < 1e-9, "requested {requested}");
            assert!(max_safe > 0.9 && max_safe < 1.2, "max_safe {max_safe}");
        }
        other => panic!("expected OffsetTooLarge decline, got {other:?}"),
    }
}

#[test]
fn thicken_open_sheet_declines_open_boundary() {
    // Manufacture an open sheet: delete one face and its loop from a
    // cube, leaving its half-edges with dangling loop references — the
    // exact free-rim shape thicken's precondition detects. The old
    // pipeline served this as two disconnected skins with an unclosed
    // rim while claiming to build side walls.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    let outer_loop = brep.faces[fh].outer_loop;
    brep.loops.remove(outer_loop);
    brep.faces.remove(fh);

    let result = popkern_shell::thicken::thicken(&brep, 0.05, &ShellOptions::default(), tol());
    match result {
        Err(ShellError::OpenBoundary { free_half_edges }) => {
            assert_eq!(
                free_half_edges, 4,
                "one deleted cube face frees 4 half-edges"
            );
        }
        other => panic!("expected OpenBoundary decline, got {other:?}"),
    }
}

#[test]
fn nurbs_offset_path_emits_approximation_warning() {
    // A cone face has no closed-form offset yet (PK-322), so its offset
    // goes through the NURBS sampler. The dispatch must SAY so: every
    // NURBS-path face carries FaceOffsetApproximated with the measured
    // deviation (the warning existed for years and was never emitted).
    use popkern_shell::classify::classify_faces;
    use popkern_shell::surface_offset::compute_offset_surfaces;

    let brep = BRep::cone(1.0, 2.0, tol()).unwrap();
    // 0.01 sits inside the cone's sampled near-apex curvature radius
    // (~0.037), so the veto stays quiet and the NURBS path runs.
    let classification = classify_faces(&brep, 0.01, &[], ShellDirection::Inward, &[]).unwrap();
    let result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
    let approx_warnings: Vec<_> = result
        .warnings
        .iter()
        .filter(|w| matches!(w, ShellWarning::FaceOffsetApproximated { .. }))
        .collect();
    assert!(
        !approx_warnings.is_empty() || !result.failed.is_empty(),
        "the cone lateral face must either warn (approximated within \
         tolerance) or fail out loud (out of tolerance) — never silence: \
         warnings {:?}, failed {:?}",
        result.warnings,
        result.failed
    );
}
