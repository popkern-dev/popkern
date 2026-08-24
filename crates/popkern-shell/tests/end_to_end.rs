//! End-to-end integration tests for popkern-shell.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Assert the Euler characteristic V - E + F = 2 for a simple solid (one shell,
/// genus 0, no inner loops beyond faces).
fn assert_euler_solid(brep: &BRep) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    assert_eq!(
        v - e + f,
        2,
        "Euler failed for solid: V={v} E={e} F={f}, V-E+F={}",
        v - e + f
    );
}

/// Assert the Euler characteristic for a compound body (e.g., shell or thicken
/// result with outer + inner skins). The expected Euler characteristic is
/// V - E + F = 2 * C where C is the number of closed connected components.
/// For a shell/thicken with no removed faces, C = 2 (outer + inner).
fn assert_euler_compound(brep: &BRep, expected_components: i64) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    let expected = 2 * expected_components;
    assert_eq!(
        v - e + f,
        expected,
        "Euler failed for compound body: V={v} E={e} F={f}, \
         V-E+F={} expected 2*{expected_components}={expected}",
        v - e + f
    );
}

/// Assert the body passes the FULL validator with zero errors. Served
/// shell results carry a structural serve gate in the pipeline (PK-321),
/// so tests hold them to the stronger validate_all bar; the historical
/// catch_unwind shim died with the stale-handle defect it tolerated.
fn assert_structurally_valid(brep: &BRep) {
    let errors = brep.validate_all();
    assert!(errors.is_empty(), "validation errors: {errors:?}");
}
#[test]
fn shell_cube_no_removal() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Outer (6) + inner (6) = 12 faces, two disjoint closed surfaces.
    assert_eq!(result.brep.faces.len(), 12);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_cube_one_face_removed() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[first_face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Outer (5) + inner (5) + side walls (4) = 14 faces.
    assert_eq!(result.brep.faces.len(), 14);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_cube_two_adjacent_removed() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces: Vec<_> = brep.faces.keys().collect();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[faces[0], faces[1]],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    // Outer (4) + inner (4) + side walls around the merged opening.
    // Two adjacent removed faces open one 6-edge rim: 4 + 4 + 6 = 14.
    assert_eq!(result.brep.faces.len(), 14);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_variable_thickness() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = shell_variable(
        &brep,
        0.1,
        &[(first_face, 0.2)],
        ShellDirection::Inward,
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 12);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_inward_outward_midplane() {
    let brep = BRep::unit_cube(tol()).unwrap();

    let r_in = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    assert_eq!(r_in.brep.faces.len(), 12);
    assert_euler_compound(&r_in.brep, 2);
    assert_structurally_valid(&r_in.brep);

    let r_out = shell(
        &brep,
        0.1,
        ShellDirection::Outward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    assert_eq!(r_out.brep.faces.len(), 12);
    assert_euler_compound(&r_out.brep, 2);
    assert_structurally_valid(&r_out.brep);

    // Midplane is honest containment now: the pipeline cannot move the
    // outer skin, so serving would deliver an inward half-shell under a
    // midplane label. The two-skin construction is its own filed task.
    let r_mid = shell(
        &brep,
        0.1,
        ShellDirection::Midplane,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    );
    assert!(
        matches!(r_mid, Err(ShellError::MidplaneUnsupported)),
        "expected MidplaneUnsupported decline, got {r_mid:?}"
    );
}

#[test]
fn offset_solid_outward() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::offset::offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
    // Topology counts must be preserved: 6 faces, 12 edges, 8 vertices.
    assert_eq!(result.brep.faces.len(), 6);
    assert_eq!(result.brep.edges.len(), 12);
    assert_eq!(result.brep.vertices.len(), 8);
    assert_euler_solid(&result.brep);
    assert_structurally_valid(&result.brep);
    // Every vertex must have moved from the original unit cube positions.
    // Original vertices have coordinates in {0.0, 1.0}; after offset no
    // coordinate should remain exactly at 0.0 or 1.0.
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - 0.0).abs() > 1e-6 && (c - 1.0).abs() > 1e-6,
                "vertex coord {c} unchanged after offset"
            );
        }
    }
}

#[test]
fn offset_solid_negative() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::offset::offset_solid(&brep, -0.1, &ShellOptions::default(), tol()).unwrap();
    // Topology counts must be preserved: 6 faces, 12 edges, 8 vertices.
    assert_eq!(result.brep.faces.len(), 6);
    assert_eq!(result.brep.edges.len(), 12);
    assert_eq!(result.brep.vertices.len(), 8);
    assert_euler_solid(&result.brep);
    assert_structurally_valid(&result.brep);
    // Every vertex must have moved from the original unit cube positions.
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - 0.0).abs() > 1e-6 && (c - 1.0).abs() > 1e-6,
                "vertex coord {c} unchanged after negative offset"
            );
        }
    }
}

#[test]
fn thicken_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::thicken::thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    // Outer (6) + inner (6) = 12 faces (no side walls for closed solid).
    assert_eq!(result.brep.faces.len(), 12);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);
}

#[test]
fn shell_empty_brep_error() {
    let brep = BRep::new(tol());
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
    assert!(matches!(result, Err(ShellError::NoFaces)));
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

// ---- New tests below ----

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
fn offset_solid_variable_per_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = popkern_shell::offset::offset_solid_variable(
        &brep,
        0.1,
        &[(first_face, 0.3)],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 6);
    assert_euler_solid(&result.brep);
    assert_structurally_valid(&result.brep);
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
fn thicken_prism_slab() {
    // Build a thin rectangular prism (acts as a slab) and thicken it.
    // The prism has 6 faces. Thicken produces outer + inner = 12 faces.
    use popkern_core::geometry::{Point3, Vec3};

    let profile = &[
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let brep = BRep::prism(profile, Vec3::new(0.0, 0.0, 1.0), 0.05, tol()).unwrap();
    // A 0.02 wall fits inside the 0.05-thick slab (0.01 inward per side).
    let result =
        popkern_shell::thicken::thicken(&brep, 0.02, &ShellOptions::default(), tol()).unwrap();
    // Thicken a closed solid: outer (6) + inner (6) = 12 faces.
    assert_eq!(result.brep.faces.len(), 12);
    assert_euler_compound(&result.brep, 2);
    assert_structurally_valid(&result.brep);

    // A 0.2 wall does NOT fit: the inward half-offset (0.1) inverts the
    // 0.05-thick slab. This used to serve an inverted "hollow" body.
    let too_thick = popkern_shell::thicken::thicken(&brep, 0.2, &ShellOptions::default(), tol());
    assert!(
        matches!(too_thick, Err(ShellError::OffsetTooLarge { .. })),
        "expected OffsetTooLarge decline, got {too_thick:?}"
    );
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
