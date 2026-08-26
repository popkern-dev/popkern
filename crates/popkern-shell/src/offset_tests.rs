use super::*;
use popkern_core::tolerance::Tolerance;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------
// offset_solid — existing tests
// ---------------------------------------------------------------

#[test]
fn offset_cube_outward() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
    // Minkowski-faithful default (PK-325): the rounded box.
    assert_eq!(result.brep.faces.len(), 26);
    // Every corner-sphere vertex sits at an original corner
    // displaced by 0.1 along one face normal: coordinates come
    // from {-0.1, 0, 1, 1.1}.
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                [-0.1, 0.0, 1.0, 1.1]
                    .iter()
                    .any(|expect| (c - expect).abs() < 1e-9),
                "vertex coord {c} unexpected after outward arc offset"
            );
        }
    }
}

#[test]
fn offset_cube_negative() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid(&brep, -0.1, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
    // Negative distance shrinks the solid inward: corners move from
    // {0,1} to {0.1, 0.9}.
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - 0.1).abs() < 0.05 || (c - 0.9).abs() < 0.05,
                "vertex coord {c} unexpected after inward offset"
            );
        }
    }
}

#[test]
fn offset_empty_brep_errors() {
    let brep = BRep::new(tol());
    let result = offset_solid(&brep, 0.1, &ShellOptions::default(), tol());
    assert!(matches!(result, Err(ShellError::NoFaces)));
}

#[test]
fn offset_variable_produces_result() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = offset_solid_variable(
        &brep,
        0.1,
        &[(first_face, 0.2)],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

// ---------------------------------------------------------------
// offset_solid_variable — dedicated tests
// ---------------------------------------------------------------

#[test]
fn variable_empty_brep_errors() {
    let brep = BRep::new(tol());
    let result = offset_solid_variable(&brep, 0.1, &[], &ShellOptions::default(), tol());
    assert!(
        matches!(result, Err(ShellError::NoFaces)),
        "variable offset on empty BRep should return NoFaces"
    );
}

#[test]
fn variable_with_no_overrides_matches_uniform() {
    // offset_solid_variable with empty face_distances should behave
    // identically to offset_solid with the same default distance.
    let brep = BRep::unit_cube(tol()).unwrap();
    let uniform = offset_solid(&brep, 0.15, &ShellOptions::default(), tol()).unwrap();
    let variable =
        offset_solid_variable(&brep, 0.15, &[], &ShellOptions::default(), tol()).unwrap();
    assert_eq!(uniform.brep.faces.len(), variable.brep.faces.len());
    assert_eq!(uniform.brep.vertices.len(), variable.brep.vertices.len());
}

#[test]
fn variable_negative_default_shrinks_inward() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid_variable(&brep, -0.1, &[], &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
    // Negative default should shrink: corners at {0.1, 0.9}.
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - 0.1).abs() < 0.05 || (c - 0.9).abs() < 0.05,
                "vertex coord {c} unexpected after negative variable offset"
            );
        }
    }
}

#[test]
fn variable_per_face_override_differs_from_default() {
    // When one face gets a larger offset, the result's vertex positions
    // should differ from a uniform offset.
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();

    let uniform = offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
    let variable = offset_solid_variable(
        &brep,
        0.1,
        &[(first_face, 0.3)],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();

    // The uniform run serves the rounded box (arc default); the
    // variable run has no uniform rolling ball and contains to the
    // sharp join (E4 downgrade) with its 6 faces preserved.
    assert_eq!(uniform.brep.faces.len(), 26);
    assert_eq!(variable.brep.faces.len(), 6);

    // At least some vertices should differ between uniform and variable
    // offset results. Collect all vertex positions as sorted tuples for
    // comparison.
    let positions = |brep: &BRep| -> Vec<[i64; 3]> {
        brep.vertices
            .iter()
            .map(|(_, v)| {
                [
                    (v.position.x * 1000.0).round() as i64,
                    (v.position.y * 1000.0).round() as i64,
                    (v.position.z * 1000.0).round() as i64,
                ]
            })
            .collect()
    };
    let uni_pos = positions(&uniform.brep);
    let var_pos = positions(&variable.brep);
    assert_ne!(
        uni_pos, var_pos,
        "variable offset with a per-face override should produce \
         different vertex positions than uniform offset"
    );
}

#[test]
fn variable_multiple_overrides() {
    // Provide per-face overrides for two faces.
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let overrides = vec![(faces[0], 0.2), (faces[1], 0.05)];
    let result =
        offset_solid_variable(&brep, 0.1, &overrides, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn variable_zero_default_with_overrides() {
    // Default distance is 0 — only faces with explicit overrides move.
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = offset_solid_variable(
        &brep,
        0.0,
        &[(first_face, 0.15)],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn variable_preserves_face_count() {
    // No matter what per-face distances we choose (within reason),
    // the face count should stay at 6 for a cube.
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let overrides: Vec<(FaceHandle, f64)> = faces
        .iter()
        .enumerate()
        .map(|(i, &fh)| (fh, 0.05 + 0.02 * i as f64))
        .collect();
    let result =
        offset_solid_variable(&brep, 0.1, &overrides, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn variable_negative_override_with_positive_default() {
    // face_distances uses .abs(), so a negative override is treated
    // as its absolute value. Verify it doesn't error.
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = offset_solid_variable(
        &brep,
        0.1,
        &[(first_face, -0.2)],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn variable_allow_partial_false_still_succeeds_on_cube() {
    // A cube with small offsets should succeed even with
    // allow_partial = false.
    let brep = BRep::unit_cube(tol()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let opts = ShellOptions {
        allow_partial: false,
        ..ShellOptions::default()
    };
    let result = offset_solid_variable(&brep, 0.1, &[(first_face, 0.2)], &opts, tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn variable_small_offset_vertices_stay_close() {
    // A very small offset should produce vertices close to the original.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid_variable(&brep, 0.01, &[], &ShellOptions::default(), tol()).unwrap();
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - (-0.01)).abs() < 0.05
                    || (c - 0.0).abs() < 0.05
                    || (c - 1.0).abs() < 0.05
                    || (c - 1.01).abs() < 0.05,
                "vertex coord {c} too far from expected range \
                 after small outward offset"
            );
        }
    }
}

#[test]
fn variable_result_has_no_warnings_on_clean_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid_variable(&brep, 0.1, &[], &ShellOptions::default(), tol()).unwrap();
    assert!(
        result.warnings.is_empty(),
        "a clean cube offset should produce no warnings, got {} warnings",
        result.warnings.len()
    );
}
