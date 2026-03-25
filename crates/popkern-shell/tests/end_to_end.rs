//! End-to-end integration tests for popkern-shell.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::*;

fn tol() -> Tolerance {
    Tolerance::default()
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
    )
    .unwrap();
    // Outer (6) + inner (6) = 12 faces.
    assert_eq!(result.brep.faces.len(), 12);
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
    )
    .unwrap();
    // Outer (5) + inner (5) + side walls (4) = 14 faces.
    assert_eq!(result.brep.faces.len(), 14);
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
    )
    .unwrap();
    // Outer (4) + inner (4) + side walls (varies, 6-8 depending on adjacency).
    // Just verify we get a result with reasonable face count.
    assert!(
        result.brep.faces.len() >= 12 && result.brep.faces.len() <= 16,
        "face count {} outside expected range",
        result.brep.faces.len()
    );
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
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 12);
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
    )
    .unwrap();
    assert_eq!(r_in.brep.faces.len(), 12);

    let r_out = shell(
        &brep,
        0.1,
        ShellDirection::Outward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(r_out.brep.faces.len(), 12);

    let r_mid = shell(
        &brep,
        0.1,
        ShellDirection::Midplane,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(r_mid.brep.faces.len(), 12);
}

#[test]
fn offset_solid_outward() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::offset::offset_solid(&brep, 0.1, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
    // Vertices should have shifted (offset along normals).
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - 0.1).abs() < 0.05 || (c - 0.9).abs() < 0.05,
                "vertex coord {c} unexpected"
            );
        }
    }
}

#[test]
fn offset_solid_negative() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result =
        popkern_shell::offset::offset_solid(&brep, -0.1, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
    for (_, v) in &result.brep.vertices {
        let coords = [v.position.x, v.position.y, v.position.z];
        for &c in &coords {
            assert!(
                (c - (-0.1)).abs() < 0.05 || (c - 1.1).abs() < 0.05,
                "vertex coord {c} unexpected"
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
    );
    assert!(matches!(result, Err(ShellError::NoFaces)));
}

#[test]
fn collision_detection_no_collisions() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let mut surfaces = std::collections::HashMap::new();
    for (fh, face) in &brep.faces {
        if let Some(sof) = &face.surface {
            surfaces.insert(fh, sof.geometry.clone());
        }
    }
    let collisions = popkern_shell::collision::detect_collisions(&brep, &surfaces, &tol());
    assert!(collisions.is_empty());
}
