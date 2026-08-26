//! Entry-point tests for [`crate::shell`]/[`crate::shell_variable`]
//! (split from `lib.rs`, 500-line ratchet).

use crate::*;
use popkern_core::tolerance::Tolerance;

#[test]
fn shell_empty_brep_returns_no_faces() {
    let brep = BRep::new(Tolerance::default());
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        Tolerance::default(),
        None,
    );
    assert!(matches!(result, Err(ShellError::NoFaces)));
}

#[test]
fn shell_cube_no_removal_produces_result() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap();
    // Outer (6) + inner (6) = 12 faces.
    assert_eq!(result.brep.faces.len(), 12);
}

#[test]
fn shell_cube_one_face_removed() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[first_face],
        &[],
        &ShellOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap();
    // Outer (5) + inner (5) + side walls (4) = 14 faces.
    assert_eq!(result.brep.faces.len(), 14);
}

#[test]
fn shell_variable_produces_result() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let first_face = brep.faces.keys().next().unwrap();
    let result = shell_variable(
        &brep,
        0.1,
        &[(first_face, 0.2)],
        ShellDirection::Inward,
        &[],
        &ShellOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap();
    assert_eq!(result.brep.faces.len(), 12);
}
