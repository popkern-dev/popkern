//! PK-306: hollow shell/thicken outputs carry honest shell records and
//! pass structural validation.
//!
//! Root cause fixed here: `assemble_shell` copied BOTH skins into one
//! shell record, so every hollow body (two disconnected closed skins,
//! χ = 4) failed the shell-aware Euler–Poincaré check ("expected 2,
//! got 4") — hollow `shell()` output had NEVER passed `validate()`.
//! The void skin is now its own shell of the solid; the side-wall case
//! (skins connected through the openings) stays one shell.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::thicken::thicken;
use popkern_shell::{ShellDirection, ShellOptions, shell};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Hollow cube shell: two shells recorded, validates clean.
#[test]
fn hollow_cube_shell_records_two_shells_and_validates() {
    let cube = BRep::unit_cube(tol()).expect("cube builds");
    let result = shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("hollow shell serves");
    assert_eq!(
        result.brep.shells.len(),
        2,
        "hollow body = outer skin + void skin, one shell each"
    );
    assert!(
        result.brep.validate().is_ok(),
        "hollow shell must validate: {:?}",
        result.brep.validate().err()
    );
}

/// Open shell (one face removed → side walls connect the skins):
/// ONE shell, edge-manifold, validating clean. Before PK-307 the wall
/// builder fabricated all four wall edges, leaving the rim half-edges
/// loopless (`BrokenCycle`) and duplicating every vertical.
#[test]
fn open_shell_is_one_connected_manifold_shell() {
    let cube = BRep::unit_cube(tol()).expect("cube builds");
    let first_face = cube.faces.keys().next().expect("cube has faces");
    let result = shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[first_face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("open shell serves");
    assert_eq!(
        result.brep.shells.len(),
        1,
        "side walls connect outer and inner skins into one shell"
    );
    // Honest entity accounting (5+5 skin faces + 4 walls; rims are
    // ADOPTED, so edges = 12 outer + 12 inner + 4 verticals, and
    // χ = V−E+F = 16−28+14 = 2 for the one genus-0 shell).
    assert_eq!(result.brep.faces.len(), 14);
    assert_eq!(result.brep.vertices.len(), 16);
    assert_eq!(
        result.brep.edges.len(),
        28,
        "wall rims must be adopted, not duplicated (12 outer + 12 inner + 4 verticals)"
    );
    assert!(
        result.brep.validate().is_ok(),
        "open shell must validate: {:?}",
        result.brep.validate().err()
    );
    // Full screen incl. self-intersection — cheap on 14 faces.
    let findings = result.brep.validate_all();
    assert!(
        findings.is_empty(),
        "open shell must pass validate_all: {findings:?}"
    );
}

/// Thicken of a closed body is a hollow two-shell solid — honest
/// record, validates clean.
#[test]
fn thickened_cube_records_two_shells_and_validates() {
    let cube = BRep::unit_cube(tol()).expect("cube builds");
    let result = thicken(&cube, 0.1, &ShellOptions::default(), tol()).expect("thicken serves");
    assert_eq!(result.brep.shells.len(), 2);
    assert!(
        result.brep.validate().is_ok(),
        "thickened cube must validate: {:?}",
        result.brep.validate().err()
    );
}
