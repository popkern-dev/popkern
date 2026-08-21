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
/// stays ONE shell. (Validation cleanliness of this path is PK-307:
/// the current side-wall builder leaves rim half-edges loopless and
/// duplicates vertical edges, a distinct pre-existing defect.)
#[test]
fn open_shell_stays_one_connected_shell() {
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
