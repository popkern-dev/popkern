use super::*;
use crate::error::ValidationError;
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, FaceSide, RegionKind};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Killing a cube's solid region leaves its surface with both sides in
/// the outside's shell, which separates space and is no valid cell
/// structure (a shell of two pieces); making the face region behind the
/// faces gives the solid back.
#[test]
fn a_cube_region_is_killed_and_made_again() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let solid_counts = NonManifoldCounts::of(&cube);
    let solid = cube.solid_regions()[0];
    cube.kill_face_region(solid).expect("kill");
    assert_eq!(
        solid_counts.delta(NonManifoldCounts::of(&cube)),
        [0, 0, 0, 0, 0, -1, -1, 0, 0]
    );
    let shell = cube
        .face_back_shell(cube.faces.keys().next().expect("face"))
        .expect("shell");
    assert!(
        cube.validate_shell_cells_all()
            .iter()
            .any(|e| matches!(e, ValidationError::ShellDisconnected { pieces: 2, .. })),
        "a closed surface in one shell is two pieces"
    );
    cube.make_face_region(shell, FaceSide::Back, RegionKind::Solid)
        .expect("fill");
    assert_eq!(NonManifoldCounts::of(&cube), solid_counts);
    assert_eq!(cube.body_kind(), BodyKind::Solid);
    let errors = cube.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}

/// The void inside a closed surface is a region too: a closed sheet,
/// void on both sides of every face.
#[test]
fn a_closed_sheet_bounds_a_void_region() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let solid = cube.solid_regions()[0];
    cube.kill_face_region(solid).expect("kill");
    let shell = cube
        .face_back_shell(cube.faces.keys().next().expect("face"))
        .expect("shell");
    let inside = cube
        .make_face_region(shell, FaceSide::Back, RegionKind::Void)
        .expect("void inside");
    assert_eq!(cube.region_kind(inside), Some(RegionKind::Void));
    assert_eq!(cube.body_kind(), BodyKind::Sheet);
    let errors = cube.validate_shell_cells_all();
    assert!(errors.is_empty(), "{errors:?}");
}

/// A shell that is not a closed surface on both sides, the outside and
/// a solid's region whose shell is not alone decline.
#[test]
fn face_regions_decline_on_open_or_one_sided_shells() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let solid = cube.solid_regions()[0];
    let shell = cube.region_outer_shell(solid).expect("shell");
    assert_eq!(
        cube.make_face_region(shell, FaceSide::Back, RegionKind::Solid),
        Err(EulerError::NotAClosedSurface { shell })
    );
    let outside = cube.outside;
    assert_eq!(
        cube.kill_face_region(outside),
        Err(EulerError::NotAFaceRegion { region: outside })
    );
}
