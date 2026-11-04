//! Sheet × solid booleans through the cell complex (PK-978): a planar
//! and a curved sheet against a box, measured against closed forms.

use popkern_bool::{BoolOperation, SheetBooleanResult, boolean_with_sheets};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};
use popkern_measure::{body_volume, face_area};

const R: f64 = 0.3;

/// An open tube of radius `R` along z from -0.5 to 1.5 through the
/// unit cube's centre: a cylinder with its caps retired.
fn tube() -> BRep {
    let tol = Tolerance::default();
    let mut c = BRep::cylinder(R, 2.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.5, 0.5, -0.5)))
        .expect("cylinder");
    let caps: Vec<_> = c
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. })
            )
        })
        .map(|(h, _)| h)
        .collect();
    assert_eq!(caps.len(), 2);
    c.retire_faces_to_boundary(&caps).expect("open tube");
    assert_eq!(c.body_kind(), BodyKind::Sheet);
    c
}

fn sheet_area(b: &BRep) -> f64 {
    let tol = Tolerance::default();
    b.faces
        .keys()
        .map(|f| face_area(b, f, &tol).expect("area"))
        .sum()
}

fn open_shell(r: Result<SheetBooleanResult, popkern_bool::BoolError>) -> BRep {
    match r.expect("sheet boolean") {
        SheetBooleanResult::OpenShell(b) => {
            assert_eq!(b.body_kind(), BodyKind::Sheet);
            let findings = b.validate_all();
            assert!(findings.is_empty(), "{findings:?}");
            b
        }
        other => panic!("expected an open shell, got {other:?}"),
    }
}

fn close(got: f64, want: f64, rel: f64) {
    assert!(
        (got - want).abs() <= rel * want.abs(),
        "{got} != {want} (rel {rel})"
    );
}

#[test]
fn a_tube_cuts_a_box_into_a_core_and_a_bored_block() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let pieces =
        match boolean_with_sheets(BoolOperation::Subtraction, &cube, &tube(), tol).expect("cut") {
            SheetBooleanResult::Pieces(p) => p,
            other => panic!("expected pieces, got {other:?}"),
        };
    let mut volumes: Vec<f64> = pieces
        .iter()
        .map(|p| {
            let findings = p.validate_all();
            assert!(findings.is_empty(), "{findings:?}");
            body_volume(p, &tol).expect("volume")
        })
        .collect();
    volumes.sort_by(f64::total_cmp);
    let core = std::f64::consts::PI * R * R;
    assert_eq!(volumes.len(), 2);
    close(volumes[0], core, 1e-3);
    close(volumes[1], 1.0 - core, 1e-3);
}

#[test]
fn a_tube_parts_inside_and_outside_a_box() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let wall = 2.0 * std::f64::consts::PI * R;
    let inside = open_shell(boolean_with_sheets(
        BoolOperation::Intersection,
        &tube(),
        &cube,
        tol,
    ));
    close(sheet_area(&inside), wall, 1e-3);
    let outside = open_shell(boolean_with_sheets(
        BoolOperation::Subtraction,
        &tube(),
        &cube,
        tol,
    ));
    close(sheet_area(&outside), wall, 1e-3);
}

#[test]
fn a_constructor_quad_parts_against_a_box_either_order() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let quad = BRep::sheet_quad(-1.0, 2.0, -1.0, 2.0, 0.5, tol).expect("quad");
    for (a, b) in [(&quad, &cube), (&cube, &quad)] {
        let inside = open_shell(boolean_with_sheets(BoolOperation::Intersection, a, b, tol));
        close(sheet_area(&inside), 1.0, 1e-9);
    }
    let outside = open_shell(boolean_with_sheets(
        BoolOperation::Subtraction,
        &quad,
        &cube,
        tol,
    ));
    close(sheet_area(&outside), 8.0, 1e-9);
    let pieces =
        match boolean_with_sheets(BoolOperation::Subtraction, &cube, &quad, tol).expect("cut") {
            SheetBooleanResult::Pieces(p) => p,
            other => panic!("expected pieces, got {other:?}"),
        };
    let mut volumes: Vec<f64> = pieces
        .iter()
        .map(|p| body_volume(p, &tol).expect("volume"))
        .collect();
    volumes.sort_by(f64::total_cmp);
    assert_eq!(volumes.len(), 2);
    close(volumes[0], 0.5, 1e-9);
    close(volumes[1], 0.5, 1e-9);
}
