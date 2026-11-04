//! Sheets in the partition as laminae (PK-978, stage 9 item 1).
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{BodySet, partition};

/// A planar sheet at height `z` over x in [-1, x1], y in [-1, 2]:
/// the core's two-face quad in the one-face sheet form.
fn sheet(x1: f64, z: f64) -> BRep {
    let quad = BRep::sheet_quad(-1.0, x1, -1.0, 2.0, z, Tolerance::default()).expect("quad");
    crate::sheet::lamina::single_face_sheet(&quad).expect("one-face sheet")
}

fn pieces(target: &BRep, sheets: &[&BRep]) -> Vec<usize> {
    let tol = Tolerance::default();
    let cut = crate::split(target, sheets, tol).expect("split");
    let mut faces: Vec<usize> = cut
        .pieces
        .iter()
        .map(|p| {
            assert!(
                p.body.validate_all().is_empty(),
                "{:?}",
                p.body.validate_all()
            );
            p.body.faces.len()
        })
        .collect();
    faces.sort_unstable();
    faces
}

#[test]
fn a_sheet_through_a_box_cuts_it_in_two() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    assert_eq!(pieces(&cube, &[&sheet(2.0, 0.5)]), vec![6, 6]);
    assert_eq!(pieces(&cube, &[&sheet(2.0, 0.3)]), vec![6, 6]);
    // Two sheets: three slabs.
    assert_eq!(
        pieces(&cube, &[&sheet(2.0, 0.25), &sheet(2.0, 0.75)]),
        vec![6, 6, 6]
    );
}

#[test]
fn a_sheet_ending_inside_a_box_leaves_it_whole() {
    // The sheet stops at x = 0.5: its records dangle in the box's y = 0
    // and y = 1 faces, kept as slits, and the cells above and below it
    // meet round its free edge.
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    assert_eq!(pieces(&cube, &[&sheet(0.5, 0.5)]), vec![6]);
}

#[test]
fn a_sheet_bounds_no_volume() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let p = partition(&[&cube, &sheet(2.0, 0.5)], tol).expect("partition");
    let regions = p.regions();
    assert_eq!(regions.len(), 2, "the box above and below the sheet");
    assert!(regions.iter().all(|r| r.bodies == BodySet::of(&[0])));
}

#[test]
fn a_sheet_floating_inside_a_box_leaves_it_whole() {
    // A sheet wholly inside the box touches no box face: its two sides
    // meet round its free edge, and the box stays one piece.
    let quad = BRep::sheet_quad(0.25, 0.75, 0.25, 0.75, 0.5, Tolerance::default()).expect("quad");
    let inner = crate::sheet::lamina::single_face_sheet(&quad).expect("one-face sheet");
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    assert_eq!(pieces(&cube, &[&inner]), vec![6]);
    // One region: the sheet's shell is an inner shell of the box's.
    let p = partition(&[&cube, &inner], Tolerance::default()).expect("partition");
    assert_eq!(p.regions().len(), 1);
}
