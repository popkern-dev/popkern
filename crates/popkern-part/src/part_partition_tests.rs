//! Partitions, splits and imprints at the Part level.

use popkern_bool::{FaceClass, SetExpr};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::part::Part;

fn cube_part(x: f64, y: f64, z: f64) -> Part {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube");
    Part::from_brep(brep, tol, "cube").expect("part")
}

/// Two unit cubes overlapping in a 0.5 × 0.7 × 0.8 box: the overlap's
/// volume is 0.28, each cube alone keeps 0.72.
#[test]
fn two_overlapping_cubes_partition_by_volume() {
    let (a, b) = (cube_part(0.0, 0.0, 0.0), cube_part(0.5, 0.3, 0.2));
    let p = Part::partition(&[&a, &b]).expect("partition");
    let mut volumes: Vec<f64> = p
        .regions()
        .expect("regions")
        .iter()
        .map(|r| r.volume)
        .collect();
    volumes.sort_by(f64::total_cmp);
    let want = [0.28, 0.72, 0.72];
    assert_eq!(volumes.len(), 3);
    for (got, want) in volumes.iter().zip(want) {
        assert!((got - want).abs() < 1e-9, "{volumes:?}");
    }
    let overlap = p
        .select(&SetExpr::All, &popkern_bool::BoolOptions::default())
        .expect("overlap");
    let v = popkern_measure::body_volume(overlap.brep(), &Tolerance::default()).expect("volume");
    assert!((v - 0.28).abs() < 1e-9, "{v}");
    assert_eq!(p.region_parts(&SetExpr::Any).expect("parts").len(), 3);
}

#[test]
fn a_part_splits_into_pieces_inside_and_outside_a_tool() {
    let target = cube_part(0.0, 0.0, 0.0);
    let tool = cube_part(0.5, 0.3, 0.2);
    let pieces = target.split(&[&tool]).expect("split");
    assert_eq!(pieces.len(), 2);
    assert_eq!(
        pieces.iter().filter(|p| p.inside_tools.contains(0)).count(),
        1
    );
}

#[test]
fn a_part_imprints_and_classifies_its_faces() {
    let target = cube_part(0.0, 0.0, 0.0);
    let tool = cube_part(0.5, 0.5, 0.5);
    let imprinted = target.imprint(&tool).expect("imprint");
    assert_eq!(imprinted.faces_classified(FaceClass::Inside).count(), 3);
    assert_eq!(imprinted.part.face_count(), 9);
}
