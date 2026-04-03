//! Tests for the coincidence engine.

use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::intersect::compute_face_aabbs;
use crate::tolerance::BoolTolerance;

use super::detect::detect_coincidences;
use super::types::CoincidenceOrientation;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn translate(brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::Transform3;
    let t = Transform3::translation(offset.x, offset.y, offset.z);
    brep.transform(&t).unwrap()
}

#[test]
fn identical_cubes_have_face_overlaps() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let map = detect_coincidences(&a, &b, &aabbs_a, &aabbs_b, &bt);

    // Two identical cubes: all 6 face pairs should be coincident
    assert!(
        map.face_overlaps.len() >= 6,
        "expected >= 6 face overlaps for identical cubes, got {}",
        map.face_overlaps.len()
    );

    // All should be same-oriented (identical cubes have matching normals)
    for fo in &map.face_overlaps {
        assert_eq!(
            fo.orientation,
            CoincidenceOrientation::Same,
            "identical cubes should have same orientation"
        );
    }
}

#[test]
fn disjoint_cubes_have_no_coincidences() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let map = detect_coincidences(&a, &b, &aabbs_a, &aabbs_b, &bt);

    assert!(
        map.face_overlaps.is_empty(),
        "disjoint cubes should have no face overlaps"
    );
}

#[test]
fn touching_cubes_have_one_face_overlap() {
    // Cube A at [0,1]^3, Cube B at [1,2] x [0,1] x [0,1]
    // They share the face at x=1
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(1.0, 0.0, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let map = detect_coincidences(&a, &b, &aabbs_a, &aabbs_b, &bt);

    // Should have at least 1 face overlap (the shared face at x=1)
    assert!(
        !map.face_overlaps.is_empty(),
        "touching cubes should have face overlap at shared face"
    );

    // The shared face should be opposite-oriented (A's x=1 face has normal +x,
    // B's x=0 face has normal -x)
    let has_opposite = map
        .face_overlaps
        .iter()
        .any(|fo| fo.orientation == CoincidenceOrientation::Opposite);
    assert!(
        has_opposite,
        "touching cubes should have opposite-oriented overlap at shared face"
    );
}

#[test]
fn vertex_merges_for_touching_cubes() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(1.0, 0.0, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let map = detect_coincidences(&a, &b, &aabbs_a, &aabbs_b, &bt);

    // 4 vertices at x=1 on A should merge with 4 vertices at x=0 on B
    assert!(
        map.vertex_merges.len() >= 4,
        "expected >= 4 vertex merges, got {}",
        map.vertex_merges.len()
    );
}

#[test]
fn coincidence_map_lookup() {
    let mut map = super::types::CoincidenceMap::default();
    assert!(map.is_empty());
    assert_eq!(map.total(), 0);

    let fa = popkern_core::topology::FaceHandle::default();
    let fb = popkern_core::topology::FaceHandle::default();

    map.face_overlaps.push(super::types::FaceOverlap {
        face_a: fa,
        face_b: fb,
        orientation: CoincidenceOrientation::Same,
        full_overlap: true,
        overlap_center: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
    });

    assert!(!map.is_empty());
    assert_eq!(map.total(), 1);

    let orient = map.face_orientation(fa, fb);
    assert_eq!(orient, Some(CoincidenceOrientation::Same));
}
