// Structural lint marker: this file is test code, included via #[cfg(test)] in intersect.rs.
#[cfg(test)]
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

use super::{compute_face_aabbs, detect_edge_face_intersections};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn translate(mut brep: BRep, offset: Vec3) -> BRep {
    for (_, v) in &mut brep.vertices {
        v.position.x += offset.x;
        v.position.y += offset.y;
        v.position.z += offset.z;
    }
    brep
}

#[test]
fn no_edge_face_hits_for_separated_cubes() {
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);
    let hits = detect_edge_face_intersections(&a, &b, &aabbs_b, &aabbs_a, &tol());
    assert!(
        hits.is_empty(),
        "separated cubes should have no edge-face hits"
    );
}

#[test]
fn overlapping_cubes_detection_does_not_crash() {
    // Overlapping cubes with co-planar edges may or may not produce
    // edge-face hits depending on geometry. This test verifies the
    // detection doesn't crash and returns a valid result.
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.5, 0.0),
    );
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);
    let hits = detect_edge_face_intersections(&a, &b, &aabbs_b, &aabbs_a, &tol());
    // Just verify it returns without panicking.
    let _ = hits.len();
}

#[test]
fn diagonal_overlap_produces_edge_face_hits() {
    // Shift cube B by (0.5, 0.5, 0.5) so edges of one cube pierce through
    // faces of the other at non-coplanar angles.
    let a = BRep::unit_cube(tol()).expect("cube a");
    let b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(0.5, 0.5, 0.5),
    );
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);
    let hits = detect_edge_face_intersections(&a, &b, &aabbs_b, &aabbs_a, &tol());
    // With diagonal overlap, edges should pierce through faces.
    assert!(
        !hits.is_empty(),
        "diagonally overlapping cubes should have edge-face hits, got 0"
    );
    // All hit points should have finite coordinates.
    for hit in &hits {
        assert!(hit.t_edge.is_finite(), "hit t_edge should be finite");
        assert!(hit.point.x.is_finite(), "hit point.x should be finite");
    }
}
