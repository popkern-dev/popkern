//! Tests for the intersection hierarchy.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::intersect::compute_face_aabbs;
use crate::tolerance::BoolTolerance;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn translate(brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::Transform3;
    let t = Transform3::translation(offset.x, offset.y, offset.z);
    brep.transform(&t).unwrap()
}

#[test]
fn overlapping_cubes_have_vv_hits() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(1.0, 0.0, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let result = super::compute_hierarchy(&a, &b, &aabbs_a, &aabbs_b, &bt);

    // Cube A has vertices at x=1, cube B has vertices at x=1 (offset by 1.0).
    // So 4 vertices of A at x=1 should coincide with 4 vertices of B at x=0+1=1.
    assert!(
        result.vv_hits.len() >= 4,
        "expected >= 4 VV hits, got {}",
        result.vv_hits.len()
    );
}

#[test]
fn disjoint_cubes_have_no_hits() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let result = super::compute_hierarchy(&a, &b, &aabbs_a, &aabbs_b, &bt);

    assert!(
        result.is_empty(),
        "expected no hits for disjoint cubes: vv={} ve={} ee={} vf={} ef={}",
        result.vv_hits.len(),
        result.ve_hits.len(),
        result.ee_hits.len(),
        result.vf_hits.len(),
        result.ef_hits.len(),
    );
}

#[test]
fn half_offset_cubes_have_ef_hits() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.5, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let result = super::compute_hierarchy(&a, &b, &aabbs_a, &aabbs_b, &bt);

    // Edges of B pierce faces of A and vice versa
    assert!(
        !result.ef_hits.is_empty(),
        "expected EF hits for half-offset cubes"
    );
}

#[test]
fn hierarchy_result_total_hits() {
    let mut r = super::result::HierarchyResult::default();
    assert_eq!(r.total_hits(), 0);
    assert!(r.is_empty());

    r.vv_hits.push(super::result::VertexVertexHit {
        vertex_a: popkern_core::topology::VertexHandle::default(),
        vertex_b: popkern_core::topology::VertexHandle::default(),
        distance: 0.0,
    });
    assert_eq!(r.total_hits(), 1);
    assert!(!r.is_empty());
}
