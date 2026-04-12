//! Tests for parallel intersection infrastructure.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::intersect::compute_face_aabbs;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn collect_candidates_overlapping() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let candidates = super::ff_parallel::collect_candidates(&aabbs_a, &aabbs_b, tol().linear());
    assert!(!candidates.is_empty(), "should find overlapping pairs");
}

#[test]
fn collect_candidates_disjoint() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let candidates = super::ff_parallel::collect_candidates(&aabbs_a, &aabbs_b, tol().linear());
    assert!(candidates.is_empty());
}

#[test]
fn process_pairs_produces_results() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let candidates = super::ff_parallel::collect_candidates(&aabbs_a, &aabbs_b, tol().linear());
    let results = super::ff_parallel::process_pairs_sequential(&candidates, &a, &b, &tol());

    let total_curves: usize = results.iter().map(|r| r.curves.len()).sum();
    assert!(
        total_curves > 0,
        "overlapping cubes should produce intersection curves"
    );
}

#[test]
fn process_pairs_matches_find_intersections() {
    // Verify parallel infrastructure produces same curve count as the main pipeline
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let candidates = super::ff_parallel::collect_candidates(&aabbs_a, &aabbs_b, tol().linear());
    let parallel_results =
        super::ff_parallel::process_pairs_sequential(&candidates, &a, &b, &tol());
    let parallel_curves: usize = parallel_results.iter().map(|r| r.curves.len()).sum();

    let main_result = crate::intersect::find_intersections(&a, &b, &tol()).unwrap();
    let main_curves = main_result.records.len();

    assert_eq!(
        parallel_curves, main_curves,
        "parallel ({parallel_curves}) should match main pipeline ({main_curves})"
    );
}
