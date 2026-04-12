use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_has_no_free_edges() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let free = find_free_edges(&brep);
    assert!(
        free.is_empty(),
        "manifold cube should have no free edges, found {}",
        free.len()
    );
}

#[test]
fn stitch_on_valid_cube_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let result = stitch_boundaries(&mut brep, 0.01, 0.005, &tol);
    assert_eq!(result.free_edges_before, 0);
    assert_eq!(result.free_edges_after, 0);
    assert_eq!(result.pairs_matched, 0);
    assert_eq!(result.pairs_merged, 0);
    assert!(result.unmatched_edges.is_empty());
    assert!(result.ambiguous_edges.is_empty());
    assert!(result.max_displacement < f64::EPSILON);
}

#[test]
fn empty_brep_has_no_free_edges() {
    let brep = BRep::new(Tolerance::default());
    let free = find_free_edges(&brep);
    assert!(free.is_empty());
}

#[test]
fn stitch_on_empty_brep_is_noop() {
    let mut brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let result = stitch_boundaries(&mut brep, 0.01, 0.005, &tol);
    assert_eq!(result.free_edges_before, 0);
    assert_eq!(result.pairs_matched, 0);
    assert!(result.ambiguous_edges.is_empty());
}

#[test]
fn find_matching_edges_empty_list_returns_empty() {
    let brep = BRep::new(Tolerance::default());
    let matches = find_matching_edges(&brep, &[], 0.01, 8);
    assert!(matches.result.is_empty());
    assert!(matches.ambiguous.is_empty());
}

#[test]
fn directed_hausdorff_identical_points_is_zero() {
    let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    let h = directed_hausdorff(&pts, &pts);
    assert!(h < f64::EPSILON);
}

#[test]
fn directed_hausdorff_known_value() {
    let pts_a = vec![Point3::new(0.0, 0.0, 0.0)];
    let pts_b = vec![Point3::new(1.0, 0.0, 0.0)];
    let h = directed_hausdorff(&pts_a, &pts_b);
    assert!((h - 1.0).abs() < f64::EPSILON);
}

#[test]
fn detect_partial_overlap_no_geometry() {
    let brep = BRep::new(Tolerance::default());
    // With no edges, should return None.
    let eh = EdgeHandle::default();
    let result = detect_partial_overlap(&brep, eh, eh, 0.01);
    assert!(result.is_none());
}
