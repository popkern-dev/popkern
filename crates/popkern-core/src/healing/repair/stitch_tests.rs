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

// -----------------------------------------------------------------------
// Helpers: build BReps with free boundary edges
// -----------------------------------------------------------------------

use crate::constructors::builder::BRepBuilder;
use crate::geometry::Vec3;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::topology::EdgeHandle;

/// Build a single-face BRep (one quad). All 4 edges are free boundaries
/// because both half-edges of each edge belong to the same face.
fn make_single_face_brep() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);
    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face).unwrap();
    b.build_unchecked()
}

/// Build two adjacent quad patches that share a geometric boundary
/// but are NOT topologically connected (separate vertices at the seam).
/// Edges have Line curve geometry so sampling works.
fn make_two_adjacent_patches() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let a0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let a1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let a2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let a3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let b0 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let b1 = b.add_vertex(Point3::new(2.0, 0.0, 0.0));
    let b2 = b.add_vertex(Point3::new(2.0, 1.0, 0.0));
    let b3 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let (ea01, ha01, _) = b.add_edge(a0, a1);
    let (ea12, ha12, _) = b.add_edge(a1, a2);
    let (ea23, ha23, _) = b.add_edge(a2, a3);
    let (ea30, ha30, _) = b.add_edge(a3, a0);
    let (eb01, hb01, _) = b.add_edge(b0, b1);
    let (eb12, hb12, _) = b.add_edge(b1, b2);
    let (eb23, hb23, _) = b.add_edge(b2, b3);
    let (eb30, hb30, _) = b.add_edge(b3, b0);

    let edges_and_verts: Vec<(crate::topology::EdgeHandle, Point3, Point3)> = vec![
        (ea01, Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (ea12, Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (ea23, Point3::new(1.0, 1.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        (ea30, Point3::new(0.0, 1.0, 0.0), Point3::new(0.0, 0.0, 0.0)),
        (eb01, Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)),
        (eb12, Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 1.0, 0.0)),
        (eb23, Point3::new(2.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (eb30, Point3::new(1.0, 1.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
    ];
    for (eh, start, end) in &edges_and_verts {
        let dir = Vec3::new(end.x - start.x, end.y - start.y, end.z - start.z);
        let len = dir.length();
        b.set_edge_curve(
            *eh,
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: *start,
                    direction: dir,
                },
                t_start: 0.0,
                t_end: if len > 0.0 { 1.0 } else { 0.0 },
            },
        );
    }

    let fa = b.add_face(shell);
    b.add_loop(&[ha01, ha12, ha23, ha30], fa).unwrap();
    let fb = b.add_face(shell);
    b.add_loop(&[hb01, hb12, hb23, hb30], fb).unwrap();
    b.build_unchecked()
}

// -----------------------------------------------------------------------
// Tests with defective geometry
// -----------------------------------------------------------------------

#[test]
fn find_free_edges_on_single_face() {
    let brep = make_single_face_brep();
    let free = find_free_edges(&brep);
    assert_eq!(
        free.len(),
        4,
        "single face should have 4 free boundary edges"
    );
}

#[test]
fn stitch_single_face_no_matches() {
    let mut brep = make_single_face_brep();
    let tol = Tolerance::default();
    let result = stitch_boundaries(&mut brep, 0.1, 0.05, &tol);
    assert_eq!(result.free_edges_before, 4);
    assert_eq!(result.pairs_matched, 0);
    assert_eq!(result.pairs_merged, 0);
}

#[test]
fn two_patches_have_8_free_edges() {
    let brep = make_two_adjacent_patches();
    let free = find_free_edges(&brep);
    assert_eq!(
        free.len(),
        8,
        "two disconnected patches should have 8 free edges"
    );
}

#[test]
fn find_matching_edges_detects_coincident_boundary() {
    let brep = make_two_adjacent_patches();
    let free = find_free_edges(&brep);
    let matches = find_matching_edges(&brep, &free, 0.01, 8);
    assert!(
        !matches.result.is_empty(),
        "coincident boundary edges should be matched"
    );
    assert_eq!(matches.result.len(), 1);
    assert!(
        matches.result[0].distance < 1e-10,
        "coincident edges should have near-zero distance, got {}",
        matches.result[0].distance
    );
}

#[test]
fn stitch_merges_coincident_boundary() {
    let mut brep = make_two_adjacent_patches();
    let tol = Tolerance::default();
    let result = stitch_boundaries(&mut brep, 0.01, 0.005, &tol);
    assert_eq!(result.free_edges_before, 8);
    assert_eq!(result.pairs_matched, 1);
    assert_eq!(result.pairs_merged, 1);
    assert!(
        result.max_displacement < 1e-10,
        "coincident boundary merge should have negligible displacement"
    );
}

#[test]
fn stitch_respects_max_displacement_limit() {
    let mut b = BRepBuilder::new(Tolerance::default());
    let a0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let a1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let a2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let a3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let b0 = b.add_vertex(Point3::new(1.5, 0.0, 0.0));
    let b1 = b.add_vertex(Point3::new(2.5, 0.0, 0.0));
    let b2 = b.add_vertex(Point3::new(2.5, 1.0, 0.0));
    let b3 = b.add_vertex(Point3::new(1.5, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let (ea01, ha01, _) = b.add_edge(a0, a1);
    let (ea12, ha12, _) = b.add_edge(a1, a2);
    let (ea23, ha23, _) = b.add_edge(a2, a3);
    let (ea30, ha30, _) = b.add_edge(a3, a0);
    let (eb01, hb01, _) = b.add_edge(b0, b1);
    let (eb12, hb12, _) = b.add_edge(b1, b2);
    let (eb23, hb23, _) = b.add_edge(b2, b3);
    let (eb30, hb30, _) = b.add_edge(b3, b0);

    let edges = [
        (ea01, Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (ea12, Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (ea23, Point3::new(1.0, 1.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        (ea30, Point3::new(0.0, 1.0, 0.0), Point3::new(0.0, 0.0, 0.0)),
        (eb01, Point3::new(1.5, 0.0, 0.0), Point3::new(2.5, 0.0, 0.0)),
        (eb12, Point3::new(2.5, 0.0, 0.0), Point3::new(2.5, 1.0, 0.0)),
        (eb23, Point3::new(2.5, 1.0, 0.0), Point3::new(1.5, 1.0, 0.0)),
        (eb30, Point3::new(1.5, 1.0, 0.0), Point3::new(1.5, 0.0, 0.0)),
    ];
    for (eh, start, end) in &edges {
        let dir = Vec3::new(end.x - start.x, end.y - start.y, end.z - start.z);
        b.set_edge_curve(
            *eh,
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: *start,
                    direction: dir,
                },
                t_start: 0.0,
                t_end: 1.0,
            },
        );
    }

    let fa = b.add_face(shell);
    b.add_loop(&[ha01, ha12, ha23, ha30], fa).unwrap();
    let fb = b.add_face(shell);
    b.add_loop(&[hb01, hb12, hb23, hb30], fb).unwrap();
    let mut brep = b.build_unchecked();

    let tol = Tolerance::default();
    let result = stitch_boundaries(&mut brep, 1.0, 0.01, &tol);
    assert_eq!(result.free_edges_before, 8);
    assert_eq!(
        result.pairs_merged, 0,
        "should not merge when gap exceeds max_displacement"
    );
}

#[test]
fn directed_hausdorff_asymmetric() {
    let pts_a = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    let pts_b = vec![Point3::new(0.5, 0.0, 0.0)];
    let ab = directed_hausdorff(&pts_a, &pts_b);
    let ba = directed_hausdorff(&pts_b, &pts_a);
    assert!((ab - 0.5).abs() < 1e-10);
    assert!((ba - 0.5).abs() < 1e-10);
}

#[test]
fn sample_edge_points_linear_fallback() {
    let brep = make_single_face_brep();
    let edge = brep.edges.keys().next().unwrap();
    let pts = sample_edge_points(&brep, edge, 3);
    assert!(pts.is_some());
    let pts = pts.unwrap();
    assert_eq!(pts.len(), 3);
}
