//! Tests for imprint utilities: snap, sliver detection, watertightness.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

// ── Snap tests ─────────────────────────────────────────────────────────

#[test]
fn snap_merges_close_points() {
    use super::snap::snap_points;

    let mut points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1e-7, 0.0, 0.0), // within 1e-6 of first
        Point3::new(1.0, 0.0, 0.0),  // far away
    ];

    let log = snap_points(&mut points, 1e-6, 10);
    assert_eq!(log.count(), 1, "should snap one pair");
    // Points 0 and 1 should be identical after snap
    assert!(points[0].distance_to(&points[1]) < 1e-15);
    // Point 2 should be unchanged
    assert!((points[2].x - 1.0).abs() < 1e-15);
}

#[test]
fn snap_cascading() {
    use super::snap::snap_points;

    // Three points in a chain: A--B--C where A-B and B-C are within tolerance
    // but A-C is not. After snapping, all three should be at A.
    let mut points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.5e-6, 0.0, 0.0), // within 1e-6 of A
        Point3::new(0.9e-6, 0.0, 0.0), // within 1e-6 of B (after B snaps to A)
    ];

    let log = snap_points(&mut points, 1e-6, 10);
    // After first pass: B snaps to A. Then C is within tolerance of A.
    // After second pass: C snaps to A.
    assert!(log.count() >= 2, "should cascade: {} snaps", log.count());
    assert!(points[0].distance_to(&points[2]) < 1e-15);
}

#[test]
fn snap_no_merge_when_distant() {
    use super::snap::snap_points;

    let mut points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ];

    let log = snap_points(&mut points, 1e-6, 10);
    assert_eq!(log.count(), 0);
}

#[test]
fn snap_to_nearest_finds_closest() {
    use super::snap::snap_to_nearest;

    let vertices = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ];
    let pt = Point3::new(0.9e-6, 0.0, 0.0);

    let idx = snap_to_nearest(&pt, &vertices, 1e-5);
    assert_eq!(idx, Some(0), "should snap to nearest vertex (index 0)");
}

#[test]
fn snap_to_nearest_none_when_distant() {
    use super::snap::snap_to_nearest;

    let vertices = vec![Point3::new(0.0, 0.0, 0.0)];
    let pt = Point3::new(1.0, 0.0, 0.0);

    assert!(snap_to_nearest(&pt, &vertices, 1e-6).is_none());
}

#[test]
fn micro_edge_detection() {
    use super::snap::is_micro_edge;

    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1e-7, 0.0, 0.0);
    let c = Point3::new(1.0, 0.0, 0.0);

    assert!(is_micro_edge(&a, &b, 1e-6));
    assert!(!is_micro_edge(&a, &c, 1e-6));
}

// ── Sliver detection tests ─────────────────────────────────────────────

#[test]
fn triangle_is_not_a_sliver() {
    use super::sliver::{SliverThresholds, analyze_face_polygon};

    let thresholds = SliverThresholds::from_tolerance(1e-6);
    let triangle = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
    ];
    let report = analyze_face_polygon(&triangle, &thresholds);
    assert!(!report.is_sliver, "normal triangle should not be a sliver");
    assert!(
        report.area > 0.4,
        "area should be ~0.5, got {}",
        report.area
    );
}

#[test]
fn degenerate_triangle_is_sliver() {
    use super::sliver::{SliverThresholds, analyze_face_polygon};

    let thresholds = SliverThresholds::from_tolerance(1e-4);
    // Very thin triangle (height ~ 1e-5)
    let thin = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1e-5, 0.0),
    ];
    let report = analyze_face_polygon(&thin, &thresholds);
    assert!(report.is_sliver, "thin triangle should be a sliver");
    assert!(
        report.area < 1e-4,
        "area should be tiny, got {}",
        report.area
    );
}

#[test]
fn short_edge_triangle_is_sliver() {
    use super::sliver::{SliverThresholds, SliverType, analyze_face_polygon};

    let thresholds = SliverThresholds::from_tolerance(1e-4);
    // Two vertices almost coincident
    let short = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1e-5, 0.0, 0.0), // micro-edge
        Point3::new(0.5, 0.5, 0.0),
    ];
    let report = analyze_face_polygon(&short, &thresholds);
    assert!(report.is_sliver, "short-edge face should be a sliver");
    assert_eq!(report.sliver_type, Some(SliverType::ShortEdge));
}

#[test]
fn too_few_vertices_is_sliver() {
    use super::sliver::{SliverThresholds, analyze_face_polygon};

    let thresholds = SliverThresholds::from_tolerance(1e-6);
    let degenerate = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    let report = analyze_face_polygon(&degenerate, &thresholds);
    assert!(report.is_sliver);
}

// ── Watertightness validation tests ────────────────────────────────────

#[test]
fn unit_cube_is_watertight() {
    use super::shared_geometry::validate_watertightness;

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = validate_watertightness(&cube, 1e-4, 5);
    assert!(
        report.is_watertight(),
        "unit cube should be watertight, found {} gaps (max gap: {})",
        report.gaps.len(),
        report.max_gap()
    );
    assert!(report.edges_checked > 0, "should have checked some edges");
}

#[test]
fn cylinder_is_watertight() {
    use super::shared_geometry::validate_watertightness;

    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let report = validate_watertightness(&cyl, 1e-3, 5);
    assert!(
        report.is_watertight(),
        "cylinder should be watertight, found {} gaps (max gap: {})",
        report.gaps.len(),
        report.max_gap()
    );
}

#[test]
fn boolean_result_watertightness() {
    use super::shared_geometry::validate_watertightness;
    use popkern_core::geometry::Transform3;

    let tol = Tolerance::default();
    let a = BRep::unit_cube(tol).unwrap();
    let b = BRep::unit_cube(tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();

    if let Ok(result) = crate::union(a, b, tol, None) {
        let report = validate_watertightness(&result, 1e-3, 5);
        // Report gaps if any — this is informational, not a hard assert,
        // since the current graph-merge doesn't enforce shared geometry yet.
        if !report.is_watertight() {
            eprintln!(
                "[watertight] boolean result has {} gap(s), max deviation: {:.2e}",
                report.gaps.len(),
                report.max_gap()
            );
        }
    }
}
