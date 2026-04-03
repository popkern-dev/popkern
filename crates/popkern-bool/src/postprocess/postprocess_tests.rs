//! Tests for post-processing operations.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ── Redundant edge tests ───────────────────────────────────────────────

#[test]
fn cube_has_no_redundant_edges() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let report = super::redundant_edge::find_redundant_edges(&cube, &tol());
    assert!(
        report.redundant_edges.is_empty(),
        "cube should have no redundant edges, found {}",
        report.redundant_edges.len()
    );
    assert_eq!(report.edges_analyzed, 12); // cube has 12 edges
}

#[test]
fn cylinder_has_no_redundant_edges() {
    let cyl = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let report = super::redundant_edge::find_redundant_edges(&cyl, &tol());
    assert!(report.redundant_edges.is_empty());
    assert!(report.edges_analyzed > 0);
}

// ── Face simplify tests ────────────────────────────────────────────────

#[test]
fn cube_has_no_merge_candidates() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let report = super::face_simplify::find_merge_candidates(&cube, &tol());
    assert!(
        report.merge_groups.is_empty(),
        "cube faces are all on different planes"
    );
    assert_eq!(report.potential_reduction(), 0);
    assert_eq!(report.faces_analyzed, 6);
}

#[test]
fn face_simplify_report_metadata() {
    let report = super::face_simplify::FaceSimplifyReport::default();
    assert_eq!(report.potential_reduction(), 0);
    assert!(report.merge_groups.is_empty());
}

// ── Tests on boolean results ───────────────────────────────────────────

#[test]
fn redundant_edges_on_boolean_result() {
    use popkern_core::geometry::Transform3;

    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    if let Ok(result) = crate::union(a, b, tol(), None) {
        let report = super::redundant_edge::find_redundant_edges(&result, &tol());
        // Log for diagnostic — boolean result may or may not have redundant edges
        // depending on whether co-surfacial faces survive assembly
        eprintln!(
            "[postprocess] boolean union: {} edges analyzed, {} redundant",
            report.edges_analyzed,
            report.redundant_edges.len()
        );
        assert!(report.edges_analyzed > 0);
    }
}

#[test]
fn face_simplify_on_boolean_result() {
    use popkern_core::geometry::Transform3;

    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    if let Ok(result) = crate::union(a, b, tol(), None) {
        let report = super::face_simplify::find_merge_candidates(&result, &tol());
        eprintln!(
            "[postprocess] boolean union: {} faces, {} merge groups, {} potential reduction",
            report.faces_analyzed,
            report.merge_groups.len(),
            report.potential_reduction()
        );
        assert!(report.faces_analyzed > 0);
    }
}
