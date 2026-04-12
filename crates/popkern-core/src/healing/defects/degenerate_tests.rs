use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_no_degenerate_defects() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    // Only run small edges, small faces (not missing geometry, since
    // unit_cube has full geometry and we want to test non-degenerate).
    check_small_edges(&brep, tol.linear(), &mut report);
    check_small_faces(&brep, &tol, &mut report);
    assert!(
        report.is_empty(),
        "Valid cube should have no small edges/faces, got {}",
        report.len()
    );
}

#[test]
fn cube_with_geometry_has_no_missing_geometry() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_missing_edge_geometry(&brep, &mut report);
    check_missing_face_geometry(&brep, &mut report);
    assert!(
        report.is_empty(),
        "Unit cube with geometry should have no missing geometry defects"
    );
}

#[test]
fn topology_only_cube_reports_missing_geometry() {
    use crate::geometry::{Point3, Vec3};

    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();

    let mut report = HealingReport::new();
    check_missing_edge_geometry(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Topology-only cube should report missing edge geometry"
    );
    assert_eq!(report.defects()[0].severity, Severity::Info);
}

#[test]
fn small_face_kind_display() {
    assert_eq!(SmallFaceKind::Pin.to_string(), "pin");
    assert_eq!(SmallFaceKind::Sliver.to_string(), "sliver");
    assert_eq!(SmallFaceKind::Bridge.to_string(), "bridge");
    assert_eq!(SmallFaceKind::Slot.to_string(), "slot");
    assert_eq!(SmallFaceKind::Micro.to_string(), "micro");
}

#[test]
fn zero_length_singularity_edges_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_zero_length_singularity_edges(&brep, &tol, &mut report);
    // A valid cube has no zero-length edges.
    assert!(
        report.is_empty(),
        "Valid cube should have no zero-length singularity edges"
    );
}

#[test]
fn empty_brep_no_degenerate_defects() {
    let brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_all(&brep, &tol, &mut report);
    assert!(report.is_empty());
}

// ── Corrupted BRep tests exercising defect-detection error paths ──

#[test]
fn small_edge_detected_via_vertex_fallback() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Remove curve and move vertices close together so the fallback
    // vertex-to-vertex distance path detects a small edge.
    let edge_h = brep.edges.keys().next().unwrap();
    brep.edges[edge_h].curve = None; // force vertex-to-vertex fallback
    let he_h = brep.edges[edge_h].half_edge;
    let he = &brep.half_edges[he_h];
    let twin_h = he.twin;
    let twin = &brep.half_edges[twin_h];
    let v2 = twin.origin;
    let v1_pos = brep.vertices[he.origin].position;
    brep.vertices[v2].position = crate::geometry::Point3::new(v1_pos.x + 1e-8, v1_pos.y, v1_pos.z);
    let mut report = HealingReport::new();
    // Use a threshold larger than the edge length.
    check_small_edges(&brep, 1e-4, &mut report);
    assert!(!report.is_empty(), "Should detect small edge");
    assert_eq!(report.defects()[0].category, DefectCategory::Degenerate);
    assert_eq!(report.defects()[0].severity, Severity::Warning);
    assert!(report.defects()[0].measured_value.is_some());
    assert!(report.defects()[0].threshold.is_some());
}

#[test]
fn small_edge_detected_without_curve() {
    use crate::geometry::{Point3, Vec3};

    // Build a topology-only BRep (no curves) with a very short edge.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1e-8, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();

    let mut report = HealingReport::new();
    check_small_edges(&brep, 1e-4, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect small edge via vertex-to-vertex fallback"
    );
}

#[test]
fn missing_edge_geometry_detected() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge_h = brep.edges.keys().next().unwrap();
    brep.edges[edge_h].curve = None;
    let mut report = HealingReport::new();
    check_missing_edge_geometry(&brep, &mut report);
    assert!(!report.is_empty(), "Should detect missing edge curve");
    assert_eq!(report.defects()[0].severity, Severity::Info);
}

#[test]
fn missing_face_geometry_detected() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    brep.faces[fh].surface = None;
    let mut report = HealingReport::new();
    check_missing_face_geometry(&brep, &mut report);
    assert!(!report.is_empty(), "Should detect missing face surface");
    assert_eq!(report.defects()[0].severity, Severity::Info);
}

#[test]
fn small_face_detected_when_vertices_collapsed() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    // Collapse all outer loop vertices to nearly the same point.
    let target = crate::geometry::Point3::new(0.5, 0.5, 0.5);
    let face_data = &brep.faces[fh];
    let loop_data = &brep.loops[face_data.outer_loop];
    let start_he = loop_data.half_edge;
    let mut cursor = start_he;
    let bound = brep.half_edges.len();
    for _ in 0..=bound {
        let origin = brep.half_edges[cursor].origin;
        let offset = brep.vertices.keys().position(|v| v == origin).unwrap_or(0) as f64;
        brep.vertices[origin].position = crate::geometry::Point3::new(
            target.x + offset * 1e-10,
            target.y + offset * 1e-10,
            target.z,
        );
        let next = brep.half_edges[cursor].next;
        if next == start_he {
            break;
        }
        cursor = next;
    }
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_small_faces(&brep, &tol, &mut report);
    assert!(!report.is_empty(), "Should detect small face");
    assert_eq!(report.defects()[0].category, DefectCategory::Degenerate);
    // The classification description should include one of the SmallFaceKind names.
    let desc = &report.defects()[0].description;
    assert!(
        desc.contains("pin")
            || desc.contains("sliver")
            || desc.contains("micro")
            || desc.contains("bridge")
            || desc.contains("slot"),
        "Description should contain a SmallFaceKind, got: {desc}"
    );
}

#[test]
fn zero_length_singularity_edge_detected() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Make an edge zero-length by moving both endpoint vertices to the same point.
    let edge_h = brep.edges.keys().next().unwrap();
    let he_h = brep.edges[edge_h].half_edge;
    let he = &brep.half_edges[he_h];
    let twin_h = he.twin;
    let twin = &brep.half_edges[twin_h];
    let v1 = he.origin;
    let v2 = twin.origin;
    let midpoint = crate::geometry::Point3::new(0.5, 0.5, 0.5);
    brep.vertices[v1].position = midpoint;
    brep.vertices[v2].position = midpoint;
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_zero_length_singularity_edges(&brep, &tol, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect zero-length singularity edge"
    );
    assert_eq!(report.defects()[0].severity, Severity::Info);
}
