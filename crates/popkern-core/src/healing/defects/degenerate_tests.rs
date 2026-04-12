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
