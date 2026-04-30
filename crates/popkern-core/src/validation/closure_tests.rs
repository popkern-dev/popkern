//! Tests for [`super::closure`].  Extracted to keep the implementation
//! file under the 500-line structural cap.

use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn make_cube() -> BRep {
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
    brep
}

#[test]
fn cube_passes_loop_closure() {
    let brep = make_cube();
    assert!(brep.validate_loop_closure().is_ok());
}

#[test]
fn cube_passes_watertight() {
    let brep = make_cube();
    assert!(brep.validate_watertight().is_ok());
}

#[test]
fn triangle_loop_is_closed() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(0.5, 1.0, 0.0))
        .unwrap();
    let _mef = brep.mef(mev2.half_edge_twin, mev1.half_edge).unwrap();
    assert!(brep.validate_loop_closure().is_ok());
}

#[test]
fn loop_closure_with_curve_geometry() {
    // Use a cube with geometry to test the curve evaluation path
    let cube = crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
    assert!(cube.validate_loop_closure().is_ok());
}

#[test]
fn watertight_with_curve_geometry() {
    let cube = crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
    assert!(cube.validate_watertight().is_ok());
}

#[test]
fn loop_closure_on_corrupted_brep_breaks_gracefully() {
    // Remove a half-edge to trigger the `None` break paths in validate_loop_closure
    let mut brep = make_cube();
    let he_handle = brep.half_edges.keys().next().unwrap();
    // Corrupt next pointer to a non-existent handle
    brep.half_edges[he_handle].next = Default::default();
    // Should not panic — just break out of walk
    let _ = brep.validate_loop_closure();
}

#[test]
fn loop_closure_twin_path_with_curve() {
    // Build a cube with geometry, then verify closure works through twin-side curve eval
    let brep = crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
    // The cube has curves on all edges. Loop closure walks all loops and
    // picks the correct curve endpoint for each half-edge by direction
    // matching. Just ensure it passes.
    assert!(brep.validate_loop_closure().is_ok());
}

#[test]
fn corrupted_edge_handle_breaks_gracefully() {
    // Remove an edge to trigger the edge.curve None path
    let mut brep = make_cube();
    let edge_handle = brep.edges.keys().next().unwrap();
    brep.edges.remove(edge_handle);
    // Should not panic
    let _ = brep.validate_loop_closure();
}

#[test]
fn corrupted_vertex_handle_breaks_gracefully() {
    // Remove a vertex to trigger the vertices.get None path
    let mut brep = make_cube();
    let v_handle = brep.vertices.keys().next().unwrap();
    brep.vertices.remove(v_handle);
    // Should not panic
    let _ = brep.validate_loop_closure();
}

#[test]
fn watertight_with_corrupted_he_skips() {
    // Corrupt an edge's half_edge to point to nonexistent handle
    let mut brep =
        crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
    let edge_handle = brep.edges.keys().next().unwrap();
    brep.edges[edge_handle].half_edge = Default::default();
    // Should not panic; just skip the corrupted edge
    let _ = brep.validate_watertight();
}

#[test]
fn watertight_with_bad_gap() {
    // Create a cube with geometry, then corrupt a vertex position to create a gap
    let mut brep =
        crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

    // Move a vertex far away to create a gap
    let v_handle = brep.vertices.keys().next().unwrap();
    brep.vertices[v_handle].position = Point3::new(100.0, 100.0, 100.0);

    let result = brep.validate_watertight();
    assert!(result.is_err());
}

#[test]
fn watertight_twin_vertex_missing_skips() {
    // Remove the twin's origin vertex to trigger the vertices.get(twin.origin) None path
    let mut brep =
        crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

    // Find an edge with a curve
    let edge_h = brep
        .edges
        .keys()
        .find(|&eh| brep.edges[eh].curve.is_some())
        .unwrap();
    let he_h = brep.edges[edge_h].half_edge;
    let twin_h = brep.half_edges[he_h].twin;
    let twin_origin = brep.half_edges[twin_h].origin;

    // Remove the twin's origin vertex
    brep.vertices.remove(twin_origin);

    // Should not panic — just skip
    let _ = brep.validate_watertight();
}

#[test]
fn watertight_twin_missing_skips() {
    // Corrupt twin pointer to trigger he.twin None path
    let mut brep =
        crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

    let edge_h = brep
        .edges
        .keys()
        .find(|&eh| brep.edges[eh].curve.is_some())
        .unwrap();
    let he_h = brep.edges[edge_h].half_edge;
    // Make twin point to nonexistent handle
    brep.half_edges[he_h].twin = Default::default();

    // Should not panic — just skip
    let _ = brep.validate_watertight();
}

#[test]
fn watertight_gap_end_detected() {
    // Create a gap at the twin (end) vertex, not the start vertex.
    // This tests the gap_end > tolerance path (line 151-153).
    let mut brep =
        crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();

    // Find an edge with a curve, then move the twin's origin far away
    // while keeping the canonical half-edge's origin in place.
    for edge_h in brep.edges.keys().collect::<Vec<_>>() {
        if brep.edges[edge_h].curve.is_some() {
            let he_h = brep.edges[edge_h].half_edge;
            let twin_h = brep.half_edges[he_h].twin;
            let twin_origin = brep.half_edges[twin_h].origin;
            brep.vertices[twin_origin].position = Point3::new(999.0, 999.0, 999.0);
            break;
        }
    }

    let result = brep.validate_watertight();
    assert!(result.is_err());
}
