use super::*;
use popkern_core::{
    geometry::{Point3, Vec3},
    topology::{BRep, Vertex},
};

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn merge_removes_coincident_vertices() {
    let mut brep = BRep::new(tol());
    let _v0 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let _v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let _v2 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });

    assert_eq!(brep.vertices.len(), 3);
    let removed = merge_coincident_vertices(&mut brep, &tol());
    assert_eq!(removed, 1);
    assert_eq!(brep.vertices.len(), 2);
}

#[test]
fn merge_noop_when_no_duplicates() {
    let mut brep = BRep::new(tol());
    let _v0 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let _v1 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });

    let removed = merge_coincident_vertices(&mut brep, &tol());
    assert_eq!(removed, 0);
    assert_eq!(brep.vertices.len(), 2);
}

#[test]
fn sew_gaps_pass_name() {
    let pass = SewGapsPass;
    assert_eq!(pass.name(), "sew_gaps");
}

#[test]
fn sew_gaps_applies_to_validate_phase() {
    let pass = SewGapsPass;
    let entry = DiagnosticEntry {
        severity: crate::step::diagnostics::Severity::Warning,
        phase: Phase::Validate,
        entity_id: None,
        message: "test".to_string(),
    };
    assert!(pass.applies_to(&entry));
}

#[test]
fn sew_gaps_does_not_apply_to_parse_phase() {
    let pass = SewGapsPass;
    let entry = DiagnosticEntry {
        severity: crate::step::diagnostics::Severity::Warning,
        phase: Phase::Parse,
        entity_id: None,
        message: "test".to_string(),
    };
    assert!(!pass.applies_to(&entry));
}

#[test]
fn reconcile_edge_curves_empty_brep() {
    let mut brep = BRep::new(tol());
    let count = reconcile_edge_curves(&mut brep);
    assert_eq!(count, 0);
}

#[test]
fn sew_gaps_attempt_repair_no_changes() {
    let pass = SewGapsPass;
    let mut brep = BRep::new(tol());
    let _ = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let entry = DiagnosticEntry {
        severity: crate::step::diagnostics::Severity::Warning,
        phase: Phase::Validate,
        entity_id: None,
        message: "gap".to_string(),
    };
    let mut log = DiagnosticLog::new();
    let result = pass.attempt_repair(&mut brep, &entry, &mut log);
    assert_eq!(result.ok(), Some(false));
}

#[test]
fn refit_line_curve_degenerate_returns_none() {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(0.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    let result = refit_line_curve(
        &curve,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        &tol(),
    );
    assert!(result.is_none());
}

#[test]
fn refit_line_curve_non_line_returns_none() {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    let result = refit_line_curve(
        &curve,
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        &tol(),
    );
    assert!(result.is_none());
}

#[test]
fn refit_line_curve_successful() {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    // Line from (0,0,0) in direction (1,0,0)
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    // Refit to start=(2,0,0) end=(5,0,0) -- these lie on the line
    let result = refit_line_curve(
        &curve,
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(5.0, 0.0, 0.0),
        &tol(),
    );
    assert!(result.is_some());
    let new_curve = result.unwrap();
    assert!((new_curve.t_start - 2.0).abs() < 1e-10);
    assert!((new_curve.t_end - 5.0).abs() < 1e-10);
}

#[test]
fn sew_gaps_attempt_repair_with_merges() {
    let pass = SewGapsPass;
    let mut brep = BRep::new(tol());
    // Insert two coincident vertices
    let _ = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let _ = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let entry = DiagnosticEntry {
        severity: crate::step::diagnostics::Severity::Warning,
        phase: Phase::Validate,
        entity_id: None,
        message: "gap".to_string(),
    };
    let mut log = DiagnosticLog::new();
    let result = pass.attempt_repair(&mut brep, &entry, &mut log);
    assert_eq!(result.ok(), Some(true));
    assert_eq!(brep.vertices.len(), 1);
}

#[test]
fn sew_gaps_applies_to_map_topology_phase() {
    let pass = SewGapsPass;
    let entry = DiagnosticEntry {
        severity: crate::step::diagnostics::Severity::Warning,
        phase: Phase::MapTopology,
        entity_id: None,
        message: "test".to_string(),
    };
    assert!(pass.applies_to(&entry));
}

#[test]
fn merge_rewires_half_edge_origins() {
    use popkern_core::topology::{Edge, HalfEdge, HalfEdgeHandle};
    let mut brep = BRep::new(tol());
    let v0 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    // v1 is coincident with v0
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });

    let edge_handle = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: None,
        tolerance: None,
        name: None,
        canonical: Default::default(),
    });
    // half-edge from v1 (will be remapped to v0)
    let he = brep.half_edges.insert(HalfEdge {
        origin: v1,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: edge_handle,
        loop_: Default::default(),
        pcurve: None,
    });
    brep.half_edges[he].twin = he; // self-twin for simplicity
    brep.edges[edge_handle].half_edge = he;
    brep.vertices[v1].outgoing = he;

    let removed = merge_coincident_vertices(&mut brep, &tol());
    assert_eq!(removed, 1);
    // The half-edge origin should now point to v0
    assert_eq!(brep.half_edges[he].origin, v0);
    assert_eq!(brep.vertices.len(), 2);
}
