//! Tests for `popkern_core::diagnostics`.

use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn unit_cube_has_no_defects() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = Diagnostics::run(&brep);
    assert_eq!(report.topology.faces, 6);
    assert_eq!(report.topology.vertices, 8);
    assert_eq!(report.topology.edges, 12);
    assert_eq!(report.topology.euler_chi, 2);
    assert!(report.defects.naked_edges.is_empty());
    assert!(report.defects.degenerate_edges.is_empty());
    assert!(report.defects.coincident_vertices.is_empty());
    assert!(report.defects.t_junctions.is_empty());
    assert!(report.defects.nurbs_defects.is_empty());
    assert!(report.defects.pcurve_escapes.is_empty());
    assert!(report.defects.face_overlaps.is_empty());
    assert_eq!(report.defects.tess_non_manifold_count, 0);
    // Cube AABB is [0,0,0] → [1,1,1] with diagonal √3.
    assert!((report.geometry.bbox_diagonal - 3.0_f64.sqrt()).abs() < 1e-9);
}

#[test]
fn surface_kinds_histogram_counts_planes() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = Diagnostics::run(&brep);
    assert_eq!(report.topology.surface_kinds.get("Plane"), Some(&6));
}

#[test]
fn coincident_vertices_detected() {
    let mut brep = BRep::new(Tolerance::default());
    let _v1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let v2 = brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(1e-12, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let report = Diagnostics::run(&brep);
    assert!(!report.defects.coincident_vertices.is_empty());
    let pair = report.defects.coincident_vertices[0];
    assert!(pair.a == v2 || pair.b == v2);
}

#[test]
fn nurbs_invalid_weight_detected() {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    use crate::geometry::nurbs_curve::NurbsCurve;

    // Build a degree-1 NURBS line with one zero weight — should
    // surface as InvalidWeight.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let edge = brep.half_edges[mev.half_edge].edge;
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            // Zero weight — invalid.
            weights: vec![1.0, 0.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    });
    let report = Diagnostics::run(&brep);
    assert!(
        report
            .defects
            .nurbs_defects
            .iter()
            .any(|d| matches!(d.kind, NurbsKind::InvalidWeight { .. })),
        "expected InvalidWeight defect, got {:?}",
        report.defects.nurbs_defects
    );
}

#[test]
fn nurbs_cp_knot_mismatch_detected() {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    use crate::geometry::nurbs_curve::NurbsCurve;

    // n+degree+1 should equal knots.len() — break it.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let edge = brep.half_edges[mev.half_edge].edge;
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            // 2 CPs + degree 1 → expect 4 knots; supply 5.
            knots: vec![0.0, 0.0, 0.5, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    });
    let report = Diagnostics::run(&brep);
    assert!(
        report
            .defects
            .nurbs_defects
            .iter()
            .any(|d| matches!(d.kind, NurbsKind::CpKnotMismatch { .. }))
    );
}

#[test]
fn t_junction_detected_when_vertex_lies_on_edge_interior() {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};

    // Build an edge from (0,0,0) to (2,0,0), then a stray vertex at
    // (1,0,0) — exactly on the edge midpoint.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(2.0, 0.0, 0.0))
        .unwrap();
    let edge = brep.half_edges[mev.half_edge].edge;
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: crate::geometry::Vec3::new(1.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 2.0,
    });
    // Stray vertex on the edge interior (no edges/loops referencing it).
    brep.vertices.insert(crate::topology::Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
        tolerance: None,
        name: None,
    });
    let report = Diagnostics::run(&brep);
    assert!(
        !report.defects.t_junctions.is_empty(),
        "expected T-junction detection"
    );
}
