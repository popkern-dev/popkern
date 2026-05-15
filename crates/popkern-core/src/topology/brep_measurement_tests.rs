//! Tests for BRep geometric deviation measurement.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn edge_deviation_on_unit_cube_is_small() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    for (eh, _) in &brep.edges {
        let dev = brep.measure_edge_deviation(eh);
        // Cube edges are linear — pcurves should match exactly.
        assert!(
            dev < 1e-6,
            "Edge deviation {dev} exceeds threshold for cube edge"
        );
    }
}

#[test]
fn edge_pcurve_deviation_report_identifies_adjacent_faces() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    for (eh, edge) in &brep.edges {
        let report = brep
            .measure_edge_pcurve_deviation(eh)
            .expect("cube edge has a 3D curve");
        assert_eq!(report.edge, eh);
        assert!(
            report.max_deviation < 1e-6,
            "Edge report deviation {} exceeds threshold for cube edge",
            report.max_deviation
        );
        assert_eq!(
            report.half_edges.len(),
            2,
            "closed cube edge should report both adjacent half-edges"
        );
        for half_edge_report in report.half_edges {
            assert_eq!(
                half_edge_report.edge_parameter_range,
                if half_edge_report.half_edge == edge.half_edge {
                    (
                        edge.curve.as_ref().unwrap().t_start,
                        edge.curve.as_ref().unwrap().t_end,
                    )
                } else {
                    (
                        edge.curve.as_ref().unwrap().t_end,
                        edge.curve.as_ref().unwrap().t_start,
                    )
                }
            );
            assert!(
                brep.faces.contains_key(half_edge_report.face),
                "reported face must belong to the BRep"
            );
            assert!(
                (0.0..=1.0).contains(&half_edge_report.worst_fraction),
                "worst sample fraction must be normalized"
            );
        }
    }
}

#[test]
fn vertex_deviation_on_unit_cube_is_small() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    for (vh, _) in &brep.vertices {
        let dev = brep.measure_vertex_deviation(vh);
        // Cube vertices should sit exactly on their edge curves.
        assert!(
            dev < 1e-10,
            "Vertex deviation {dev} exceeds threshold for cube vertex"
        );
    }
}

#[test]
fn edge_deviation_returns_zero_for_invalid_handle() {
    let brep = BRep::new(Tolerance::default());
    // Use a default (invalid) handle.
    let eh = slotmap::KeyData::from_ffi(u64::MAX).into();
    assert_eq!(brep.measure_edge_deviation(eh), 0.0);
    assert!(brep.measure_edge_pcurve_deviation(eh).is_none());
}

#[test]
fn vertex_deviation_returns_zero_for_invalid_handle() {
    let brep = BRep::new(Tolerance::default());
    let vh = slotmap::KeyData::from_ffi(u64::MAX).into();
    assert_eq!(brep.measure_vertex_deviation(vh), 0.0);
}
