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

/// PK-287 direction contract: a pcurve tracing the right point set in
/// the WRONG direction must be flagged by design, not by the sampling
/// floor that accidentally caught it before (floor scaled with edge
/// length, so short reversed pcurves were blessed and long ones
/// flagged — incoherent at model scale).
#[test]
fn reversed_pcurve_is_flagged_by_direction_contract() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let heh = brep.half_edges.keys().next().unwrap();
    let pcurve = brep.half_edges[heh].pcurve.as_mut().unwrap();
    let nurbs = pcurve.as_nurbs_mut().expect("cube pcurves are NURBS");
    nurbs.control_points.reverse();
    nurbs.weights.reverse();
    let report = brep.measure_half_edge_pcurve_deviation(heh).unwrap();
    assert!(
        !report.trace_direction_agrees,
        "reversed pcurve must fail the direction contract"
    );
    // The trace itself is coincident — only the direction gate keeps
    // the deviation honest.
    assert!(
        report.max_trace_deviation < 1e-9,
        "reversed trace is still geometrically coincident, measured {:.3e}",
        report.max_trace_deviation
    );
    assert!(
        report.max_deviation > 0.1,
        "accepted deviation must reflect the reversal, measured {:.3e}",
        report.max_deviation
    );
}

/// Direction classifier: forward-monotone sequences pass, a single
/// seam wrap on a closed trace passes, mid-band backward steps (a
/// reversed traversal) fail, and degenerate ranges pass trivially.
#[test]
fn trace_direction_classifier_bands() {
    use super::trace_parameters_advance_forward;
    let range = (0.0, 1.0);
    let forward = [0.0, 0.25, 0.5, 0.75, 1.0];
    assert!(trace_parameters_advance_forward(
        forward.iter().copied(),
        range,
        4
    ));
    // Seam wrap: nearest parameter jumps from the high end back to
    // the low end once (closed trace with offset sampling).
    let wrapped = [0.6, 0.85, 0.1, 0.35, 0.6];
    assert!(trace_parameters_advance_forward(
        wrapped.iter().copied(),
        range,
        4
    ));
    // Reversed traversal: every step backward by ~range/samples.
    let reversed = [1.0, 0.75, 0.5, 0.25, 0.0];
    assert!(!trace_parameters_advance_forward(
        reversed.iter().copied(),
        range,
        4
    ));
    // Noise-scale jitter is not a violation.
    let jitter = [0.0, 0.25, 0.249, 0.5, 1.0];
    assert!(trace_parameters_advance_forward(
        jitter.iter().copied(),
        range,
        4
    ));
    // Degenerate range: trivially forward.
    assert!(trace_parameters_advance_forward(
        forward.iter().copied(),
        (0.5, 0.5),
        4
    ));
}
