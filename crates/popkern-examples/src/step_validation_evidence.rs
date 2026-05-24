use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, BrepExportQualityOptions,
    analyze_brep_edge_surface_compatibility,
};

use super::{
    StepRoundTripEdgePcurveEvidence, StepRoundTripEdgeSurfaceEvidence,
    StepRoundTripHalfEdgePcurveEvidence, StepRoundTripMissingEdgeCurveEvidence,
    StepRoundTripMissingEdgeCurveHalfEdgeEvidence,
};

pub(super) fn missing_edge_curve_evidence(
    body: &BRep,
) -> Vec<StepRoundTripMissingEdgeCurveEvidence> {
    let mut evidence: Vec<_> = body
        .edges
        .iter()
        .filter(|(_, edge)| edge.curve.is_none())
        .map(|(edge, _)| {
            let mut half_edges: Vec<_> = body
                .half_edges
                .iter()
                .filter(|(_, half_edge)| half_edge.edge == edge)
                .map(|(half_edge, half_edge_data)| {
                    let face = body
                        .loops
                        .get(half_edge_data.loop_)
                        .map(|loop_data| loop_data.face)
                        .unwrap_or_default();
                    let surface_present = body
                        .faces
                        .get(face)
                        .is_some_and(|face_data| face_data.surface.is_some());
                    StepRoundTripMissingEdgeCurveHalfEdgeEvidence {
                        half_edge,
                        face,
                        pcurve_present: half_edge_data.pcurve.is_some(),
                        surface_present,
                    }
                })
                .collect();
            half_edges.sort_by_key(|half_edge| format!("{:?}", half_edge.half_edge));
            StepRoundTripMissingEdgeCurveEvidence { edge, half_edges }
        })
        .collect();
    evidence.sort_by_key(|edge| format!("{:?}", edge.edge));
    evidence
}

pub(super) fn pcurve_evidence(
    body: &BRep,
    tolerance: Option<f64>,
) -> Vec<StepRoundTripEdgePcurveEvidence> {
    let mut evidence: Vec<_> = body
        .edges
        .keys()
        .filter_map(|edge| body.measure_edge_pcurve_deviation(edge))
        .filter(|report| {
            tolerance
                .map(|tolerance| report.max_deviation > tolerance)
                .unwrap_or(report.max_deviation > 0.0)
        })
        .map(|report| StepRoundTripEdgePcurveEvidence {
            edge: report.edge,
            max_deviation: report.max_deviation,
            half_edges: report
                .half_edges
                .into_iter()
                .map(|half_edge| StepRoundTripHalfEdgePcurveEvidence {
                    half_edge: half_edge.half_edge,
                    face: half_edge.face,
                    max_deviation: half_edge.max_deviation,
                    same_fraction_deviation: half_edge.max_same_fraction_deviation,
                    trace_deviation: half_edge.max_trace_deviation,
                    worst_fraction: half_edge.worst_fraction,
                })
                .collect(),
        })
        .collect();
    evidence.sort_by(|a, b| b.max_deviation.total_cmp(&a.max_deviation));
    evidence
}

pub(super) fn edge_surface_evidence(
    body: &BRep,
    quality: BrepExportQualityOptions,
) -> Vec<StepRoundTripEdgeSurfaceEvidence> {
    let report = analyze_brep_edge_surface_compatibility(
        body,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: quality.max_edge_surface_deviation.unwrap_or(f64::INFINITY),
            max_anchor_motion: 0.0,
            samples: quality.trim.samples,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    let mut evidence: Vec<_> = report
        .findings
        .into_iter()
        .map(|finding| StepRoundTripEdgeSurfaceEvidence {
            edge: finding.edge,
            half_edge: finding.half_edge,
            face: finding.face,
            max_deviation: finding.max_deviation,
            max_curve_parameter: finding.max_curve_parameter,
        })
        .collect();
    evidence.sort_by(|a, b| b.max_deviation.total_cmp(&a.max_deviation));
    evidence
}
