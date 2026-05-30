//! Candidate-space cap/support continuity measurement for multipatch routing.

use popkern_core::geometry::surface::SurfaceEval;

use crate::types::ApexConstructionDiagnostic;

use super::super::multipatch::ExactBoundaryMultiPatchCandidate;
use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::sector::ExactBoundaryCapSectorInput;

pub(super) fn candidate_selected_strip_continuity_diagnostic(
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
    pcurves: &PreparedExactBoundaryPcurves,
) -> Option<ApexConstructionDiagnostic> {
    let pair = pcurves.pairs.iter().find(|pair| {
        pair.sector_name == sector.sector_name && pair.side == candidate.strip_inventory.side
    })?;
    let evidence = measure_candidate_selected_strip_continuity(candidate, pair);
    if evidence.worst_g0_gap <= 5.0e-3 && evidence.worst_g1_coplanarity <= 5.0e-2 {
        return None;
    }
    let target = crate::types::CapContinuityTarget::native_cap_g1(5.0e-3, 5.0e-2);
    let continuity = if evidence.worst_g0_gap > 5.0e-3 {
        crate::types::CapContinuityReport::declined_g0(
            target,
            1,
            evidence.worst_g0_gap,
            evidence.worst_g1_coplanarity,
        )
    } else {
        crate::types::CapContinuityReport::declined_g1(
            target,
            1,
            evidence.worst_g0_gap,
            evidence.worst_g1_coplanarity,
        )
    };
    let sample_parameter = if evidence.worst_g0_gap > 5.0e-3 {
        evidence.worst_g0_parameter
    } else {
        evidence.worst_g1_parameter
    };
    Some(ApexConstructionDiagnostic::CapSupportContinuity(
        crate::types::CapSupportContinuityDiagnostic {
            continuity,
            cap_face: sector.cap_face,
            support_face: pair.support_face,
            edge: pair.edge,
            measurement_scope: "multipatch_candidate_selected_strip".to_string(),
            sample_parameter: Some(sample_parameter),
            worst_g1_coplanarity: evidence.worst_g1_coplanarity,
            gate: 5.0e-2,
            edges_measured: 1,
        },
    ))
}

#[derive(Debug, Clone, Copy)]
struct CandidateStripContinuityEvidence {
    worst_g0_gap: f64,
    worst_g0_parameter: f64,
    worst_g1_coplanarity: f64,
    worst_g1_parameter: f64,
}

fn measure_candidate_selected_strip_continuity(
    candidate: &ExactBoundaryMultiPatchCandidate,
    pair: &PreparedExactBoundaryPcurvePair,
) -> CandidateStripContinuityEvidence {
    let cap_pcurve = &candidate.strip_candidate.outer_pcurve;
    let support_pcurve = &pair.support_pcurve;
    let cap_forward_point = candidate_strip_point(candidate, cap_pcurve.t_start);
    let cap_end_point = candidate_strip_point(candidate, cap_pcurve.t_end);
    let support_start = support_pcurve_point(pair, support_pcurve.t_start);
    let support_end = support_pcurve_point(pair, support_pcurve.t_end);
    let support_forward = cap_forward_point.distance_to(&support_start)
        + cap_end_point.distance_to(&support_end)
        <= cap_forward_point.distance_to(&support_end) + cap_end_point.distance_to(&support_start);
    let mut evidence = CandidateStripContinuityEvidence {
        worst_g0_gap: 0.0,
        worst_g0_parameter: 0.0,
        worst_g1_coplanarity: 0.0,
        worst_g1_parameter: 0.0,
    };
    for index in 0..100 {
        measure_sample(candidate, pair, support_forward, index, &mut evidence);
    }
    evidence
}

fn measure_sample(
    candidate: &ExactBoundaryMultiPatchCandidate,
    pair: &PreparedExactBoundaryPcurvePair,
    support_forward: bool,
    index: usize,
    evidence: &mut CandidateStripContinuityEvidence,
) {
    let cap_pcurve = &candidate.strip_candidate.outer_pcurve;
    let support_pcurve = &pair.support_pcurve;
    let parameter = 0.05 + 0.90 * index as f64 / 99.0;
    let support_parameter = if support_forward {
        parameter
    } else {
        1.0 - parameter
    };
    let cap_t = cap_pcurve.t_start + (cap_pcurve.t_end - cap_pcurve.t_start) * parameter;
    let support_t = support_pcurve.t_start
        + (support_pcurve.t_end - support_pcurve.t_start) * support_parameter;
    let cap_uv = cap_pcurve.point_at(cap_t);
    let support_uv = support_pcurve.point_at(support_t);
    let cap_point = candidate
        .strip_candidate
        .surface
        .point_at(cap_uv.x, cap_uv.y);
    let support_point = pair.support_surface.point_at(support_uv.x, support_uv.y);
    let g0_gap = cap_point.distance_to(&support_point);
    if g0_gap > evidence.worst_g0_gap {
        evidence.worst_g0_gap = g0_gap;
        evidence.worst_g0_parameter = parameter;
    }
    let cap_normal = candidate
        .strip_candidate
        .surface
        .partial_u(cap_uv.x, cap_uv.y)
        .cross(
            &candidate
                .strip_candidate
                .surface
                .partial_v(cap_uv.x, cap_uv.y),
        );
    let support_normal = pair
        .support_surface
        .partial_u(support_uv.x, support_uv.y)
        .cross(&pair.support_surface.partial_v(support_uv.x, support_uv.y));
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length > 1.0e-300 && support_length > 1.0e-300 {
        let dot = cap_normal.x * support_normal.x
            + cap_normal.y * support_normal.y
            + cap_normal.z * support_normal.z;
        let g1_coplanarity = 1.0 - (dot / (cap_length * support_length)).abs();
        if g1_coplanarity > evidence.worst_g1_coplanarity {
            evidence.worst_g1_coplanarity = g1_coplanarity;
            evidence.worst_g1_parameter = parameter;
        }
    }
}

fn candidate_strip_point(
    candidate: &ExactBoundaryMultiPatchCandidate,
    parameter: f64,
) -> popkern_core::geometry::Point3 {
    let uv = candidate.strip_candidate.outer_pcurve.point_at(parameter);
    candidate.strip_candidate.surface.point_at(uv.x, uv.y)
}

fn support_pcurve_point(
    pair: &PreparedExactBoundaryPcurvePair,
    parameter: f64,
) -> popkern_core::geometry::Point3 {
    let uv = pair.support_pcurve.point_at(parameter);
    pair.support_surface.point_at(uv.x, uv.y)
}
