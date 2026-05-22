//! Pcurve-trace repair candidate for edge/surface compatibility.

use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairDecision, BrepEdgeSurfaceCompatibilityRepairOutcome,
};
use crate::topology::brep_edge_surface_compatibility_diagnostics::BrepEdgeSurfaceRepairEvidence;
use crate::topology::brep_edge_surface_compatibility_helpers::residuals_for_candidate;
use crate::topology::{
    BRep, EdgeHandle, HalfEdgeHandle, SharedTrimEdgeInput, commit_shared_trim_edge_transactionally,
};

#[allow(clippy::too_many_arguments)]
pub(super) fn try_commit_pcurve_trace_edge(
    brep: &mut BRep,
    edge: EdgeHandle,
    he_a: HalfEdgeHandle,
    he_b: HalfEdgeHandle,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    options: BrepEdgeSurfaceCompatibilityOptions,
    mut evidence: BrepEdgeSurfaceRepairEvidence,
) -> Option<BrepEdgeSurfaceCompatibilityRepairDecision> {
    let candidates = [
        pcurve_trace_curve(brep, he_a, surf_a),
        pcurve_trace_curve(brep, he_b, surf_b),
    ];
    let (curve, residual_a, residual_b) = candidates
        .into_iter()
        .flatten()
        .filter_map(|curve| {
            let (residual_a, residual_b) = residuals_for_candidate(
                &curve.geometry,
                (curve.t_start, curve.t_end),
                surf_a,
                surf_b,
                options.samples,
            );
            (residual_a.max(residual_b) <= options.tolerance)
                .then_some((curve, residual_a, residual_b))
        })
        .min_by(|left, right| left.1.max(left.2).total_cmp(&right.1.max(right.2)))?;
    let pcurve_a = brep.half_edges.get(he_a)?.pcurve.clone()?;
    let pcurve_b = brep.half_edges.get(he_b)?.pcurve.clone()?;
    evidence.residual_deviation_a = residual_a;
    evidence.residual_deviation_b = residual_b;
    let rollback_snapshot = brep.clone();
    if let Err(failure) = commit_shared_trim_edge_transactionally(
        brep,
        SharedTrimEdgeInput {
            edge,
            half_edge_a: he_a,
            half_edge_b: he_b,
            curve,
            pcurve_a,
            pcurve_b,
            tolerance: options.tolerance,
            validate_brep: options.validate_output,
            provenance: "edge_surface_compatibility_pcurve_trace_repair",
        },
    ) {
        *brep = rollback_snapshot;
        evidence.shared_trim_failure = Some(failure.reason);
        evidence.shared_trim_curve_surface_deviation =
            failure.validation.max_curve_surface_deviation;
        evidence.shared_trim_side_agreement_deviation =
            failure.validation.max_side_agreement_deviation;
        evidence.shared_trim_endpoint_gap = failure.validation.max_endpoint_gap;
        return Some(BrepEdgeSurfaceCompatibilityRepairDecision {
            edge,
            outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(
                BrepEdgeSurfaceCompatibilityDeclineReason::SharedTrimValidationFailed,
            ),
            endpoint_gap: 0.0,
            residual_deviation: residual_a.max(residual_b),
            evidence,
            shared_fact_ids: Vec::new(),
        });
    }
    Some(BrepEdgeSurfaceCompatibilityRepairDecision {
        edge,
        outcome: BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary,
        endpoint_gap: 0.0,
        residual_deviation: residual_a.max(residual_b),
        evidence,
        shared_fact_ids: Vec::new(),
    })
}

fn pcurve_trace_curve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<CurveOnEdge> {
    let pcurve = brep.half_edges.get(half_edge)?.pcurve.as_ref()?;
    let (t_start, t_end) = pcurve.parameter_range();
    let samples = 17usize;
    let mut points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let parameter = t_start + (t_end - t_start) * fraction;
        let uv = pcurve.point_at(parameter);
        points.push(surface.point_at(uv.x, uv.y));
        parameters.push(parameter);
    }
    let degree = 3.min(points.len().saturating_sub(1));
    let geometry = crate::geometry::fitting::interpolate_curve(&points, &parameters, degree)
        .ok()
        .map(CurveGeometry::Nurbs)?;
    Some(CurveOnEdge {
        geometry,
        t_start,
        t_end,
    })
}
