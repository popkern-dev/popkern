//! Shared trim endpoint agreement facts for BRep repair knowledge.

use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;
use crate::topology::{BRep, BrepLoopAdjacencyFacts, HalfEdgeHandle};

/// Shared facts describing geometric endpoint agreement across adjacent trims.
#[derive(Clone, Debug, PartialEq)]
pub struct BrepTrimEndpointAgreementFacts {
    /// First trim half-edge.
    pub half_edge: HalfEdgeHandle,
    /// Adjacent trim half-edge across the topological twin.
    pub adjacent_half_edge: HalfEdgeHandle,
    /// Agreement status.
    pub status: BrepTrimEndpointAgreementStatus,
    /// Max endpoint disagreement using the direct endpoint order.
    pub direct_disagreement: f64,
    /// Max endpoint disagreement using reversed adjacent endpoint order.
    pub reconciled_disagreement: f64,
    /// Best disagreement selected for downstream diagnostics.
    pub selected_disagreement: f64,
}

/// Status for trim endpoint agreement evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BrepTrimEndpointAgreementStatus {
    /// Direct endpoint order agrees within local BRep tolerance.
    AcceptedDirect,
    /// Reversed adjacent endpoint order agrees within local BRep tolerance.
    AcceptedReconciled,
    /// Both endpoint orderings exceed local BRep tolerance.
    SideDisagreement,
    /// Owning face has no surface.
    MissingSurface,
    /// Half-edge has no pcurve.
    MissingPcurve,
    /// Topology needed for evaluation is missing.
    MissingTopology,
    /// Surface or pcurve evaluation produced non-finite points.
    EvaluationFailed,
}

impl BrepTrimEndpointAgreementFacts {
    /// Migration-compatible compact diagnostic label (pre-render-key scheme).
    pub fn diagnostic_label(&self) -> String {
        format!(
            "{:?}<->{:?} status={:?} direct={:.3e} reconciled={:.3e} selected={:.3e}",
            self.half_edge,
            self.adjacent_half_edge,
            self.status,
            self.direct_disagreement,
            self.reconciled_disagreement,
            self.selected_disagreement
        )
    }
}

/// Collect trim endpoint agreement facts for adjacent loop pairs.
pub fn collect_brep_trim_endpoint_agreements(
    brep: &BRep,
    adjacencies: &[BrepLoopAdjacencyFacts],
) -> Vec<BrepTrimEndpointAgreementFacts> {
    collect_brep_trim_endpoint_agreements_with_tolerance(brep, adjacencies, brep.tolerance.linear())
}

/// Collect trim endpoint agreement facts using an explicit local tolerance.
pub fn collect_brep_trim_endpoint_agreements_with_tolerance(
    brep: &BRep,
    adjacencies: &[BrepLoopAdjacencyFacts],
    tolerance: f64,
) -> Vec<BrepTrimEndpointAgreementFacts> {
    adjacencies
        .iter()
        .map(|adjacency| trim_endpoint_agreement(brep, adjacency, tolerance))
        .collect()
}

fn trim_endpoint_agreement(
    brep: &BRep,
    adjacency: &BrepLoopAdjacencyFacts,
    tolerance: f64,
) -> BrepTrimEndpointAgreementFacts {
    let result = measure_trim_endpoint_agreement(brep, adjacency, tolerance);
    BrepTrimEndpointAgreementFacts {
        half_edge: adjacency.half_edge(),
        adjacent_half_edge: adjacency.adjacent_half_edge(),
        status: result.status,
        direct_disagreement: result.direct_disagreement,
        reconciled_disagreement: result.reconciled_disagreement,
        selected_disagreement: result
            .direct_disagreement
            .min(result.reconciled_disagreement),
    }
}

struct TrimEndpointAgreementMeasurement {
    status: BrepTrimEndpointAgreementStatus,
    direct_disagreement: f64,
    reconciled_disagreement: f64,
}

fn measure_trim_endpoint_agreement(
    brep: &BRep,
    adjacency: &BrepLoopAdjacencyFacts,
    tolerance: f64,
) -> TrimEndpointAgreementMeasurement {
    let points_a = match trim_endpoint_points(brep, adjacency.half_edge()) {
        Ok(points) => points,
        Err(status) => return unavailable(status),
    };
    let points_b = match trim_endpoint_points(brep, adjacency.adjacent_half_edge()) {
        Ok(points) => points,
        Err(status) => return unavailable(status),
    };
    let direct = max_pair_gap(points_a, points_b);
    let reconciled = max_pair_gap(points_a, [points_b[1], points_b[0]]);
    let status = if !direct.is_finite() || !reconciled.is_finite() {
        BrepTrimEndpointAgreementStatus::EvaluationFailed
    } else if direct <= tolerance {
        BrepTrimEndpointAgreementStatus::AcceptedDirect
    } else if reconciled <= tolerance {
        BrepTrimEndpointAgreementStatus::AcceptedReconciled
    } else {
        BrepTrimEndpointAgreementStatus::SideDisagreement
    };
    TrimEndpointAgreementMeasurement {
        status,
        direct_disagreement: direct,
        reconciled_disagreement: reconciled,
    }
}

fn trim_endpoint_points(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Result<[Point3; 2], BrepTrimEndpointAgreementStatus> {
    let half_edge_data = brep
        .half_edges
        .get(half_edge)
        .ok_or(BrepTrimEndpointAgreementStatus::MissingTopology)?;
    let loop_data = brep
        .loops
        .get(half_edge_data.loop_)
        .ok_or(BrepTrimEndpointAgreementStatus::MissingTopology)?;
    let surface = &brep
        .faces
        .get(loop_data.face)
        .ok_or(BrepTrimEndpointAgreementStatus::MissingTopology)?
        .surface
        .as_ref()
        .ok_or(BrepTrimEndpointAgreementStatus::MissingSurface)?
        .geometry;
    let pcurve = half_edge_data
        .pcurve
        .as_ref()
        .ok_or(BrepTrimEndpointAgreementStatus::MissingPcurve)?;
    let (t_start, t_end) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(t_start);
    let end_uv = pcurve.point_at(t_end);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    (point_is_finite(start) && point_is_finite(end))
        .then_some([start, end])
        .ok_or(BrepTrimEndpointAgreementStatus::EvaluationFailed)
}

fn max_pair_gap(a: [Point3; 2], b: [Point3; 2]) -> f64 {
    a[0].distance_to(&b[0]).max(a[1].distance_to(&b[1]))
}

fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn unavailable(status: BrepTrimEndpointAgreementStatus) -> TrimEndpointAgreementMeasurement {
    TrimEndpointAgreementMeasurement {
        status,
        direct_disagreement: f64::INFINITY,
        reconciled_disagreement: f64::INFINITY,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;

    #[test]
    fn empty_adjacency_list_collects_no_trim_endpoint_facts() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        assert!(collect_brep_trim_endpoint_agreements(&brep, &[]).is_empty());
    }

    #[test]
    fn trim_endpoint_label_contains_compact_fields() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let half_edge = brep.half_edges.keys().next().expect("cube has half-edge");
        let adjacent_half_edge = brep.half_edges[half_edge].twin;
        let facts = BrepTrimEndpointAgreementFacts {
            half_edge,
            adjacent_half_edge,
            status: BrepTrimEndpointAgreementStatus::AcceptedDirect,
            direct_disagreement: 0.0,
            reconciled_disagreement: 1.0,
            selected_disagreement: 0.0,
        };

        let label = facts.diagnostic_label();
        assert!(label.contains("status=AcceptedDirect"), "{label}");
        assert!(label.contains("direct=0.000e0"), "{label}");
        assert!(label.contains("selected=0.000e0"), "{label}");
    }

    #[test]
    fn migration_matrix_non_apex_trim_endpoint_fixture_reads_shared_fact() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let half_edge = brep.half_edges.keys().next().expect("cube has half-edge");
        let half_edge_data = &brep.half_edges[half_edge];
        let adjacent_half_edge = half_edge_data.twin;
        let adjacent_data = &brep.half_edges[adjacent_half_edge];
        let adjacency = BrepLoopAdjacencyFacts {
            face: brep.loops[half_edge_data.loop_].face,
            loop_handle: half_edge_data.loop_,
            half_edge,
            adjacent_face: brep.loops[adjacent_data.loop_].face,
            adjacent_loop: adjacent_data.loop_,
            adjacent_half_edge,
        };

        let facts = collect_brep_trim_endpoint_agreements(&brep, &[adjacency]);

        assert_eq!(facts.len(), 1);
        assert_eq!(facts[0].half_edge, half_edge);
        assert_eq!(facts[0].adjacent_half_edge, adjacent_half_edge);
        assert!(
            facts[0].selected_disagreement.is_finite()
                || matches!(
                    facts[0].status,
                    BrepTrimEndpointAgreementStatus::MissingSurface
                        | BrepTrimEndpointAgreementStatus::MissingPcurve
                        | BrepTrimEndpointAgreementStatus::MissingTopology
                        | BrepTrimEndpointAgreementStatus::EvaluationFailed
                )
        );
    }
}
