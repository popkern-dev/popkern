//! Transactional construction primitive for shared trim-edge geometry.
//!
//! A shared trim edge is valid only when the topological edge curve and both
//! adjacent half-edge pcurves describe the same 3-D locus on their owning
//! surfaces. Construction clients use this module to commit those three
//! geometric representations atomically after validation.

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

/// Input for a construction-time shared trim edge commit.
#[derive(Debug, Clone)]
pub struct SharedTrimEdgeInput {
    /// Topological edge to update.
    pub edge: EdgeHandle,
    /// First half-edge side to update.
    pub half_edge_a: HalfEdgeHandle,
    /// Second half-edge side to update. Must be the twin of `half_edge_a`.
    pub half_edge_b: HalfEdgeHandle,
    /// Authoritative 3-D curve for the shared edge.
    pub curve: CurveOnEdge,
    /// Pcurve for `half_edge_a`, parameterized in that half-edge's traversal direction.
    pub pcurve_a: Pcurve,
    /// Pcurve for `half_edge_b`, parameterized in that half-edge's traversal direction.
    pub pcurve_b: Pcurve,
    /// Local construction tolerance for endpoint and surface agreement.
    pub tolerance: f64,
    /// Whether the whole BRep must pass validation after the local commit.
    pub validate_brep: bool,
    /// Human-readable provenance label included in diagnostics.
    pub provenance: &'static str,
}

/// Result of a successful shared trim edge commit.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedTrimEdgeCommit {
    /// Updated edge.
    pub edge: EdgeHandle,
    /// Maximum sampled curve-to-pcurve/surface deviation.
    pub max_curve_surface_deviation: f64,
    /// Maximum sampled gap between side A and side B surface images.
    pub max_side_agreement_deviation: f64,
    /// Maximum endpoint gap against topology vertices.
    pub max_endpoint_gap: f64,
}

/// Machine-readable validation/commit failure.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedTrimEdgeFailure {
    /// Decline reason.
    pub reason: SharedTrimEdgeFailureReason,
    /// Provenance label from the input.
    pub provenance: &'static str,
    /// Best available validation summary.
    pub validation: SharedTrimEdgeValidation,
}

/// Reason a shared trim edge candidate was declined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SharedTrimEdgeFailureReason {
    /// Edge, half-edge, loop, face, vertex, or surface topology was missing.
    MissingTopology,
    /// The supplied half-edges are not the two sides of the supplied edge.
    HalfEdgeEdgeMismatch,
    /// The supplied half-edges are not twins.
    HalfEdgesNotTwins,
    /// Curve endpoints do not agree with edge topology.
    EndpointMismatch,
    /// One or both pcurves do not agree with the 3-D curve on their surfaces.
    CurveSurfaceMismatch,
    /// The two pcurve/surface images do not describe the same 3-D locus.
    SideAgreementMismatch,
    /// Repaired clone failed BRep validation.
    ValidationFailed,
}

/// Read-only validation summary for a shared trim edge candidate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SharedTrimEdgeValidation {
    /// Maximum sampled curve-to-pcurve/surface deviation.
    pub max_curve_surface_deviation: f64,
    /// Maximum sampled side-A to side-B pcurve/surface deviation.
    pub max_side_agreement_deviation: f64,
    /// Maximum endpoint gap against topology vertices.
    pub max_endpoint_gap: f64,
}

impl Default for SharedTrimEdgeValidation {
    fn default() -> Self {
        Self {
            max_curve_surface_deviation: f64::INFINITY,
            max_side_agreement_deviation: f64::INFINITY,
            max_endpoint_gap: f64::INFINITY,
        }
    }
}

/// Validate and atomically commit a shared trim edge.
///
/// On failure the input BRep is restored to its pre-call state.
pub fn commit_shared_trim_edge_transactionally(
    brep: &mut BRep,
    input: SharedTrimEdgeInput,
) -> Result<SharedTrimEdgeCommit, SharedTrimEdgeFailure> {
    let snapshot = brep.clone();
    match validate_shared_trim_edge(brep, &input) {
        Ok(validation) => {
            brep.edges[input.edge].curve = Some(input.curve.clone());
            brep.half_edges[input.half_edge_a].pcurve = Some(input.pcurve_a.clone());
            brep.half_edges[input.half_edge_b].pcurve = Some(input.pcurve_b.clone());
            if input.validate_brep
                && let Err(_err) = brep.validate()
            {
                *brep = snapshot;
                return Err(failure(
                    SharedTrimEdgeFailureReason::ValidationFailed,
                    input.provenance,
                    validation,
                ));
            }
            Ok(SharedTrimEdgeCommit {
                edge: input.edge,
                max_curve_surface_deviation: validation.max_curve_surface_deviation,
                max_side_agreement_deviation: validation.max_side_agreement_deviation,
                max_endpoint_gap: validation.max_endpoint_gap,
            })
        }
        Err(err) => {
            *brep = snapshot;
            Err(err)
        }
    }
}

/// Validate a shared trim edge candidate without mutating the BRep.
pub fn validate_shared_trim_edge(
    brep: &BRep,
    input: &SharedTrimEdgeInput,
) -> Result<SharedTrimEdgeValidation, SharedTrimEdgeFailure> {
    let Some(edge_data) = brep.edges.get(input.edge) else {
        return Err(failure(
            SharedTrimEdgeFailureReason::MissingTopology,
            input.provenance,
            SharedTrimEdgeValidation::default(),
        ));
    };
    let Some(he_a) = brep.half_edges.get(input.half_edge_a) else {
        return Err(missing(input.provenance));
    };
    let Some(he_b) = brep.half_edges.get(input.half_edge_b) else {
        return Err(missing(input.provenance));
    };
    if he_a.edge != input.edge || he_b.edge != input.edge {
        return Err(failure(
            SharedTrimEdgeFailureReason::HalfEdgeEdgeMismatch,
            input.provenance,
            SharedTrimEdgeValidation::default(),
        ));
    }
    if he_a.twin != input.half_edge_b || he_b.twin != input.half_edge_a {
        return Err(failure(
            SharedTrimEdgeFailureReason::HalfEdgesNotTwins,
            input.provenance,
            SharedTrimEdgeValidation::default(),
        ));
    }
    let Some(surface_a) = surface_for_half_edge(brep, input.half_edge_a) else {
        return Err(missing(input.provenance));
    };
    let Some(surface_b) = surface_for_half_edge(brep, input.half_edge_b) else {
        return Err(missing(input.provenance));
    };
    let Some((start, end)) = canonical_edge_vertices(brep, edge_data.half_edge) else {
        return Err(missing(input.provenance));
    };

    let endpoint_gap = edge_endpoint_gap(&input.curve, start, end);
    let mut validation = SharedTrimEdgeValidation {
        max_endpoint_gap: endpoint_gap,
        ..SharedTrimEdgeValidation::default()
    };
    if endpoint_gap > input.tolerance {
        return Err(failure(
            SharedTrimEdgeFailureReason::EndpointMismatch,
            input.provenance,
            validation,
        ));
    }

    let forward_a = input.half_edge_a == edge_data.half_edge;
    let forward_b = input.half_edge_b == edge_data.half_edge;
    let curve_surface_a =
        curve_surface_deviation(&input.curve, &input.pcurve_a, surface_a, forward_a);
    let curve_surface_b =
        curve_surface_deviation(&input.curve, &input.pcurve_b, surface_b, forward_b);
    validation.max_curve_surface_deviation = curve_surface_a.max(curve_surface_b);
    validation.max_side_agreement_deviation =
        side_agreement_deviation(&input.pcurve_a, surface_a, &input.pcurve_b, surface_b);
    if validation.max_curve_surface_deviation > input.tolerance {
        return Err(failure(
            SharedTrimEdgeFailureReason::CurveSurfaceMismatch,
            input.provenance,
            validation,
        ));
    }
    if validation.max_side_agreement_deviation > input.tolerance {
        return Err(failure(
            SharedTrimEdgeFailureReason::SideAgreementMismatch,
            input.provenance,
            validation,
        ));
    }
    Ok(validation)
}

fn failure(
    reason: SharedTrimEdgeFailureReason,
    provenance: &'static str,
    validation: SharedTrimEdgeValidation,
) -> SharedTrimEdgeFailure {
    SharedTrimEdgeFailure {
        reason,
        provenance,
        validation,
    }
}

fn missing(provenance: &'static str) -> SharedTrimEdgeFailure {
    failure(
        SharedTrimEdgeFailureReason::MissingTopology,
        provenance,
        SharedTrimEdgeValidation::default(),
    )
}

fn surface_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<&SurfaceGeometry> {
    let loop_ = brep.loops.get(brep.half_edges.get(half_edge)?.loop_)?;
    brep.faces
        .get(loop_.face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)
}

fn canonical_edge_vertices(
    brep: &BRep,
    canonical_half_edge: HalfEdgeHandle,
) -> Option<(Point3, Point3)> {
    let he = brep.half_edges.get(canonical_half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some((
        brep.vertices.get(he.origin)?.position,
        brep.vertices.get(twin.origin)?.position,
    ))
}

fn edge_endpoint_gap(curve: &CurveOnEdge, start: Point3, end: Point3) -> f64 {
    curve
        .geometry
        .point_at(curve.t_start)
        .distance_to(&start)
        .max(curve.geometry.point_at(curve.t_end).distance_to(&end))
}

fn curve_surface_deviation(
    curve: &CurveOnEdge,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    forward: bool,
) -> f64 {
    let (pc_start, pc_end) = pcurve.parameter_range();
    let mut worst = 0.0_f64;
    for i in 0..9 {
        let frac = i as f64 / 8.0;
        let curve_frac = if forward { frac } else { 1.0 - frac };
        let curve_t = curve.t_start + curve_frac * (curve.t_end - curve.t_start);
        let pc_t = pc_start + frac * (pc_end - pc_start);
        let uv = pcurve.point_at(pc_t);
        let p_surface = surface.point_at(uv.x, uv.y);
        worst = worst.max(curve.geometry.point_at(curve_t).distance_to(&p_surface));
    }
    worst
}

fn side_agreement_deviation(
    pcurve_a: &Pcurve,
    surface_a: &SurfaceGeometry,
    pcurve_b: &Pcurve,
    surface_b: &SurfaceGeometry,
) -> f64 {
    let (a_start, a_end) = pcurve_a.parameter_range();
    let (b_start, b_end) = pcurve_b.parameter_range();
    let mut worst = 0.0_f64;
    for i in 0..9 {
        let frac = i as f64 / 8.0;
        let a_t = a_start + frac * (a_end - a_start);
        let b_t = b_start + (1.0 - frac) * (b_end - b_start);
        let a_uv = pcurve_a.point_at(a_t);
        let b_uv = pcurve_b.point_at(b_t);
        let a_point = surface_a.point_at(a_uv.x, a_uv.y);
        let b_point = surface_b.point_at(b_uv.x, b_uv.y);
        worst = worst.max(a_point.distance_to(&b_point));
    }
    worst
}

#[cfg(test)]
#[path = "shared_trim_edge_tests.rs"]
mod tests;
