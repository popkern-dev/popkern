//! Blend-boundary seam re-anchoring for N-quad cap rings.
//!
//! The N-gon corner face inherits its outer seam vertices and edge
//! curves from the mitred corner construction, which places them up
//! to ~4e-2 (hexagon6_variable, PK-087 seam census) away from the
//! `BlendBoundary` segments that every promoted cap family
//! interpolates. Downstream boundary-owned seam authoring then
//! declines on `endpoint_gap` and cap/support acceptance measures
//! the mitred-vs-ideal offset instead of the cap fit quality.
//!
//! This pass re-anchors the live ring to the blend boundary:
//! 1. match each cap outer seam edge to its boundary segment
//!    (complete bijection required — partial re-anchoring would mix
//!    mitred and ideal frames);
//! 2. snap each ring vertex to the segment junction (incident NURBS
//!    curves get an endpoint-CP drag, incident analytic curves a
//!    closest-point re-trim);
//! 3. write the exact segment curve into each seam edge and refit
//!    the support-side pcurve.
//!
//! Transactional: the pass measures the worst vertex-to-junction gap
//! first and no-ops when the ring is already anchored (within the
//! boundary-authoring tolerance); any postcondition failure restores
//! the entry snapshot. All declines are recorded as structured
//! `ApexConstructionDiagnostic` facts by [`reanchor_diagnostic`].
//!
//! Tolerances: segment identity is structural (global cyclic ring
//! alignment — see the match module), guarded only by a
//! quarter-perimeter sanity bound; junction coherence must reach
//! twice the boundary-authoring tolerance (vertices land on the
//! junction midpoint, so each curve endpoint is half a spread from
//! its vertex), and the final endpoint postcondition re-measures
//! every rewired edge at the boundary-authoring tolerance.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::fitting;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

#[path = "three_quad_blend_seam_reanchor_edit.rs"]
mod three_quad_blend_seam_reanchor_edit;
#[path = "three_quad_blend_seam_reanchor_match.rs"]
mod three_quad_blend_seam_reanchor_match;
use three_quad_blend_seam_reanchor_edit::{
    edge_curve_endpoint_gap, rewire_seam_edge, snap_vertex_with_retrim,
};
use three_quad_blend_seam_reanchor_match::match_ring_to_segments;

/// Measured outcome of one re-anchoring pass.
pub(crate) struct BlendSeamReanchorOutcome {
    pub matched_jobs: usize,
    pub expected_jobs: usize,
    pub vertices_snapped: usize,
    pub edges_rewired: usize,
    pub worst_gap_before: f64,
    pub worst_gap_after: f64,
    pub worst_junction_spread: f64,
    pub decline: Option<&'static str>,
}

struct ReanchorJob {
    half_edge: HalfEdgeHandle,
    edge: EdgeHandle,
    segment_index: usize,
    reversed: bool,
    v_origin: VertexHandle,
    v_dest: VertexHandle,
    target_origin: Point3,
    target_dest: Point3,
}

fn authoring_tolerance(brep: &BRep) -> f64 {
    brep.tolerance.linear().max(2.0e-4)
}

fn segment_curve(segment: &BoundarySegment) -> &NurbsCurve {
    match segment {
        BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => curve,
    }
}

/// Re-anchor the outer seam ring of an N-quad cap decomposition to
/// the blend-boundary segments. See the module docs for the model;
/// returns a measured outcome (never panics, never partially
/// applies — on any decline the BRep is exactly the entry state).
pub(crate) fn reanchor_blend_seams_to_boundary(
    brep: &mut BRep,
    boundary: &BlendBoundary,
    sub_faces: &[FaceHandle],
) -> BlendSeamReanchorOutcome {
    let outcome = reanchor_inner(brep, boundary, sub_faces);
    trace_outcome(&outcome);
    outcome
}

fn reanchor_inner(
    brep: &mut BRep,
    boundary: &BlendBoundary,
    sub_faces: &[FaceHandle],
) -> BlendSeamReanchorOutcome {
    let expected_jobs = boundary.segments.len();
    let mut outcome = BlendSeamReanchorOutcome {
        matched_jobs: 0,
        expected_jobs,
        vertices_snapped: 0,
        edges_rewired: 0,
        worst_gap_before: 0.0,
        worst_gap_after: 0.0,
        worst_junction_spread: 0.0,
        decline: None,
    };
    if expected_jobs != 2 * sub_faces.len() || expected_jobs == 0 {
        outcome.decline = Some("segment_count_mismatch");
        return outcome;
    }
    let matched = match match_ring_to_segments(brep, boundary, sub_faces) {
        Ok(matched) => matched,
        Err(reason) => {
            outcome.decline = Some(reason);
            return outcome;
        }
    };
    let jobs = matched.jobs;
    outcome.matched_jobs = jobs.len();
    outcome.worst_junction_spread = matched.worst_junction_spread;
    outcome.worst_gap_before = jobs
        .iter()
        .flat_map(|job| {
            [
                brep.vertices[job.v_origin]
                    .position
                    .distance_to(&job.target_origin),
                brep.vertices[job.v_dest]
                    .position
                    .distance_to(&job.target_dest),
            ]
        })
        .fold(0.0_f64, f64::max);
    outcome.worst_gap_after = outcome.worst_gap_before;
    let tolerance = authoring_tolerance(brep);
    if outcome.worst_gap_before <= tolerance {
        // Ring already anchored — nothing to do (the common case for
        // fixtures whose corner topology was authored on the boundary).
        return outcome;
    }

    // Junction pre-gate: vertices land on the junction midpoint, so
    // each rewired curve endpoint sits half a spread from its vertex
    // — past twice the authoring tolerance the postcondition cannot
    // hold and the decline reason should say why.
    if outcome.worst_junction_spread > 2.0 * tolerance {
        outcome.decline = Some("junction_spread_exceeds_gate");
        return outcome;
    }

    // Both jobs at a shared vertex carry the identical junction
    // midpoint (assigned in the match pass), so a plain first-claim
    // dedupe is exact.
    let mut vertex_targets: Vec<(VertexHandle, Point3)> = Vec::new();
    for job in &jobs {
        for (v, target) in [
            (job.v_origin, job.target_origin),
            (job.v_dest, job.target_dest),
        ] {
            if !vertex_targets.iter().any(|(vh, _)| *vh == v) {
                vertex_targets.push((v, target));
            }
        }
    }

    let snapshot = brep.clone();
    let seam_edges: Vec<EdgeHandle> = jobs.iter().map(|job| job.edge).collect();
    for (v, target) in &vertex_targets {
        if let Err(reason) = snap_vertex_with_retrim(brep, *v, *target, tolerance, &seam_edges) {
            *brep = snapshot;
            outcome.decline = Some(reason);
            return outcome;
        }
        outcome.vertices_snapped += 1;
    }
    for job in &jobs {
        if let Err(reason) = rewire_seam_edge(brep, boundary, job) {
            *brep = snapshot;
            outcome.decline = Some(reason);
            outcome.vertices_snapped = 0;
            outcome.edges_rewired = 0;
            return outcome;
        }
        outcome.edges_rewired += 1;
    }

    // Postcondition: every rewired edge's curve endpoints coincide
    // with its (snapped) vertices at trim precision.
    let mut worst_after = 0.0_f64;
    for job in &jobs {
        let Some(gap) = edge_curve_endpoint_gap(brep, job) else {
            *brep = snapshot;
            outcome.decline = Some("postcondition_edge_unreadable");
            outcome.vertices_snapped = 0;
            outcome.edges_rewired = 0;
            return outcome;
        };
        worst_after = worst_after.max(gap);
    }
    if worst_after > tolerance {
        *brep = snapshot;
        outcome.decline = Some("postcondition_endpoint_gap");
        outcome.vertices_snapped = 0;
        outcome.edges_rewired = 0;
        return outcome;
    }
    outcome.worst_gap_after = worst_after;
    outcome
}

fn trace_outcome(outcome: &BlendSeamReanchorOutcome) {
    if crate::direct_apex_trace::cap_g0_enabled() {
        eprintln!(
            "[direct-apex-trace cap-g0] blend-seam reanchor: decline={:?} jobs={}/{} snapped={} \
             rewired={} gap {:.3e} -> {:.3e} junction_spread={:.3e}",
            outcome.decline,
            outcome.matched_jobs,
            outcome.expected_jobs,
            outcome.vertices_snapped,
            outcome.edges_rewired,
            outcome.worst_gap_before,
            outcome.worst_gap_after,
            outcome.worst_junction_spread,
        );
    }
}

/// Structured fact for the diagnostics pipeline — one per apply,
/// covering applied, no-op, and declined outcomes.
pub(crate) fn reanchor_diagnostic(
    outcome: &BlendSeamReanchorOutcome,
    n_sides: usize,
) -> ApexConstructionDiagnostic {
    let stage_outcome = match outcome.decline {
        Some(reason) => reason,
        None if outcome.edges_rewired > 0 => "applied",
        None => "already_anchored",
    };
    ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
        stage: format!("blend_seam_reanchor.{stage_outcome}"),
        code: "blend_seam_reanchor.fact".to_string(),
        detail_code: Some(format!(
            "matched={};expected={};snapped={};rewired={};gap_before={:.3e};gap_after={:.3e};junction_spread={:.3e}",
            outcome.matched_jobs,
            outcome.expected_jobs,
            outcome.vertices_snapped,
            outcome.edges_rewired,
            outcome.worst_gap_before,
            outcome.worst_gap_after,
            outcome.worst_junction_spread,
        )),
        sector_name: format!("blend_seam_reanchor_ring_n{n_sides}"),
        edge: None,
        support_half_edge: None,
        support_face: None,
        residual: Some(outcome.worst_gap_before),
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
        shared_fact_ids: Vec::new(),
    })
}

fn average_point(points: &[Point3]) -> Point3 {
    let n = points.len().max(1) as f64;
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for p in points {
        x += p.x;
        y += p.y;
        z += p.z;
    }
    Point3::new(x / n, y / n, z / n)
}
