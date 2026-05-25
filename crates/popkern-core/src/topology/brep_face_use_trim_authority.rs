//! Face-use trim authority diagnostics for shared-edge topology.

use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle};

/// Options for face-use trim authority diagnostics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepFaceUseTrimAuthorityOptions {
    /// Maximum accepted 3D gap for pcurve endpoint and sampled trace checks.
    pub tolerance: f64,
    /// Number of samples used for curve/surface trace agreement.
    pub samples: usize,
}

/// Aggregate report for face-use trim authority diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFaceUseTrimAuthorityReport {
    /// Number of face loop half-edge uses inspected.
    pub inspected_face_uses: usize,
    /// Number of face uses whose loop ownership does not match the half-edge authority.
    pub loop_ownership_mismatches: usize,
    /// Number of face uses with no pcurve authority.
    pub missing_pcurves: usize,
    /// Number of face uses whose pcurve is not valid on the face-use surface.
    pub pcurve_surface_mismatches: usize,
    /// Number of face uses resolved by explicit face-use-local pcurve authority.
    pub explicit_face_use_authorities: usize,
    /// Maximum measured endpoint gap.
    pub max_endpoint_gap: f64,
    /// Maximum measured sampled trace gap.
    pub max_sampled_gap: f64,
    /// Per-face-use findings.
    pub findings: Vec<BrepFaceUseTrimAuthorityFinding>,
}

impl BrepFaceUseTrimAuthorityReport {
    /// True when every inspected face use has coherent trim authority.
    pub fn is_clean(&self) -> bool {
        self.findings.is_empty()
    }
}

/// One face-use trim authority finding.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepFaceUseTrimAuthorityFinding {
    /// Face that references the loop use.
    pub face: FaceHandle,
    /// Loop referenced by the face.
    pub loop_: LoopHandle,
    /// Half-edge used by the loop.
    pub half_edge: HalfEdgeHandle,
    /// Edge referenced by the half-edge.
    pub edge: EdgeHandle,
    /// Finding kind.
    pub kind: BrepFaceUseTrimAuthorityFindingKind,
    /// Endpoint gap measured on the face-use surface.
    pub endpoint_gap: f64,
    /// Sampled curve/surface trace gap measured on the face-use surface.
    pub sampled_gap: f64,
}

/// Machine-readable face-use trim authority finding kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepFaceUseTrimAuthorityFindingKind {
    /// Face references a loop or half-edge whose stored ownership points elsewhere.
    LoopOwnershipMismatch,
    /// The face-use half-edge has no pcurve.
    MissingPcurve,
    /// The pcurve does not reproduce the edge on the face-use surface.
    PcurveSurfaceMismatch,
}

/// Analyze whether each face loop use has trim authority valid for that face.
pub fn analyze_brep_face_use_trim_authority(
    brep: &BRep,
    options: BrepFaceUseTrimAuthorityOptions,
) -> BrepFaceUseTrimAuthorityReport {
    let mut report = BrepFaceUseTrimAuthorityReport {
        inspected_face_uses: 0,
        loop_ownership_mismatches: 0,
        missing_pcurves: 0,
        pcurve_surface_mismatches: 0,
        explicit_face_use_authorities: 0,
        max_endpoint_gap: 0.0,
        max_sampled_gap: 0.0,
        findings: Vec::new(),
    };
    for (face_handle, face) in &brep.faces {
        inspect_loop(brep, face_handle, face.outer_loop, options, &mut report);
        for &inner_loop in &face.inner_loops {
            inspect_loop(brep, face_handle, inner_loop, options, &mut report);
        }
    }
    report
}

fn inspect_loop(
    brep: &BRep,
    face: FaceHandle,
    loop_: LoopHandle,
    options: BrepFaceUseTrimAuthorityOptions,
    report: &mut BrepFaceUseTrimAuthorityReport,
) {
    let Some(loop_data) = brep.loops.get(loop_) else {
        return;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return;
    }
    let mut cursor = start;
    let bound = brep.half_edges.len().max(1);
    for _ in 0..=bound {
        let Some(half_edge) = brep.half_edges.get(cursor) else {
            break;
        };
        report.inspected_face_uses += 1;
        let explicit_authority = brep.face_use_pcurve(face, cursor).is_some();
        if explicit_authority {
            report.explicit_face_use_authorities += 1;
        }
        if !explicit_authority && (half_edge.loop_ != loop_ || loop_data.face != face) {
            push_finding(
                report,
                face,
                loop_,
                cursor,
                half_edge.edge,
                BrepFaceUseTrimAuthorityFindingKind::LoopOwnershipMismatch,
                (0.0, 0.0),
            );
        }
        inspect_half_edge_face_use(brep, face, loop_, cursor, options, report);
        cursor = half_edge.next;
        if cursor == start {
            break;
        }
    }
}

fn inspect_half_edge_face_use(
    brep: &BRep,
    face: FaceHandle,
    loop_: LoopHandle,
    half_edge: HalfEdgeHandle,
    options: BrepFaceUseTrimAuthorityOptions,
    report: &mut BrepFaceUseTrimAuthorityReport,
) {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return;
    };
    let Some(pcurve) = brep
        .face_use_pcurve(face, half_edge)
        .or(half_edge_data.pcurve.as_ref())
    else {
        push_finding(
            report,
            face,
            loop_,
            half_edge,
            half_edge_data.edge,
            BrepFaceUseTrimAuthorityFindingKind::MissingPcurve,
            (f64::INFINITY, f64::INFINITY),
        );
        return;
    };
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return;
    };
    let Some(curve) = brep
        .edges
        .get(half_edge_data.edge)
        .and_then(|edge| edge.curve.as_ref())
    else {
        return;
    };
    let endpoint_gap =
        pcurve_endpoint_gap(brep, half_edge, pcurve, surface).unwrap_or(f64::INFINITY);
    let sampled_gap = pcurve_sampled_gap(brep, half_edge, pcurve, surface, curve, options.samples)
        .unwrap_or(f64::INFINITY);
    report.max_endpoint_gap = report.max_endpoint_gap.max(endpoint_gap);
    report.max_sampled_gap = report.max_sampled_gap.max(sampled_gap);
    if endpoint_gap.max(sampled_gap) > options.tolerance {
        push_finding(
            report,
            face,
            loop_,
            half_edge,
            half_edge_data.edge,
            BrepFaceUseTrimAuthorityFindingKind::PcurveSurfaceMismatch,
            (endpoint_gap, sampled_gap),
        );
    }
}

fn pcurve_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
) -> Option<f64> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let next = brep.half_edges.get(half_edge_data.next)?;
    let start = brep.vertices.get(half_edge_data.origin)?.position;
    let end = brep.vertices.get(next.origin)?.position;
    let start_uv = pcurve.point_at(pcurve.t_start);
    let end_uv = pcurve.point_at(pcurve.t_end);
    Some(
        surface
            .point_at(start_uv.x, start_uv.y)
            .distance_to(&start)
            .max(surface.point_at(end_uv.x, end_uv.y).distance_to(&end)),
    )
}

fn pcurve_sampled_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    curve: &CurveOnEdge,
    samples: usize,
) -> Option<f64> {
    let forward = half_edge_curve_forward(brep, half_edge, curve)?;
    let samples = samples.max(2);
    let mut max_gap = 0.0_f64;
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let pcurve_t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * fraction;
        let uv = pcurve.point_at(pcurve_t);
        max_gap = max_gap.max(
            curve
                .geometry
                .point_at(curve_t)
                .distance_to(&surface.point_at(uv.x, uv.y)),
        );
    }
    Some(max_gap)
}

fn half_edge_curve_forward(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
) -> Option<bool> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    Some(
        origin.distance_to(&curve.geometry.point_at(curve.t_start))
            <= origin.distance_to(&curve.geometry.point_at(curve.t_end)),
    )
}

fn push_finding(
    report: &mut BrepFaceUseTrimAuthorityReport,
    face: FaceHandle,
    loop_: LoopHandle,
    half_edge: HalfEdgeHandle,
    edge: EdgeHandle,
    kind: BrepFaceUseTrimAuthorityFindingKind,
    gaps: (f64, f64),
) {
    match kind {
        BrepFaceUseTrimAuthorityFindingKind::LoopOwnershipMismatch => {
            report.loop_ownership_mismatches += 1
        }
        BrepFaceUseTrimAuthorityFindingKind::MissingPcurve => report.missing_pcurves += 1,
        BrepFaceUseTrimAuthorityFindingKind::PcurveSurfaceMismatch => {
            report.pcurve_surface_mismatches += 1
        }
    }
    report.findings.push(BrepFaceUseTrimAuthorityFinding {
        face,
        loop_,
        half_edge,
        edge,
        kind,
        endpoint_gap: gaps.0,
        sampled_gap: gaps.1,
    });
}

#[cfg(test)]
#[path = "brep_face_use_trim_authority_tests.rs"]
mod tests;
