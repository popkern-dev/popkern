//! Shared-trim commit for exact-boundary cap/cap spoke edges.

use std::collections::BTreeSet;

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, SharedTrimEdgeInput,
    commit_shared_trim_edge_transactionally,
};

use super::commit::{ExactBoundaryCommitDeclineReason, ExactBoundaryCommitDiagnostic};
use super::sector::ExactBoundaryCapSectorInput;

const PCURVE_SAMPLES: usize = 9;
const COMMIT_TOLERANCE: f64 = 1.0e-6;

#[allow(clippy::result_large_err)]
pub(super) fn commit_cap_spoke_shared_trims(
    brep: &mut BRep,
    sectors: &[ExactBoundaryCapSectorInput],
) -> Result<usize, ExactBoundaryCommitDiagnostic> {
    let sector_faces = sectors
        .iter()
        .map(|sector| sector.cap_face)
        .collect::<BTreeSet<_>>();
    let outer_edges = sectors
        .iter()
        .flat_map(|sector| std::iter::once(sector.fillet_end_edge).chain(sector.face_edge))
        .collect::<BTreeSet<_>>();
    let mut committed = BTreeSet::new();
    for sector in sectors {
        for half_edge in face_outer_loop_half_edges(brep, sector.cap_face) {
            let Some(edge) = brep
                .half_edges
                .get(half_edge)
                .map(|half_edge| half_edge.edge)
            else {
                continue;
            };
            if outer_edges.contains(&edge) || committed.contains(&edge) {
                continue;
            }
            let Some(twin) = brep
                .half_edges
                .get(half_edge)
                .map(|half_edge| half_edge.twin)
            else {
                continue;
            };
            let Some(twin_face) = face_for_half_edge(brep, twin) else {
                continue;
            };
            if !sector_faces.contains(&twin_face) {
                continue;
            }
            commit_one_cap_spoke(brep, sector, edge, half_edge, twin)?;
            committed.insert(edge);
        }
    }
    Ok(committed.len())
}

#[allow(clippy::result_large_err)]
fn commit_one_cap_spoke(
    brep: &mut BRep,
    sector: &ExactBoundaryCapSectorInput,
    edge: EdgeHandle,
    half_edge_a: HalfEdgeHandle,
    half_edge_b: HalfEdgeHandle,
) -> Result<(), ExactBoundaryCommitDiagnostic> {
    let curve = brep
        .edges
        .get(edge)
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| diagnostic(sector, edge, half_edge_a, half_edge_b, None))?;
    let surface_a = surface_for_half_edge(brep, half_edge_a)
        .ok_or_else(|| diagnostic(sector, edge, half_edge_a, half_edge_b, None))?;
    let surface_b = surface_for_half_edge(brep, half_edge_b)
        .ok_or_else(|| diagnostic(sector, edge, half_edge_a, half_edge_b, None))?;
    let pcurve_a = project_curve_to_half_edge_pcurve(brep, edge, half_edge_a, &curve, surface_a)
        .ok_or_else(|| diagnostic(sector, edge, half_edge_a, half_edge_b, None))?;
    let pcurve_b = project_curve_to_half_edge_pcurve(brep, edge, half_edge_b, &curve, surface_b)
        .ok_or_else(|| diagnostic(sector, edge, half_edge_a, half_edge_b, None))?;
    if let Err(failure) = commit_shared_trim_edge_transactionally(
        brep,
        SharedTrimEdgeInput {
            edge,
            half_edge_a,
            half_edge_b,
            curve,
            pcurve_a,
            pcurve_b,
            tolerance: COMMIT_TOLERANCE,
            validate_brep: false,
            provenance: "exact_boundary_cap_spoke_shared_trim",
        },
    ) {
        return Err(diagnostic(
            sector,
            edge,
            half_edge_a,
            half_edge_b,
            Some(
                failure
                    .validation
                    .max_curve_surface_deviation
                    .max(failure.validation.max_side_agreement_deviation)
                    .max(failure.validation.max_endpoint_gap),
            ),
        ));
    }
    Ok(())
}

fn project_curve_to_half_edge_pcurve(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    surface: &SurfaceGeometry,
) -> Option<Pcurve> {
    let forward = half_edge == brep.edges.get(edge)?.half_edge;
    let mut uv_points = Vec::with_capacity(PCURVE_SAMPLES);
    let mut parameters = Vec::with_capacity(PCURVE_SAMPLES);
    for index in 0..PCURVE_SAMPLES {
        let fraction = index as f64 / (PCURVE_SAMPLES - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + curve_fraction * (curve.t_end - curve.t_start);
        let point = curve.geometry.point_at(curve_t);
        let (u, v) = surface.closest_point(&point)?;
        uv_points.push(Point3::new(u, v, 0.0));
        parameters.push(fraction);
    }
    let degree = 3.min(uv_points.len().saturating_sub(1));
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &parameters, degree)
        .ok()
        .or_else(|| linear_uv_curve(&uv_points))
        .map(Pcurve::from_nurbs)
}

fn linear_uv_curve(points: &[Point3]) -> Option<NurbsCurve> {
    Some(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![*points.first()?, *points.last()?],
        weights: vec![1.0, 1.0],
    })
}

fn face_outer_loop_half_edges(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let Some(loop_handle) = brep.faces.get(face).map(|face| face.outer_loop) else {
        return Vec::new();
    };
    collect_loop_half_edges(brep, loop_handle)
}

fn collect_loop_half_edges(brep: &BRep, loop_handle: LoopHandle) -> Vec<HalfEdgeHandle> {
    let Some(start) = brep.loops.get(loop_handle).map(|loop_| loop_.half_edge) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        let Some(next) = brep.half_edges.get(cursor).map(|half_edge| half_edge.next) else {
            return Vec::new();
        };
        cursor = next;
        if cursor == start {
            return out;
        }
    }
    Vec::new()
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    brep.half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_data| loop_data.face)
}

fn surface_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<&SurfaceGeometry> {
    let face = face_for_half_edge(brep, half_edge)?;
    brep.faces
        .get(face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)
}

fn diagnostic(
    sector: &ExactBoundaryCapSectorInput,
    edge: EdgeHandle,
    half_edge_a: HalfEdgeHandle,
    half_edge_b: HalfEdgeHandle,
    residual: Option<f64>,
) -> ExactBoundaryCommitDiagnostic {
    ExactBoundaryCommitDiagnostic {
        reason: ExactBoundaryCommitDeclineReason::SharedTrimValidationFailed,
        code: ExactBoundaryCommitDeclineReason::SharedTrimValidationFailed.code(),
        sector_name: Some(sector.sector_name.clone()),
        face: Some(sector.cap_face),
        edge: Some(edge),
        cap_half_edge: Some(half_edge_a),
        support_half_edge: Some(half_edge_b),
        shared_trim_failure: None,
        residual,
        outer_side: None,
        sample_parameter: None,
        tolerance: None,
    }
}
