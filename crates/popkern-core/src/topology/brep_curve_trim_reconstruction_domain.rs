//! Bounded-surface domain repair for curve-backed trim reconstruction.

use crate::geometry::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::topology::brep_curve_trim_reconstruction::{
    BrepCurveTrimReconstructionOptions, ExistingTrimStatus, existing_trim_status,
};
use crate::topology::brep_curve_trim_reconstruction_helpers::project_edge_curve_to_surface_pcurve;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};

pub(super) struct DomainExpansionCommit {
    pub(super) face: FaceHandle,
    pub(super) previous_surface: Option<SurfaceOnFace>,
}

pub(super) fn expand_repaired_face_domain(
    original: &BRep,
    repaired: &mut BRep,
    half_edge: HalfEdgeHandle,
    options: BrepCurveTrimReconstructionOptions,
) -> Option<DomainExpansionCommit> {
    let candidate = project_edge_curve_to_surface_pcurve(original, half_edge, options.samples)?;
    let face = original
        .half_edges
        .get(half_edge)
        .and_then(|he| original.loops.get(he.loop_))
        .map(|loop_| loop_.face)?;
    let (u_range, v_range) = expanded_face_uv_ranges(original, face, &candidate)?;
    let source_surface = original.faces[face].surface.as_ref()?.geometry.clone();
    let expanded = source_surface.approximate_to_nurbs_in_range(u_range, v_range);
    let previous_surface = {
        let face_data = repaired.faces.get_mut(face)?;
        let previous_surface = face_data.surface.clone();
        face_data.surface.as_mut()?.geometry = SurfaceGeometry::Nurbs(expanded);
        previous_surface
    };
    if validate_preserved_face_trims(original, repaired, face, half_edge, options) {
        Some(DomainExpansionCommit {
            face,
            previous_surface,
        })
    } else {
        repaired.faces.get_mut(face)?.surface = previous_surface;
        None
    }
}

fn expanded_face_uv_ranges(
    brep: &BRep,
    face: FaceHandle,
    candidate: &crate::geometry::pcurve::Pcurve,
) -> Option<((f64, f64), (f64, f64))> {
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let ((mut u_min, mut u_max), (mut v_min, mut v_max)) = surface.parameter_range();
    accumulate_pcurve_range(candidate, &mut u_min, &mut u_max, &mut v_min, &mut v_max);
    for half_edge in face_half_edges(brep, face)? {
        if let Some(pcurve) = brep.half_edges[half_edge].pcurve.as_ref() {
            accumulate_pcurve_range(pcurve, &mut u_min, &mut u_max, &mut v_min, &mut v_max);
        }
    }
    let u_pad = ((u_max - u_min).abs() * 1.0e-9).max(1.0e-9);
    let v_pad = ((v_max - v_min).abs() * 1.0e-9).max(1.0e-9);
    Some((
        (u_min - u_pad, u_max + u_pad),
        (v_min - v_pad, v_max + v_pad),
    ))
}

fn accumulate_pcurve_range(
    pcurve: &crate::geometry::pcurve::Pcurve,
    u_min: &mut f64,
    u_max: &mut f64,
    v_min: &mut f64,
    v_max: &mut f64,
) {
    let samples = 9usize;
    for index in 0..samples {
        let s = index as f64 / (samples - 1) as f64;
        let uv = pcurve.point_at(pcurve.t_start + (pcurve.t_end - pcurve.t_start) * s);
        *u_min = u_min.min(uv.x);
        *u_max = u_max.max(uv.x);
        *v_min = v_min.min(uv.y);
        *v_max = v_max.max(uv.y);
    }
}

fn validate_preserved_face_trims(
    original: &BRep,
    repaired: &BRep,
    face: FaceHandle,
    excluded_half_edge: HalfEdgeHandle,
    options: BrepCurveTrimReconstructionOptions,
) -> bool {
    let Some(half_edges) = face_half_edges(repaired, face) else {
        return false;
    };
    half_edges
        .into_iter()
        .filter(|half_edge| {
            *half_edge != excluded_half_edge
                && !is_stale_reconstruction_target(original, *half_edge, options)
        })
        .all(|half_edge| {
            !matches!(
                existing_trim_status(repaired, half_edge, options),
                ExistingTrimStatus::Stale
            )
        })
}

fn is_stale_reconstruction_target(
    original: &BRep,
    half_edge: HalfEdgeHandle,
    options: BrepCurveTrimReconstructionOptions,
) -> bool {
    options.replace_stale_existing
        && matches!(
            existing_trim_status(original, half_edge, options),
            ExistingTrimStatus::Stale
        )
}

fn face_half_edges(brep: &BRep, face: FaceHandle) -> Option<Vec<HalfEdgeHandle>> {
    let start = brep.loops.get(brep.faces.get(face)?.outer_loop)?.half_edge;
    let mut cursor = start;
    let mut out = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}
