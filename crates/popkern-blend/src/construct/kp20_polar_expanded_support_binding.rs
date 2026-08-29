//! Surface binding contracts for expanded KP20 support topology.
#![allow(dead_code)]

use popkern_core::error::CoreError;
use popkern_core::geometry::edge_triple::BoundarySide;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::BRep;
use popkern_core::topology::FaceHandle;

#[path = "kp20_polar_expanded_support_unit_domain.rs"]
mod expanded_support_unit_domain;
use expanded_support_unit_domain::unit_domain_loop_candidates;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Kp20ExpandedSupportSurfaceSlot {
    pub(crate) face: FaceHandle,
    pub(crate) boundary_edge_count: usize,
    pub(crate) assignment: Kp20ExpandedSupportSurfaceAssignment,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Kp20ExpandedSupportSurfaceAssignment {
    FullDomainRectangularPatch,
    TrimmedMacroFace,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Kp20ExpandedSupportSurfaceBindingPlan {
    pub(crate) bindings: Vec<Kp20ExpandedSupportSurfaceBinding>,
    pub(crate) full_domain_binding_count: usize,
    pub(crate) trimmed_macro_binding_count: usize,
    pub(crate) reused_source_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Kp20ExpandedSupportSurfaceBinding {
    pub(crate) face: FaceHandle,
    pub(crate) support_surface_index: usize,
    pub(crate) assignment: Kp20ExpandedSupportSurfaceAssignment,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Kp20ExpandedSupportSharedBindingSummary {
    pub(crate) topology_slot_count: usize,
    pub(crate) support_surface_count: usize,
    pub(crate) trimmed_slot_count: usize,
    pub(crate) reused_source_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Kp20ExpandedSupportCommitShape {
    pub(crate) face_count: usize,
    pub(crate) pentagon_count: usize,
    pub(crate) alignment: &'static str,
}

pub(crate) fn expanded_support_commit_shape(
    sector_count: usize,
    target_count: usize,
    rectangular_surface_count: usize,
) -> Kp20ExpandedSupportCommitShape {
    if target_count != sector_count || !sector_count.is_multiple_of(2) {
        return Kp20ExpandedSupportCommitShape {
            face_count: 0,
            pentagon_count: 0,
            alignment: "expanded_support_topology_shape_unmeasured",
        };
    }
    let face_count = sector_count * 15 / 2;
    let pentagon_count = sector_count;
    let alignment = if face_count == rectangular_surface_count && pentagon_count == 0 {
        "expanded_support_topology_rectangular_surfaces_aligned"
    } else if summarize_kp20_expanded_support_shared_trim_bindings(
        face_count,
        pentagon_count,
        rectangular_surface_count,
    )
    .is_ok()
    {
        "expanded_support_topology_shared_trim_binding_ready"
    } else {
        "expanded_support_topology_requires_trimmed_macro_face_assignment"
    };
    Kp20ExpandedSupportCommitShape {
        face_count,
        pentagon_count,
        alignment,
    }
}

pub(crate) fn summarize_kp20_expanded_support_shared_trim_bindings(
    topology_slot_count: usize,
    trimmed_slot_count: usize,
    support_surface_count: usize,
) -> Result<Kp20ExpandedSupportSharedBindingSummary, CoreError> {
    if support_surface_count > topology_slot_count {
        return Ok(Kp20ExpandedSupportSharedBindingSummary {
            topology_slot_count,
            support_surface_count,
            trimmed_slot_count,
            reused_source_count: 0,
        });
    }
    let deficit = topology_slot_count - support_surface_count;
    if deficit > trimmed_slot_count {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support shared binding needs {deficit} trimmed-source reuses but only {trimmed_slot_count} trimmed macro slots exist",
        )));
    }
    Ok(Kp20ExpandedSupportSharedBindingSummary {
        topology_slot_count,
        support_surface_count,
        trimmed_slot_count,
        reused_source_count: deficit,
    })
}

pub(crate) fn plan_kp20_expanded_support_surface_bindings(
    slots: &[Kp20ExpandedSupportSurfaceSlot],
    support_surface_count: usize,
) -> Result<Kp20ExpandedSupportSurfaceBindingPlan, CoreError> {
    if support_surface_count < slots.len() {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support surface binding requires {} support surfaces for {} topology slots, got {support_surface_count}",
            slots.len(),
            slots.len(),
        )));
    }
    let source_map = (0..slots.len()).collect::<Vec<_>>();
    plan_kp20_expanded_support_surface_bindings_with_sources(
        slots,
        support_surface_count,
        &source_map,
    )
}

pub(crate) fn plan_kp20_expanded_support_shared_trim_bindings(
    slots: &[Kp20ExpandedSupportSurfaceSlot],
    support_surface_count: usize,
) -> Result<Kp20ExpandedSupportSurfaceBindingPlan, CoreError> {
    if support_surface_count > slots.len() {
        return plan_kp20_expanded_support_surface_bindings(slots, support_surface_count);
    }
    let trimmed_slots = slots
        .iter()
        .filter(|slot| slot.assignment == Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace)
        .count();
    let summary = summarize_kp20_expanded_support_shared_trim_bindings(
        slots.len(),
        trimmed_slots,
        support_surface_count,
    )?;
    let mut source_map = Vec::with_capacity(slots.len());
    let mut next_source = 0;
    let mut remaining_reuses = summary.reused_source_count;
    for slot in slots {
        if slot.assignment == Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace
            && remaining_reuses > 0
            && next_source > 0
        {
            source_map.push(next_source - 1);
            remaining_reuses -= 1;
        } else {
            source_map.push(next_source);
            next_source += 1;
        }
    }
    if next_source != support_surface_count || remaining_reuses != 0 {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support shared binding could not consume source inventory: used={next_source}, available={support_surface_count}, remaining_reuses={remaining_reuses}",
        )));
    }
    plan_kp20_expanded_support_surface_bindings_with_sources(
        slots,
        support_surface_count,
        &source_map,
    )
}

pub(crate) fn plan_kp20_expanded_support_surface_bindings_with_sources(
    slots: &[Kp20ExpandedSupportSurfaceSlot],
    support_surface_count: usize,
    source_map: &[usize],
) -> Result<Kp20ExpandedSupportSurfaceBindingPlan, CoreError> {
    if source_map.len() != slots.len() {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 expanded support shared binding expected {} source refs, got {}",
            slots.len(),
            source_map.len(),
        )));
    }
    let mut used = vec![false; support_surface_count];
    let mut reused_source_count = 0;
    let mut full_domain_binding_count = 0;
    let mut trimmed_macro_binding_count = 0;
    let mut bindings = Vec::with_capacity(slots.len());
    for (slot, &support_surface_index) in slots.iter().zip(source_map.iter()) {
        if support_surface_index >= support_surface_count {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 expanded support binding source {support_surface_index} out of range for {support_surface_count} surfaces",
            )));
        }
        if used[support_surface_index] {
            if slot.assignment != Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace {
                return Err(CoreError::PreconditionFailed(
                    "kp20 expanded support binding may only reuse sources for trimmed macro faces"
                        .into(),
                ));
            }
            reused_source_count += 1;
        }
        used[support_surface_index] = true;
        match slot.assignment {
            Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch => {
                full_domain_binding_count += 1;
            }
            Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace => {
                trimmed_macro_binding_count += 1;
            }
        }
        bindings.push(Kp20ExpandedSupportSurfaceBinding {
            face: slot.face,
            support_surface_index,
            assignment: slot.assignment,
        });
    }
    Ok(Kp20ExpandedSupportSurfaceBindingPlan {
        bindings,
        full_domain_binding_count,
        trimmed_macro_binding_count,
        reused_source_count,
    })
}

pub(crate) fn attach_kp20_expanded_support_surfaces(
    brep: &mut BRep,
    binding_plan: &Kp20ExpandedSupportSurfaceBindingPlan,
    support_surfaces: &[GregoryPatch],
) -> Result<(), CoreError> {
    let before = brep.clone();
    match attach_kp20_expanded_support_surfaces_inner(brep, binding_plan, support_surfaces) {
        Ok(()) => Ok(()),
        Err(error) => {
            *brep = before;
            Err(error)
        }
    }
}

fn attach_kp20_expanded_support_surfaces_inner(
    brep: &mut BRep,
    binding_plan: &Kp20ExpandedSupportSurfaceBindingPlan,
    support_surfaces: &[GregoryPatch],
) -> Result<(), CoreError> {
    for binding in &binding_plan.bindings {
        let Some(surface) = support_surfaces.get(binding.support_surface_index) else {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 expanded support surface attachment source {} out of range for {} surfaces",
                binding.support_surface_index,
                support_surfaces.len(),
            )));
        };
        let Some(face) = brep.faces.get_mut(binding.face) else {
            return Err(CoreError::InvalidHandle(
                "kp20 expanded support surface attachment face invalid".into(),
            ));
        };
        face.surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Gregory(surface.clone()),
            orientation: true,
        });
    }
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 expanded support surface attachment produced invalid BRep: {error}"
        ))
    })
}

pub(crate) fn author_kp20_expanded_support_pcurves(
    brep: &mut BRep,
    binding_plan: &Kp20ExpandedSupportSurfaceBindingPlan,
) -> Result<(), CoreError> {
    let before = brep.clone();
    match author_kp20_expanded_support_pcurves_inner(brep, binding_plan) {
        Ok(()) => Ok(()),
        Err(error) => {
            *brep = before;
            Err(error)
        }
    }
}

fn author_kp20_expanded_support_pcurves_inner(
    brep: &mut BRep,
    binding_plan: &Kp20ExpandedSupportSurfaceBindingPlan,
) -> Result<(), CoreError> {
    for binding in &binding_plan.bindings {
        if author_unit_domain_boundary_pcurves(brep, binding.face)? {
            verify_outer_loop_pcurves(brep, binding.face)?;
            continue;
        }
        clear_outer_loop_pcurves(brep, binding.face)?;
        brep.recompute_face_pcurves(binding.face)?;
        verify_outer_loop_pcurves(brep, binding.face)?;
    }
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 expanded support pcurve authoring produced invalid BRep: {error}"
        ))
    })
}

fn author_unit_domain_boundary_pcurves(
    brep: &mut BRep,
    face: FaceHandle,
) -> Result<bool, CoreError> {
    let half_edges = face_outer_half_edges(brep, face)?;
    if half_edges.len() < 4 {
        return Ok(false);
    }
    let Some(SurfaceGeometry::Gregory(surface)) = brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface_on_face| &surface_on_face.geometry)
    else {
        return Ok(false);
    };
    let Some(vertex_uvs) = best_unit_domain_loop_uvs(brep, surface, &half_edges) else {
        return Ok(false);
    };
    for (index, &half_edge) in half_edges.iter().enumerate() {
        let start = vertex_uvs[index];
        let end = vertex_uvs[(index + 1) % vertex_uvs.len()];
        brep.half_edges[half_edge].pcurve = Some(Pcurve::linear(start, end, 0.0, 1.0));
        record_unit_domain_boundary_hint(brep, half_edge, start, end);
    }
    Ok(true)
}

fn record_unit_domain_boundary_hint(
    brep: &mut BRep,
    half_edge: popkern_core::topology::HalfEdgeHandle,
    start: (f64, f64),
    end: (f64, f64),
) {
    if let Some(boundary) = unit_domain_boundary_side(start, end) {
        brep.set_surface_boundary_hint(half_edge, boundary);
    } else {
        brep.clear_surface_boundary_hint(half_edge);
    }
}

fn unit_domain_boundary_side(start: (f64, f64), end: (f64, f64)) -> Option<BoundarySide> {
    let span_u = (start.0 - end.0).abs();
    let span_v = (start.1 - end.1).abs();
    if close_unit(start.0, 0.0) && close_unit(end.0, 0.0) && span_v > 1.0e-9 {
        Some(BoundarySide::UMin)
    } else if close_unit(start.0, 1.0) && close_unit(end.0, 1.0) && span_v > 1.0e-9 {
        Some(BoundarySide::UMax)
    } else if close_unit(start.1, 0.0) && close_unit(end.1, 0.0) && span_u > 1.0e-9 {
        Some(BoundarySide::VMin)
    } else if close_unit(start.1, 1.0) && close_unit(end.1, 1.0) && span_u > 1.0e-9 {
        Some(BoundarySide::VMax)
    } else {
        None
    }
}

fn close_unit(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}

fn best_unit_domain_loop_uvs(
    brep: &BRep,
    surface: &GregoryPatch,
    half_edges: &[popkern_core::topology::HalfEdgeHandle],
) -> Option<Vec<(f64, f64)>> {
    let candidates = unit_domain_loop_candidates(half_edges.len());
    let mut best = None;
    for candidate in candidates {
        let mut residual = 0.0;
        for (&half_edge, &(u, v)) in half_edges.iter().zip(candidate.iter()) {
            let position = brep.vertices[brep.half_edges[half_edge].origin].position;
            residual += position.distance_to(&surface.point_at(u, v));
        }
        if best
            .as_ref()
            .is_none_or(|(_, best_residual)| residual < *best_residual)
        {
            best = Some((candidate, residual));
        }
    }
    best.map(|(candidate, _)| candidate)
}

fn face_outer_half_edges(
    brep: &BRep,
    face: FaceHandle,
) -> Result<Vec<popkern_core::topology::HalfEdgeHandle>, CoreError> {
    let face_data = brep.faces.get(face).ok_or_else(|| {
        CoreError::InvalidHandle("kp20 expanded support pcurve face invalid".into())
    })?;
    let loop_data = brep.loops.get(face_data.outer_loop).ok_or_else(|| {
        CoreError::InvalidHandle("kp20 expanded support pcurve loop invalid".into())
    })?;
    let start = loop_data.half_edge;
    let mut half_edges = Vec::new();
    let mut cursor = start;
    for _ in 0..brep.half_edges.len().max(1) {
        let half_edge = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::InvalidHandle("kp20 expanded support pcurve half-edge invalid".into())
        })?;
        half_edges.push(cursor);
        cursor = half_edge.next;
        if cursor == start {
            return Ok(half_edges);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 expanded support pcurve loop did not close".into(),
    ))
}

fn clear_outer_loop_pcurves(brep: &mut BRep, face: FaceHandle) -> Result<(), CoreError> {
    let face_data = brep.faces.get(face).ok_or_else(|| {
        CoreError::InvalidHandle("kp20 expanded support pcurve face invalid".into())
    })?;
    let loop_data = brep.loops.get(face_data.outer_loop).ok_or_else(|| {
        CoreError::InvalidHandle("kp20 expanded support pcurve loop invalid".into())
    })?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..brep.half_edges.len().max(1) {
        let next = {
            let half_edge = brep.half_edges.get_mut(cursor).ok_or_else(|| {
                CoreError::InvalidHandle("kp20 expanded support pcurve half-edge invalid".into())
            })?;
            half_edge.pcurve = None;
            half_edge.next
        };
        cursor = next;
        if cursor == start {
            return Ok(());
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 expanded support pcurve loop did not close".into(),
    ))
}

fn verify_outer_loop_pcurves(brep: &BRep, face: FaceHandle) -> Result<(), CoreError> {
    let face_data = brep.faces.get(face).ok_or_else(|| {
        CoreError::InvalidHandle("kp20 expanded support pcurve face invalid".into())
    })?;
    let loop_data = brep.loops.get(face_data.outer_loop).ok_or_else(|| {
        CoreError::InvalidHandle("kp20 expanded support pcurve loop invalid".into())
    })?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..brep.half_edges.len().max(1) {
        let half_edge = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::InvalidHandle("kp20 expanded support pcurve half-edge invalid".into())
        })?;
        if half_edge.pcurve.is_none() {
            return Err(CoreError::PreconditionFailed(format!(
                "kp20 expanded support pcurve authoring left half-edge {cursor:?} without pcurve",
            )));
        }
        cursor = half_edge.next;
        if cursor == start {
            return Ok(());
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 expanded support pcurve loop did not close".into(),
    ))
}
