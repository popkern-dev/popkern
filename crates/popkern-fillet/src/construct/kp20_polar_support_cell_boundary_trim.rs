//! Boundary-trim enumeration for KP20 support-cell patchwork.

use popkern_core::error::CoreError;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot;

#[derive(Debug, Clone, Copy)]
pub(crate) struct FullBoundaryTrim {
    pub(super) edge: EdgeHandle,
    pub(super) half_edge_a: HalfEdgeHandle,
    pub(super) face_a: FaceHandle,
    pub(super) side_a: BoundarySide,
    pub(super) start_a: f64,
    pub(super) end_a: f64,
    pub(super) half_edge_b: HalfEdgeHandle,
    pub(super) face_b: FaceHandle,
    pub(super) side_b: BoundarySide,
    pub(super) start_b: f64,
    pub(super) end_b: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum BoundarySide {
    UMin,
    UMax,
    VMin,
    VMax,
}

pub(crate) fn full_boundary_trims(
    brep: &BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<Vec<FullBoundaryTrim>, CoreError> {
    let mut out = Vec::new();
    for (slot_a, slot) in slots.iter().enumerate() {
        for &half_edge in &face_outer_half_edges(brep, slot.face)? {
            let twin = brep.half_edges[half_edge].twin;
            let Some(face_b) = brep
                .loops
                .get(brep.half_edges[twin].loop_)
                .map(|loop_| loop_.face)
            else {
                continue;
            };
            let Some(slot_b) = slots.iter().position(|slot| slot.face == face_b) else {
                continue;
            };
            if slot_b <= slot_a {
                continue;
            }
            let Some(pcurve_a) = brep.half_edges[half_edge].pcurve.as_ref() else {
                continue;
            };
            let Some(pcurve_b) = brep.half_edges[twin].pcurve.as_ref() else {
                continue;
            };
            let Some((side_a, start_a, end_a)) = boundary_span(pcurve_a) else {
                continue;
            };
            let Some((side_b, start_b, end_b)) = boundary_span(pcurve_b) else {
                continue;
            };
            out.push(FullBoundaryTrim {
                edge: brep.half_edges[half_edge].edge,
                half_edge_a: half_edge,
                face_a: slot.face,
                side_a,
                start_a,
                end_a,
                half_edge_b: twin,
                face_b,
                side_b,
                start_b,
                end_b,
            });
        }
    }
    Ok(out)
}

/// Count of slot face-pair edges that could NOT become trims,
/// by reason — every skipped edge is a hole in the chain
/// network's coverage (its rows are never glued), so the counts
/// belong in the construction facts.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct SkippedBoundaryTrims {
    pub(crate) missing_pcurve: usize,
    pub(crate) unrecognized_boundary_span: usize,
}

/// Mirror of [`full_boundary_trims`]' enumeration that counts the
/// candidate edges its filters reject.
pub(crate) fn count_skipped_boundary_trims(
    brep: &BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<SkippedBoundaryTrims, CoreError> {
    let mut skipped = SkippedBoundaryTrims::default();
    for (slot_a, slot) in slots.iter().enumerate() {
        for &half_edge in &face_outer_half_edges(brep, slot.face)? {
            let twin = brep.half_edges[half_edge].twin;
            let Some(face_b) = brep
                .loops
                .get(brep.half_edges[twin].loop_)
                .map(|loop_| loop_.face)
            else {
                continue;
            };
            let Some(slot_b) = slots.iter().position(|slot| slot.face == face_b) else {
                continue;
            };
            if slot_b <= slot_a {
                continue;
            }
            let (Some(pcurve_a), Some(pcurve_b)) = (
                brep.half_edges[half_edge].pcurve.as_ref(),
                brep.half_edges[twin].pcurve.as_ref(),
            ) else {
                skipped.missing_pcurve += 1;
                continue;
            };
            if boundary_span(pcurve_a).is_none() || boundary_span(pcurve_b).is_none() {
                skipped.unrecognized_boundary_span += 1;
            }
        }
    }
    Ok(skipped)
}

pub(crate) fn face_outer_half_edges(
    brep: &BRep,
    face: FaceHandle,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "kp20 support-cell boundary reconciliation face invalid".into(),
            )
        })?
        .outer_loop;
    let start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| {
            CoreError::InvalidHandle(
                "kp20 support-cell boundary reconciliation loop invalid".into(),
            )
        })?
        .half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 support-cell boundary reconciliation loop did not close".into(),
    ))
}

pub(crate) fn boundary_span(pcurve: &Pcurve) -> Option<(BoundarySide, f64, f64)> {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    if close(start.x, 0.0) && close(end.x, 0.0) && usable_span(start.y, end.y) {
        Some((BoundarySide::UMin, start.y, end.y))
    } else if close(start.x, 1.0) && close(end.x, 1.0) && usable_span(start.y, end.y) {
        Some((BoundarySide::UMax, start.y, end.y))
    } else if close(start.y, 0.0) && close(end.y, 0.0) && usable_span(start.x, end.x) {
        Some((BoundarySide::VMin, start.x, end.x))
    } else if close(start.y, 1.0) && close(end.y, 1.0) && usable_span(start.x, end.x) {
        Some((BoundarySide::VMax, start.x, end.x))
    } else {
        None
    }
}

fn usable_span(a: f64, b: f64) -> bool {
    (a - b).abs() > 1.0e-9 && (0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}
