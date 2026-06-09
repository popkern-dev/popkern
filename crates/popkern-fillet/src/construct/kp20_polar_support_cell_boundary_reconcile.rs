//! Segment-aware G0 reconciliation for KP20 support-cell shared full boundaries.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot;

#[path = "kp20_polar_support_cell_boundary_row_group.rs"]
mod boundary_row_group;
#[path = "kp20_polar_support_cell_group_reconcile.rs"]
mod group_reconcile;
#[path = "kp20_polar_support_cell_row_fit.rs"]
mod row_fit;
#[path = "kp20_polar_support_cell_segment_transfer.rs"]
mod segment_transfer;

#[derive(Debug, Clone, Copy)]
struct FullBoundaryTrim {
    face_a: FaceHandle,
    side_a: BoundarySide,
    start_a: f64,
    end_a: f64,
    face_b: FaceHandle,
    side_b: BoundarySide,
    start_b: f64,
    end_b: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum BoundarySide {
    UMin,
    UMax,
    VMin,
    VMax,
}

pub(crate) fn reconcile_support_cell_full_boundaries(
    brep: &mut BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<(), CoreError> {
    let trims = full_boundary_trims(brep, slots)?;
    for trim in &trims {
        reconcile_full_boundary_trim(brep, *trim)?;
    }
    group_reconcile::reconcile_grouped_boundary_rows(brep, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell boundary reconciliation produced invalid BRep: {error}"
        ))
    })
}

pub(crate) fn reconcile_assigned_support_cell_full_boundaries(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<(), CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    reconcile_support_cell_full_boundaries(brep, &slots)
}

pub(crate) fn reconcile_assigned_support_cell_worst_full_boundary(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<(), CoreError> {
    reconcile_assigned_support_cell_top_full_boundaries(brep, faces, 1)
}

pub(crate) fn reconcile_assigned_support_cell_top_full_boundaries(
    brep: &mut BRep,
    faces: &[FaceHandle],
    count: usize,
) -> Result<(), CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    let top_trims = top_boundary_trims(brep, &trims, count)?;
    for &trim in &top_trims {
        reconcile_full_boundary_trim(brep, trim)?;
    }
    group_reconcile::reconcile_grouped_boundary_rows(brep, &top_trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell targeted-boundary reconciliation produced invalid BRep: {error}"
        ))
    })
}

pub(crate) fn reconcile_assigned_support_cell_worst_boundary_row_group(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<(), CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    let row_group = boundary_row_group::worst_boundary_row_group_trims(brep, &trims)?;
    group_reconcile::reconcile_grouped_boundary_rows(brep, &row_group)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell row-group reconciliation produced invalid BRep: {error}"
        ))
    })
}

fn assigned_support_cell_slots(
    brep: &BRep,
    faces: &[FaceHandle],
) -> Vec<Kp20ExpandedSupportSurfaceSlot> {
    let mut slots = Vec::with_capacity(faces.len());
    for &face in faces {
        let Ok(half_edges) = face_outer_half_edges(brep, face) else {
            continue;
        };
        let boundary_edge_count = half_edges.len();
        let assignment = if boundary_edge_count == 4 {
            crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch
        } else {
            crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace
        };
        slots.push(Kp20ExpandedSupportSurfaceSlot {
            face,
            boundary_edge_count,
            assignment,
        });
    }
    slots
}

fn top_boundary_trims(
    brep: &BRep,
    trims: &[FullBoundaryTrim],
    count: usize,
) -> Result<Vec<FullBoundaryTrim>, CoreError> {
    let mut ranked = Vec::with_capacity(trims.len());
    for &trim in trims {
        let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
        let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
        let gap = trim_gap(&patch_a, &patch_b, trim);
        ranked.push((trim, gap));
    }
    ranked.sort_by(|(_, a), (_, b)| b.total_cmp(a));
    Ok(ranked
        .into_iter()
        .take(count)
        .map(|(trim, _)| trim)
        .collect())
}

fn full_boundary_trims(
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
                face_a: slot.face,
                side_a,
                start_a,
                end_a,
                face_b,
                side_b,
                start_b,
                end_b,
            });
        }
    }
    Ok(out)
}

fn reauthor_support_cell_shared_edge_curves(
    brep: &mut BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<(), CoreError> {
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
            reauthor_edge_curve_from_topology_vertices(brep, half_edge)?;
        }
    }
    Ok(())
}

fn reauthor_edge_curve_from_topology_vertices(
    brep: &mut BRep,
    half_edge: HalfEdgeHandle,
) -> Result<(), CoreError> {
    let edge = brep.half_edges[half_edge].edge;
    let canonical = brep.edges[edge].half_edge;
    let twin = brep.half_edges[canonical].twin;
    let start = brep.vertices[brep.half_edges[canonical].origin].position;
    let end = brep.vertices[brep.half_edges[twin].origin].position;
    let tolerance = brep.tolerance.linear().max(1.0e-12);
    let Some(curve) = CurveOnEdge::line_segment(start, end, tolerance) else {
        return Err(CoreError::PreconditionFailed(
            "kp20 support-cell boundary reconciliation found degenerate shared edge".into(),
        ));
    };
    brep.edges[edge].curve = Some(curve);
    Ok(())
}

fn reconcile_full_boundary_trim(brep: &mut BRep, trim: FullBoundaryTrim) -> Result<(), CoreError> {
    let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
    let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
    let before_gap = trim_gap(&patch_a, &patch_b, trim);
    if let Some((candidate_a, candidate_b)) = segment_transfer::exact_segment_transfer_candidate(
        &patch_a,
        trim.side_a,
        trim.start_a,
        trim.end_a,
        &patch_b,
        trim.side_b,
        trim.start_b,
        trim.end_b,
    ) {
        let after_gap = trim_gap(&candidate_a, &candidate_b, trim);
        if after_gap + 1.0e-12 < before_gap {
            set_rectangular_binonic_on_face(brep, trim.face_a, candidate_a)?;
            set_rectangular_binonic_on_face(brep, trim.face_b, candidate_b)?;
            return Ok(());
        }
    }
    let mut candidate_a = patch_a.clone();
    let mut candidate_b = patch_b.clone();
    let mut params_a = [0.0; 10];
    let mut params_b = [0.0; 10];
    let mut targets = [Point3::new(0.0, 0.0, 0.0); 10];
    for index in 0..10 {
        let fraction = index as f64 / 9.0;
        params_a[index] = lerp_scalar(trim.start_a, trim.end_a, fraction);
        params_b[index] = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
        targets[index] = average_point(
            row_fit::eval_boundary(&patch_a, trim.side_a, params_a[index]),
            row_fit::eval_boundary(&patch_b, trim.side_b, params_b[index]),
        );
    }
    row_fit::fit_boundary_samples(
        &mut candidate_a,
        trim.side_a,
        &boundary_fit_samples(&params_a, &targets),
    );
    row_fit::fit_boundary_samples(
        &mut candidate_b,
        trim.side_b,
        &boundary_fit_samples(&params_b, &targets),
    );
    let after_gap = trim_gap(&candidate_a, &candidate_b, trim);
    if after_gap + 1.0e-12 < before_gap
        && row_fit::max_boundary_displacement(&patch_a, &candidate_a, trim.side_a) <= 0.1
        && row_fit::max_boundary_displacement(&patch_b, &candidate_b, trim.side_b) <= 0.1
    {
        set_rectangular_binonic_on_face(brep, trim.face_a, candidate_a)?;
        set_rectangular_binonic_on_face(brep, trim.face_b, candidate_b)?;
    }
    Ok(())
}

fn boundary_fit_samples(
    params: &[f64; 10],
    targets: &[Point3; 10],
) -> [row_fit::BoundaryRowFitSample; 10] {
    std::array::from_fn(|index| row_fit::BoundaryRowFitSample {
        parameter: params[index],
        target: targets[index],
    })
}

fn rectangular_binonic_on_face(
    brep: &BRep,
    face: FaceHandle,
) -> Result<popkern_core::geometry::gregory::RectangularGregoryBiNonic, CoreError> {
    let Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(patch))) = brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return Err(CoreError::PreconditionFailed(
            "kp20 support-cell boundary reconciliation requires rectangular bi-nonic surfaces"
                .into(),
        ));
    };
    Ok((**patch).clone())
}

fn set_rectangular_binonic_on_face(
    brep: &mut BRep,
    face: FaceHandle,
    patch: popkern_core::geometry::gregory::RectangularGregoryBiNonic,
) -> Result<(), CoreError> {
    let Some(face_data) = brep.faces.get_mut(face) else {
        return Err(CoreError::InvalidHandle(
            "kp20 support-cell boundary reconciliation face invalid".into(),
        ));
    };
    face_data.surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(Box::new(patch))),
        orientation: true,
    });
    Ok(())
}

fn boundary_span(pcurve: &Pcurve) -> Option<(BoundarySide, f64, f64)> {
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

fn trim_gap(
    patch_a: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    patch_b: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    trim: FullBoundaryTrim,
) -> f64 {
    let mut worst = 0.0_f64;
    for index in 0..9 {
        let fraction = index as f64 / 8.0;
        let param_a = lerp_scalar(trim.start_a, trim.end_a, fraction);
        let param_b = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
        worst = worst.max(
            row_fit::eval_boundary(patch_a, trim.side_a, param_a)
                .distance_to(&row_fit::eval_boundary(patch_b, trim.side_b, param_b)),
        );
    }
    worst
}

fn average_point(a: Point3, b: Point3) -> Point3 {
    a + (b - a) * 0.5
}

fn lerp_scalar(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, CoreError> {
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
