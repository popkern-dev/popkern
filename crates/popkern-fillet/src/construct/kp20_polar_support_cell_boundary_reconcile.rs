//! Segment-aware G0 reconciliation for KP20 support-cell shared full boundaries.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot;

#[path = "kp20_polar_support_cell_band_derive.rs"]
mod band_derive;
#[path = "kp20_polar_support_cell_band_field.rs"]
mod band_field;
#[path = "kp20_polar_support_cell_leg_class.rs"]
mod leg_class;
#[path = "kp20_polar_support_cell_leg_direction.rs"]
mod leg_direction;
#[path = "kp20_polar_support_cell_surface_access.rs"]
mod surface_access;
pub(self) use surface_access::{
    rectangular_binonic_on_face, set_rectangular_binonic_dual_twist_on_face,
    set_rectangular_binonic_on_face,
};
#[path = "kp20_polar_support_cell_band_predictor.rs"]
mod band_predictor;
#[path = "kp20_polar_support_cell_boundary_cluster.rs"]
mod boundary_cluster;
#[path = "kp20_polar_support_cell_boundary_composed.rs"]
pub(crate) mod boundary_composed;
#[path = "kp20_polar_support_cell_boundary_prefix.rs"]
mod boundary_prefix;
#[path = "kp20_polar_support_cell_boundary_row_group.rs"]
mod boundary_row_group;
#[path = "kp20_polar_support_cell_boundary_trim.rs"]
mod boundary_trim;
#[path = "kp20_polar_support_cell_cell_reauthor.rs"]
mod cell_reauthor;
#[path = "kp20_polar_support_cell_chain_discovery.rs"]
mod chain_discovery;
#[path = "kp20_polar_support_cell_chain_fit.rs"]
mod chain_fit;
#[path = "kp20_polar_support_cell_chain_nodes.rs"]
mod chain_nodes;
#[path = "kp20_polar_support_cell_chain_ribbon.rs"]
mod chain_ribbon;
#[path = "kp20_polar_support_cell_chain_ribbon_field.rs"]
mod chain_ribbon_field;
#[path = "kp20_polar_support_cell_chain_solver.rs"]
mod chain_solver;
#[path = "kp20_polar_support_cell_chain_spacing.rs"]
mod chain_spacing;
#[path = "kp20_polar_support_cell_exact_boundary_row.rs"]
mod exact_boundary_row;
#[path = "kp20_polar_support_cell_group_reconcile.rs"]
mod group_reconcile;
#[path = "kp20_polar_support_cell_row_fit.rs"]
mod row_fit;
#[path = "kp20_polar_support_cell_row_segment.rs"]
mod row_segment;
#[path = "kp20_polar_support_cell_segment_topology.rs"]
mod segment_topology;
#[path = "kp20_polar_support_cell_segment_transfer.rs"]
mod segment_transfer;
#[path = "kp20_polar_support_cell_tangent_row.rs"]
mod tangent_row;
#[path = "kp20_polar_support_cell_trim_ownership.rs"]
mod trim_ownership;
#[path = "kp20_polar_support_cell_twist_resolution.rs"]
mod twist_resolution;

pub(super) use boundary_trim::{
    BoundarySide, FullBoundaryTrim, count_skipped_boundary_trims, face_outer_half_edges,
    full_boundary_trims,
};

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

/// Canonical chain-curve network: author one degree-9 curve per
/// discoverable boundary sub-chain and assign every member row
/// its exact restriction, so all trims of the chain (full-span,
/// one-sided, and both-partial) become CP identities. Declined
/// chains keep their current rows; the status string carries the
/// authored/declined counts and first decline reason.
pub(crate) fn reconcile_assigned_chain_curve_network_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let status = chain_fit::author_chain_curve_network(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell chain-curve network produced invalid BRep: {error}"
        ))
    })?;
    Ok(status)
}

/// Chain-curve network plus shared cross-tangent ribbons: chain
/// curves are tangent-plane-constrained at every pinned node
/// normal and each member row's row-1 controls are projected
/// onto the shared plane field — exact G0 by shared curves,
/// measured G1 by shared tangent planes.
pub(crate) fn reconcile_assigned_chain_network_cross_tangent_ribbon_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let status = chain_fit::author_chain_network_with_cross_tangent_ribbons(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell chain ribbon network produced invalid BRep: {error}"
        ))
    })?;
    Ok(status)
}

/// Chain network plus constructive cell re-authoring: row-1
/// bands DERIVED from the shared chain curves and transverse
/// ribbon fields (linear blends matched to the exact end legs),
/// twist conflicts averaged, interiors Coons-refilled — exact
/// G0 by shared curves, G1 by shared tangent-plane fields
/// rather than projection of disagreeing geometry.
pub(crate) fn reconcile_assigned_chain_constructive_cells_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    // Two fixed-point passes: the first pass's node normals are
    // averaged from the incoherent pre-construction sheets, so
    // the two boundaries' twist demands at each corner disagree
    // by the sheet mismatch (the bi-9 twist-compatibility
    // residual). After one pass the geometry is coherent;
    // recomputing the normals and re-deriving collapses the
    // demand disagreement and with it the corner G1 residual.
    let first_pass_status = chain_fit::author_chain_network_with_constructive_cells(brep, &slots)?;
    let status = chain_fit::author_chain_network_with_constructive_cells(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell constructive cell re-authoring produced invalid BRep: {error}"
        ))
    })?;
    let first_pass_twist = first_pass_status
        .split("max_twist_conflict_distance=")
        .nth(1)
        .and_then(|tail| tail.split(';').next())
        .unwrap_or("unknown");
    Ok(format!(
        "{status};constructive_pass_count=2;first_pass_max_twist_conflict_distance={first_pass_twist}"
    ))
}

/// Two-pass constructive re-authoring on the DUAL-TWIST patch
/// substrate: identical fixed-point structure to
/// [`reconcile_assigned_chain_constructive_cells_status`], but
/// each boundary's corner band coefficient lands in its own
/// rational candidate slot, so the corner twist conflicts the
/// plain bi-9 had to average cease to exist.
pub(crate) fn reconcile_assigned_chain_constructive_dual_twist_cells_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<String, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let first_pass_status =
        chain_fit::author_chain_network_with_constructive_dual_twist_cells(brep, &slots)?;
    let status = chain_fit::author_chain_network_with_constructive_dual_twist_cells(brep, &slots)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell dual-twist constructive re-authoring produced invalid BRep: {error}"
        ))
    })?;
    let first_pass_twist = first_pass_status
        .split("max_dual_twist_pair_separation=")
        .nth(1)
        .and_then(|tail| tail.split(';').next())
        .unwrap_or("unknown");
    Ok(format!(
        "{status};constructive_pass_count=2;first_pass_max_dual_twist_pair_separation={first_pass_twist}"
    ))
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

pub(crate) fn reconcile_assigned_support_cell_ranked_boundary_prefix(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<(), CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    boundary_prefix::reconcile_ranked_boundary_prefix(brep, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell ranked-prefix reconciliation produced invalid BRep: {error}"
        ))
    })
}

pub(crate) fn reconcile_assigned_support_cell_worst_boundary_neighborhood(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<(), CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    boundary_cluster::reconcile_worst_boundary_neighborhood(brep, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell worst-neighborhood reconciliation produced invalid BRep: {error}"
        ))
    })
}

pub(crate) fn reconcile_assigned_support_cell_worst_endpoint_cluster_with_status(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<&'static str, CoreError> {
    let slots = assigned_support_cell_slots(brep, faces);
    let trims = full_boundary_trims(brep, &slots)?;
    let status = boundary_cluster::reconcile_worst_endpoint_cluster(brep, &trims)?;
    reauthor_support_cell_shared_edge_curves(brep, &slots)?;
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell worst-endpoint cluster reconciliation produced invalid BRep: {error}"
        ))
    })?;
    Ok(status)
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
