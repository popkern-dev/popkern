//! Shared directional skin ladder: classify → offset surfaces → clone
//! skin → recompute vertices/edges (+ collapse gate) → apply geometry.
//!
//! One ladder produces one moved skin. `thicken` runs two (outward +
//! inward over the whole body); midplane shells run two with removed
//! faces; the classic inward/outward shell runs the equivalent stages
//! inline in `shell_impl`.

use std::collections::HashMap;

use crate::classify::classify_faces;
use crate::edge_recompute::recompute_edges_with_skips;
use crate::error::{ShellError, ShellWarning};
use crate::inner_skin::{apply_offset_geometry, build_inner_skin};
use crate::join_classification::{classify_joins, intersection_gap_postcheck};
use crate::skin_map::{JoinRecord, SkinMap};
use crate::surface_offset::{OffsetResult, compute_offset_surfaces};
use crate::types::{JoinStyle, ShellDirection};
use crate::vertex_recompute::recompute_vertices_with_gap_corners;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// What the join stage of one skin ladder produced.
#[derive(Debug, Clone, Default)]
pub(crate) struct JoinOutcome {
    /// Join faces constructed on this skin.
    pub records: Vec<JoinRecord>,
    /// Containment warnings (style downgrades, sharp cone apexes).
    pub warnings: Vec<ShellWarning>,
    /// Rim-corner splits: (new SKIN vertex, OUTER body rim vertex).
    pub rim_splits: Vec<(
        popkern_core::topology::VertexHandle,
        popkern_core::topology::VertexHandle,
    )>,
}

/// Stages one ladder run emits (for PK-276 progress totals).
pub(crate) const SKIN_LADDER_STAGE_COUNT: usize = 5;

/// Run one directional skin ladder.
///
/// `removed` lists faces excluded from offsetting; when
/// `delete_removed_from_skin` is true they are deleted from the skin
/// (the void-skin shape the classic shell pipeline assembles against),
/// otherwise they stay present-but-unclassified so a later assembly can
/// discover the opening rims on THIS skin (the midplane outer skin).
///
/// `requested_scale` converts the ladder's per-side distance into the
/// caller's parameter units in `OffsetTooLarge` declines (2.0 when the
/// caller speaks total thickness and the ladder runs half per side).
#[allow(clippy::too_many_arguments)]
pub(crate) fn offset_skin_ladder(
    brep: &BRep,
    distance: f64,
    face_distances: &[(FaceHandle, f64)],
    direction: ShellDirection,
    removed: &[FaceHandle],
    delete_removed_from_skin: bool,
    allow_partial: bool,
    join_style: JoinStyle,
    serve_rim_corners: bool,
    stage_tol: &Tolerance,
    requested_scale: f64,
    labels: [&'static str; SKIN_LADDER_STAGE_COUNT],
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<(BRep, SkinMap, OffsetResult, JoinOutcome), ShellError> {
    use std::time::Instant;
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(ShellError::Cancelled);
            }
        };
    }
    stage_gate!();
    let t0 = Instant::now();
    let classification = classify_faces(brep, distance, face_distances, direction, removed)?;
    tracker.finish_stage(labels[0], t0, brep.faces.len());

    stage_gate!();
    let t0 = Instant::now();
    let offset = compute_offset_surfaces(brep, &classification, allow_partial, stage_tol)?;
    tracker.finish_stage(labels[1], t0, offset.surfaces.len());

    stage_gate!();
    let t0 = Instant::now();
    let deleted: Vec<FaceHandle> = if delete_removed_from_skin {
        let mut all = classification.removed_faces.clone();
        all.extend_from_slice(&offset.auto_removed);
        all
    } else {
        Vec::new()
    };
    let (mut skin, skin_map) = build_inner_skin(brep, &deleted)?;
    tracker.finish_stage(labels[2], t0, skin.faces.len());

    stage_gate!();
    let t0 = Instant::now();
    // Face-collapse resolution (PK-326): features thinner than the
    // offset dissolve out of the skin before the join plan is made.
    let collapse_warnings = crate::collapse_resolve::resolve_face_collapse(
        brep,
        &mut skin,
        &skin_map,
        &offset.surfaces,
        stage_tol,
    )?;
    // Join classification (PK-325): measured per-edge gap/smooth/
    // resolvable classes against the offset surfaces, planned before
    // any recompute so gap corners take their retained-run positions
    // and gap edges keep their clone curves for canal anchoring.
    let join_plan = classify_joins(
        &skin,
        &skin_map,
        &offset.surfaces,
        stage_tol,
        join_style,
        serve_rim_corners,
    )?;
    recompute_vertices_with_gap_corners(
        &mut skin,
        &offset.surfaces,
        &skin_map,
        brep,
        stage_tol,
        &join_plan.gap_corners,
    )?;
    let undone = recompute_edges_with_skips(
        &mut skin,
        &offset.surfaces,
        &skin_map,
        brep,
        stage_tol,
        &join_plan.arc_skip_edges(),
    )?;
    intersection_gap_postcheck(&join_plan, &undone)?;
    // Collapse gate (part of the recompute stage): a distance past the
    // local feature size inverts this skin. Declines report the
    // caller's parameter units via `requested_scale`.
    crate::collapse_check::check_skin_collapse(brep, &skin, &skin_map, distance, stage_tol)
        .map_err(|e| match e {
            ShellError::OffsetTooLarge {
                requested,
                max_safe,
                binding,
            } => ShellError::OffsetTooLarge {
                requested: requested * requested_scale,
                max_safe: max_safe * requested_scale,
                binding,
            },
            other => other,
        })?;
    tracker.finish_stage(labels[3], t0, skin.edges.len());

    stage_gate!();
    let t0 = Instant::now();
    apply_offset_geometry(&mut skin, &offset.surfaces, &skin_map)?;
    // Join construction (PK-325): Euler surgery inserts the rolling-
    // ball strips and sphere corner windows on the recomputed skin.
    let arc = crate::arc_join::construct_arc_joins(
        &mut skin,
        &join_plan,
        &offset.surfaces,
        &skin_map,
        brep,
        stage_tol,
    )?;
    tracker.finish_stage(labels[4], t0, skin.faces.len());

    let mut warnings = collapse_warnings;
    warnings.extend(join_plan.warnings);
    crate::arc_join_apex::retire_served_apex_warnings(&mut warnings, &arc.served_apexes);
    let outcome = JoinOutcome {
        records: arc.records,
        warnings,
        rim_splits: arc.rim_splits,
    };
    Ok((skin, skin_map, offset, outcome))
}

/// Combine two original→skin maps into one skinA→skinB map.
pub(crate) fn combine_skin_maps(map_a: &SkinMap, map_b: &SkinMap) -> SkinMap {
    let reverse_a_vertices = map_a.reverse_vertex_map();
    let reverse_a_edges = map_a.reverse_edge_map();
    let reverse_a_he: HashMap<_, _> = map_a.half_edges.iter().map(|(&k, &v)| (v, k)).collect();
    let reverse_a_faces = map_a.reverse_face_map();

    let mut vertices = HashMap::new();
    let mut edges = HashMap::new();
    let mut half_edges = HashMap::new();
    let mut faces = HashMap::new();

    for (&a_vh, &orig_vh) in &reverse_a_vertices {
        if let Some(&b_vh) = map_b.vertices.get(&orig_vh) {
            vertices.insert(a_vh, b_vh);
        }
    }
    for (&a_eh, &orig_eh) in &reverse_a_edges {
        if let Some(&b_eh) = map_b.edges.get(&orig_eh) {
            edges.insert(a_eh, b_eh);
        }
    }
    for (&a_heh, &orig_heh) in &reverse_a_he {
        if let Some(&b_heh) = map_b.half_edges.get(&orig_heh) {
            half_edges.insert(a_heh, b_heh);
        }
    }
    for (&a_fh, &orig_fh) in &reverse_a_faces {
        if let Some(&b_fh) = map_b.faces.get(&orig_fh) {
            faces.insert(a_fh, b_fh);
        }
    }

    SkinMap {
        vertices,
        edges,
        half_edges,
        faces,
    }
}
