//! Mixed cap-boundary recovery from materialized target-cylinder replacements.

use super::*;

pub(super) fn materialized_cap_boundary_from_target_cylinder_replacements(
    brep: &mut BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    fillet_ends: &[FilletEndOwnership],
    cylinder_loop_replacements: &[PendingCylinderLoopReplacement],
    reuse_shared_cap_cylinder_seam: bool,
) -> Result<Option<PartitionBoundaryOutcome>, ApexBuildError> {
    if cylinder_loop_replacements.is_empty() {
        return Ok(None);
    }
    let replacements_by_fillet = cylinder_loop_replacements
        .iter()
        .map(|replacement| (replacement.fillet_index, replacement))
        .collect::<std::collections::BTreeMap<_, _>>();
    let mut cap_cycle = Vec::with_capacity(2 * fillet_ends.len());
    let mut fillet_end_cap_sides = std::collections::BTreeMap::new();
    let mut strip_side_overrides = std::collections::BTreeMap::new();
    for (fillet_end, face_edge) in fillet_ends.iter().zip(face_edges) {
        let (cap_side, opposite_side) = match replacements_by_fillet.get(&fillet_end.fillet_index) {
            Some(replacement)
                if reuse_shared_cap_cylinder_seam
                    && replacement.kind == CylinderLoopOrientationKind::CurrentSharedFilletEnd =>
            {
                // Strict V-free mode: the replacement loop keeps the
                // current shared FilletEnd cylinder side, so the FilletEnd
                // edge itself must be the cap-to-cylinder seam — the cap
                // cycle takes that side's live twin. Cloning a parallel
                // edge here leaves the cap and the re-trimmed cylinder on
                // two coincident but topologically disconnected edges
                // (Stage-D contract: exactly one cap side and one true
                // cylinder side). Default mode keeps the clone: its
                // opposite side is the strip side the V-spur apex cells
                // consume via the strip-loop overrides.
                (
                    reused_cap_side_from_shared_cylinder_seam(brep, fillet_end, replacement)?,
                    None,
                )
            }
            Some(replacement) => {
                cloned_cap_side_from_cylinder_replacement(brep, fillet_end, replacement)?
            }
            None => (
                live_cap_side_for_edge(brep, cap_face, fillet_end.edge, fillet_end.cap_half_edge)?,
                None,
            ),
        };
        let face_edge_cap_side = live_cap_side_for_edge(
            brep,
            cap_face,
            face_edge.edge,
            face_edge.face_edge_cap_half_edge,
        )?;
        cap_cycle.push(cap_side);
        cap_cycle.push(face_edge_cap_side);
        fillet_end_cap_sides.insert(fillet_end.fillet_index, cap_side);
        if let Some(opposite_side) = opposite_side {
            strip_side_overrides.insert(fillet_end.fillet_index, opposite_side);
        }
    }
    // Strict-mode junction bridging (PK-126, measured on cube_mixed
    // N=3): the reused FilletEnd cap side and the following face-edge
    // cap side can be separated by intervening live cap-loop
    // half-edges (split and pocket structures the N=3 execution
    // order leaves between them; N >= 4 has none). Bridge each
    // junction with the CURRENT `.next` chain so the reauthored
    // cycle is contiguous; an unbridgeable junction falls through to
    // the validator's original diagnostic.
    let cap_cycle = bridge_cap_cycle_junction_gaps(brep, cap_cycle);
    validate_half_edge_cycle(
        brep,
        &cap_cycle,
        "materialized target-cylinder cap boundary",
        cap_face,
    )?;
    Ok(Some(PartitionBoundaryOutcome::from_reauthored_cap_cycle(
        cap_cycle,
        fillet_end_cap_sides,
        strip_side_overrides,
    )))
}

/// Bridge non-adjacent junctions in a reauthored cap cycle — now a
/// delegation to the CENTRAL junction settlement (PK-129): the same
/// splice serves every planned cycle inside `wire_existing_loop`,
/// and this pre-wiring validation sees the settled cycle.
fn bridge_cap_cycle_junction_gaps(brep: &BRep, cycle: Vec<HalfEdgeHandle>) -> Vec<HalfEdgeHandle> {
    crate::apex_construction::repeated_apex_partition_topology::settle_planned_loop_junction_gaps(
        brep, &cycle,
    )
}

/// The live twin of the replacement loop's shared FilletEnd cylinder
/// side — the cap side of the reused cap/cylinder seam edge.
fn reused_cap_side_from_shared_cylinder_seam(
    brep: &BRep,
    fillet_end: &FilletEndOwnership,
    replacement: &PendingCylinderLoopReplacement,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let Some(shared_cylinder_side) = replacement.cycle.first().copied() else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot reuse the shared \
                 cap/cylinder seam: selected target-cylinder replacement loop is empty",
                fillet_end.fillet_index
            ),
        });
    };
    let cap_side = brep
        .half_edges
        .get(shared_cylinder_side)
        .map(|half_edge| half_edge.twin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot reuse the shared \
                 cap/cylinder seam: shared cylinder side {shared_cylinder_side:?} is stale",
                fillet_end.fillet_index
            ),
        })?;
    if !brep.half_edges.contains_key(cap_side) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot reuse the shared \
                 cap/cylinder seam: cap-side twin {cap_side:?} of {shared_cylinder_side:?} is \
                 dead",
                fillet_end.fillet_index
            ),
        });
    }
    Ok(cap_side)
}

fn cloned_cap_side_from_cylinder_replacement(
    brep: &mut BRep,
    fillet_end: &FilletEndOwnership,
    replacement: &PendingCylinderLoopReplacement,
) -> Result<(HalfEdgeHandle, Option<HalfEdgeHandle>), ApexBuildError> {
    let Some(shared_cylinder_side) = replacement.cycle.first().copied() else {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot derive materialized cap \
                 boundary: selected target-cylinder replacement loop is empty",
                fillet_end.fillet_index
            ),
        });
    };
    let cylinder_data = brep
        .half_edges
        .get(shared_cylinder_side)
        .cloned()
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot derive materialized cap \
                     boundary: selected shared cylinder side {shared_cylinder_side:?} is stale",
                fillet_end.fillet_index
            ),
        })?;
    let cylinder_end =
        support_split_half_edge_end_vertex(brep, shared_cylinder_side).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cannot derive materialized cap \
                 boundary: selected shared cylinder side {shared_cylinder_side:?} has no \
                 recoverable endpoint",
                    fillet_end.fillet_index
                ),
            }
        })?;
    let cap_start = cylinder_end;
    let cap_end = cylinder_data.origin;
    let curve = brep
        .edges
        .get(cylinder_data.edge)
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot derive materialized cap \
                 boundary: selected shared cylinder side {shared_cylinder_side:?} has no exact \
                 curve to clone",
                fillet_end.fillet_index
            ),
        })?;
    let (edge, cap_side, opposite_side) = insert_unwired_edge_pair(brep, cap_start, cap_end);
    brep.edges[edge].curve = Some(curve);
    if half_edge_end_vertex(brep, cap_side) != Some(cap_end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot derive materialized cap \
                 boundary: cloned cap side {cap_side:?} does not span {:?}->{:?}",
                fillet_end.fillet_index, cap_start, cap_end
            ),
        });
    }
    Ok((cap_side, Some(opposite_side)))
}

fn live_cap_side_for_edge(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
    recorded: HalfEdgeHandle,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    if half_edge_is_on_face(brep, recorded, cap_face) {
        return Ok(recorded);
    }
    let cap_loop = brep
        .faces
        .get(cap_face)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot refresh cap boundary: cap face \
                 {cap_face:?} is stale"
            ),
        })?
        .outer_loop;
    for half_edge in collect_loop_half_edges(brep, cap_loop)? {
        if brep
            .half_edges
            .get(half_edge)
            .is_some_and(|data| data.edge == edge)
        {
            return Ok(half_edge);
        }
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition cannot refresh cap boundary edge {edge:?}: no live \
             cap half-edge on {cap_face:?}"
        ),
    })
}

fn half_edge_is_on_face(brep: &BRep, half_edge: HalfEdgeHandle, face: FaceHandle) -> bool {
    brep.half_edges
        .get(half_edge)
        .and_then(|data| brep.loops.get(data.loop_))
        .is_some_and(|loop_data| loop_data.face == face)
}
