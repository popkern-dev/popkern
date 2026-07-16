use super::partition_spilled_corner::normalize_spilled_fillet_corners;
use super::*;

#[path = "repeated_apex_partition_support_loop_residual.rs"]
mod support_loop_residual;

pub(super) use support_loop_residual::materialize_authored_support_boundary_residual_cycles;
use support_loop_residual::materialize_support_boundary_residual_cycles;

pub(super) fn materialize_closed_support_trim_cycles(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
    route_face_edges_to_coplanar_facets: bool,
) -> Result<(), ApexBuildError> {
    let rollback = brep.clone();
    let result = materialize_closed_support_trim_cycles_candidate(
        brep,
        face_edges,
        trim_sides,
        route_face_edges_to_coplanar_facets,
    );
    if let Err(err) = result {
        *brep = rollback;
        return Err(err);
    }
    validate_live_half_edge_cycles_after_support_rewrite(brep)?;
    Ok(())
}

fn validate_live_half_edge_cycles_after_support_rewrite(brep: &BRep) -> Result<(), ApexBuildError> {
    for (loop_handle, loop_data) in &brep.loops {
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut visited = std::collections::BTreeSet::new();
        for _ in 0..=brep.half_edges.len() {
            let half_edge =
                brep.half_edges
                    .get(cursor)
                    .ok_or_else(|| ApexBuildError::PlanMismatch {
                        reason: format!(
                            "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: missing half-edge {cursor:?}"
                        ),
                    })?;
            if half_edge.loop_ != loop_handle {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: half-edge {cursor:?} belongs to \
                         {:?}",
                        half_edge.loop_
                    ),
                });
            }
            if !visited.insert(cursor) {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition support ownership rewrite produced \
                         invalid live loop {loop_handle:?}: repeated half-edge {cursor:?}"
                    ),
                });
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
        if cursor != start {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition support ownership rewrite produced invalid \
                     live loop {loop_handle:?}: cycle did not close from {start:?}"
                ),
            });
        }
    }
    Ok(())
}

pub(super) fn repair_support_source_loops_after_cylinder_rewire(
    brep: &mut BRep,
    trim_sides: &[SupportOwnedTrimSide],
) -> Result<Vec<SupportSourceLoopTransfer>, ApexBuildError> {
    let mut transfers = Vec::new();
    for trim_side in trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        let Some(source) = trim_side.support_source else {
            continue;
        };
        if !brep.loops.contains_key(source.loop_) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot repair \
                     support loop after cylinder trim ownership transfer: stale source loop \
                     {loop_handle:?}",
                    fillet_index = trim_side.fillet_index,
                    loop_handle = source.loop_
                ),
            });
        }
        let mut source = source;
        if !brep.half_edges.contains_key(source.prev)
            || !brep.half_edges.contains_key(source.next)
            || !brep.half_edges.contains_key(trim_side.half_edge)
        {
            let Some(moved) = brep.half_edges.get(trim_side.half_edge) else {
                continue;
            };
            if moved.loop_ != source.loop_ {
                continue;
            }
            source.prev = moved.prev;
            source.next = moved.next;
        }
        let prev_still_points_to_moved = brep.half_edges[source.prev].next == trim_side.half_edge;
        let next_still_points_to_moved = brep.half_edges[source.next].prev == trim_side.half_edge;
        if !prev_still_points_to_moved || !next_still_points_to_moved {
            if !support_loop_contains_half_edge(brep, source.loop_, trim_side.half_edge) {
                continue;
            }
            let moved = &brep.half_edges[trim_side.half_edge];
            if moved.loop_ != source.loop_ {
                continue;
            }
            source.prev = moved.prev;
            source.next = moved.next;
        }
        let replacement = author_support_source_loop_replacement_span(brep, trim_side, source)?;
        brep.half_edges[source.prev].next = replacement;
        brep.half_edges[source.next].prev = replacement;
        brep.half_edges[replacement].prev = source.prev;
        brep.half_edges[replacement].next = source.next;
        brep.half_edges[replacement].loop_ = source.loop_;
        if brep.loops[source.loop_].half_edge == trim_side.half_edge {
            brep.loops[source.loop_].half_edge = replacement;
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/support-loop-repair] FilletEnd[{}] replaced stolen cylinder \
                 side {:?} in source support loop {:?} with boundary span {:?} between {:?}->{:?}",
                trim_side.fillet_index,
                trim_side.half_edge,
                source.loop_,
                replacement,
                source.prev,
                source.next
            );
        }
        transfers.push(SupportSourceLoopTransfer {
            loop_: source.loop_,
            replacement,
        });
    }
    Ok(transfers)
}

pub(super) struct SupportSourceLoopTransfer {
    pub(super) loop_: LoopHandle,
    pub(super) replacement: HalfEdgeHandle,
}

/// Remove a side from its current loop's chain before another loop
/// adopts it: predecessor and successor join, and the loop re-anchors
/// if the side was its anchor. The donor is retiring scaffolding.
fn splice_side_out_of_donor_loop(brep: &mut BRep, side: HalfEdgeHandle) {
    let Some(data) = brep.half_edges.get(side).cloned() else {
        return;
    };
    if data.prev != side
        && data.next != side
        && brep.half_edges.contains_key(data.prev)
        && brep.half_edges.contains_key(data.next)
    {
        brep.half_edges[data.prev].next = data.next;
        brep.half_edges[data.next].prev = data.prev;
    }
    if brep
        .loops
        .get(data.loop_)
        .is_some_and(|loop_data| loop_data.half_edge == side)
        && data.next != side
    {
        brep.loops[data.loop_].half_edge = data.next;
    }
}

fn face_edge_endpoints_fit_face_plane(brep: &BRep, face: FaceHandle, side: HalfEdgeHandle) -> bool {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    if !matches!(surface.geometry, SurfaceGeometry::Plane { .. }) {
        return false;
    }
    let Some(data) = brep.half_edges.get(side) else {
        return false;
    };
    let Some(end) = brep
        .half_edges
        .get(data.twin)
        .map(|twin| twin.origin)
        .or_else(|| brep.half_edges.get(data.next).map(|next| next.origin))
    else {
        return false;
    };
    let fit = brep.tolerance.linear().max(1.0e-9) * 10.0;
    [data.origin, end].into_iter().all(|vertex| {
        brep.vertices.get(vertex).is_some_and(|v| {
            surface
                .geometry
                .closest_point(&v.position)
                .map(|(u, w)| surface.geometry.point_at(u, w).distance_to(&v.position) <= fit)
                .unwrap_or(false)
        })
    })
}

/// Explicit extended-facet cycle: [genA-twin, top, genB-twin,
/// span(anchor->old-tangent, the base side's displaced twin),
/// kept base-contact, span] — every side identified semantically; no
/// chain search. Returns Ok(None) for groups that are not the
/// two-generator facet shape.
fn assemble_extended_facet_cycle(
    brep: &BRep,
    support_face: FaceHandle,
    group: &[HalfEdgeHandle],
) -> Result<Option<Vec<HalfEdgeHandle>>, ApexBuildError> {
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace §14.9/assembly] {support_face:?}: group {} side(s)",
            group.len()
        );
    }
    if group.len() != 3 {
        return Ok(None);
    }
    let end_of = |he: HalfEdgeHandle| -> Option<VertexHandle> {
        let data = brep.half_edges.get(he)?;
        brep.half_edges.get(data.twin).map(|twin| twin.origin)
    };
    // Classify: the top side's endpoints are shared with the other two
    // sides' T-endpoints; generators' other endpoints are the anchors.
    let mut top = None;
    for (index, side) in group.iter().enumerate() {
        let (Some(start), Some(end)) = (
            brep.half_edges.get(*side).map(|data| data.origin),
            end_of(*side),
        ) else {
            return Ok(None);
        };
        let others: Vec<HalfEdgeHandle> = group
            .iter()
            .enumerate()
            .filter_map(|(other_index, other)| (other_index != index).then_some(*other))
            .collect();
        let touches = |vertex: VertexHandle| {
            others.iter().any(|other| {
                brep.half_edges
                    .get(*other)
                    .is_some_and(|data| data.origin == vertex)
                    || end_of(*other) == Some(vertex)
            })
        };
        if touches(start) && touches(end) {
            top = Some((*side, start, end, others));
            break;
        }
    }
    let Some((top_side, top_start, top_end, generators)) = top else {
        trace_assembly_decline(support_face, "no top side classified");
        return Ok(None);
    };
    // Order generators: genA ends at top_start; genB starts at top_end.
    let gen_a = generators
        .iter()
        .copied()
        .find(|side| end_of(*side) == Some(top_start));
    let gen_b = generators.iter().copied().find(|side| {
        brep.half_edges
            .get(*side)
            .is_some_and(|d| d.origin == top_end)
    });
    let (Some(gen_a), Some(gen_b)) = (gen_a, gen_b) else {
        trace_assembly_decline(support_face, "generator orientation mismatch");
        return Ok(None);
    };
    let anchor_a = brep
        .half_edges
        .get(gen_a)
        .map(|data| data.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: "explicit facet assembly lost generator A".to_owned(),
        })?;
    let anchor_b = end_of(gen_b).ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: "explicit facet assembly lost generator B".to_owned(),
    })?;
    // Spans: the displaced twins of the base's tangent spans, one out of
    // each anchor, whose mutual twin lives in a production Plane face
    // other than this facet.
    let facet_vertices: std::collections::BTreeSet<VertexHandle> = brep
        .faces
        .get(support_face)
        .and_then(|face| collect_loop_half_edges(brep, face.outer_loop).ok())
        .map(|cycle| {
            cycle
                .iter()
                .filter_map(|side| brep.half_edges.get(*side).map(|data| data.origin))
                .collect()
        })
        .unwrap_or_default();
    // The facet's corner spans are the TWINS of the base-loop tangent
    // spans touching each anchor (the base occupies one direction; the
    // displaced twin is the facet's side of the seam).
    let span_for = |anchor: VertexHandle, want_origin_at_anchor: bool| -> Option<HalfEdgeHandle> {
        brep.half_edges.iter().find_map(|(base_side, data)| {
            let twin = brep.half_edges.get(data.twin)?;
            if twin.twin != base_side {
                return None;
            }
            let base_face = brep.loops.get(data.loop_)?.face;
            if base_face == support_face {
                return None;
            }
            let plane_production = brep.faces.get(base_face).is_some_and(|face| {
                face.role != TopologyRole::ConstructionResidual
                    && matches!(
                        face.surface.as_ref().map(|surface| &surface.geometry),
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. })
                    )
            });
            if !plane_production {
                return None;
            }
            let base_start = data.origin;
            let base_end = twin.origin;
            let candidate = data.twin;
            if want_origin_at_anchor {
                // Facet span anchor -> tangent: twin of base side running
                // tangent -> anchor.
                (base_end == anchor && facet_vertices.contains(&base_start)).then_some(candidate)
            } else {
                // Facet span tangent -> anchor: twin of base side running
                // anchor -> tangent.
                (base_start == anchor && facet_vertices.contains(&base_end)).then_some(candidate)
            }
        })
    };
    let span_b = span_for(anchor_b, true);
    let span_a = span_for(anchor_a, false);
    let (Some(span_b), Some(span_a)) = (span_b, span_a) else {
        trace_assembly_decline(support_face, "no base tangent-span twins at the anchors");
        return Ok(None);
    };
    let old_tangent_b = end_of(span_b).ok_or_else(|| ApexBuildError::PlanMismatch {
        reason: "explicit facet assembly lost span B endpoint".to_owned(),
    })?;
    let old_tangent_a = brep
        .half_edges
        .get(span_a)
        .map(|data| data.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: "explicit facet assembly lost span A origin".to_owned(),
        })?;
    // Kept base-contact: current facet-loop side old_tangent_b -> old_tangent_a.
    let facet_loop = brep.faces.get(support_face).map(|face| face.outer_loop);
    let base_contact = facet_loop.and_then(|loop_handle| {
        collect_loop_half_edges(brep, loop_handle)
            .ok()?
            .into_iter()
            .find(|side| {
                brep.half_edges
                    .get(*side)
                    .is_some_and(|data| data.origin == old_tangent_b)
                    && end_of(*side) == Some(old_tangent_a)
            })
    });
    let Some(base_contact) = base_contact else {
        trace_assembly_decline(support_face, "no kept base-contact between old tangents");
        return Ok(None);
    };
    Ok(Some(vec![
        gen_a,
        top_side,
        gen_b,
        span_b,
        base_contact,
        span_a,
    ]))
}

fn trace_assembly_decline(support_face: FaceHandle, reason: &str) {
    if crate::direct_apex_trace::enabled() {
        eprintln!("[direct-apex-trace §14.9/assembly] {support_face:?} declined: {reason}");
    }
}

fn support_loop_contains_half_edge(
    brep: &BRep,
    loop_handle: LoopHandle,
    half_edge: HalfEdgeHandle,
) -> bool {
    collect_loop_half_edges(brep, loop_handle)
        .ok()
        .is_some_and(|cycle| cycle.contains(&half_edge))
}

fn author_support_source_loop_replacement_span(
    brep: &mut BRep,
    trim_side: &SupportOwnedTrimSide,
    source: SupportOwnedTrimSource,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let start = brep
        .half_edges
        .get(trim_side.half_edge)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: stale moved half-edge",
                trim_side.fillet_index, trim_side.half_edge
            ),
        })?;
    let end = brep
        .half_edges
        .get(source.next)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: stale next half-edge {:?}",
                trim_side.fillet_index, trim_side.half_edge, source.next
            ),
        })?;
    let curve = brep
        .half_edges
        .get(trim_side.half_edge)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot replace stolen support \
                 loop trim {:?}: moved edge has no curve",
                trim_side.fillet_index, trim_side.half_edge
            ),
        })?;
    let (edge, replacement) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].curve = Some(curve);
    brep.edges[edge].role = TopologyRole::Production;
    if support_split_half_edge_end_vertex(brep, replacement) != Some(end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] authored support replacement \
                 span with wrong endpoint: expected {start:?}->{end:?}, got {replacement:?}",
                trim_side.fillet_index
            ),
        });
    }
    Ok(replacement)
}

fn materialize_closed_support_trim_cycles_candidate(
    brep: &mut BRep,
    face_edges: &[FaceEdgeOwnership],
    trim_sides: &[SupportOwnedTrimSide],
    route_face_edges_to_coplanar_facets: bool,
) -> Result<(), ApexBuildError> {
    trace_topology_counts(brep, "support-loop-entry");
    let source_loop_transfers =
        repair_support_source_loops_after_cylinder_rewire(brep, trim_sides)?;
    let mut authored_boundary_spans = source_loop_transfers
        .iter()
        .map(|transfer| transfer.replacement)
        .collect::<Vec<_>>();
    let mut by_support_face: std::collections::BTreeMap<FaceHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    let mut generated_trim_count_by_support_face =
        std::collections::BTreeMap::<FaceHandle, usize>::new();
    let mut support_faces_with_existing_boundary = std::collections::BTreeSet::new();
    for trim_side in trim_sides {
        if !trim_side.requires_support_ownership {
            continue;
        }
        if brep
            .half_edges
            .get(trim_side.half_edge)
            .and_then(|half_edge| brep.loops.get(half_edge.loop_))
            .is_none_or(|loop_data| loop_data.face != trim_side.cylinder_face)
        {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot plan \
                     support ownership for {cylinder_half_edge:?}: it is not owned by target \
                     cylinder face {cylinder_face:?}",
                    fillet_index = trim_side.fillet_index,
                    cylinder_half_edge = trim_side.half_edge,
                    cylinder_face = trim_side.cylinder_face
                ),
            });
        }
        let twin = brep
            .half_edges
            .get(trim_side.half_edge)
            .map(|half_edge| half_edge.twin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{fillet_index}] cannot plan \
                     support ownership for stale half-edge {cylinder_half_edge:?}",
                    fillet_index = trim_side.fillet_index,
                    cylinder_half_edge = trim_side.half_edge
                ),
            })?;
        if trim_side.support_side_owned {
            continue;
        }
        by_support_face
            .entry(trim_side.support_face)
            .or_default()
            .push(twin);
        *generated_trim_count_by_support_face
            .entry(trim_side.support_face)
            .or_default() += 1;
    }
    for face_edge in face_edges {
        // Strict V-free routing: the FaceEdge's neighbor side belongs to
        // the support group whose facet plane carries both endpoints —
        // the retiring K-P fragment it currently sits in is merely
        // coplanar scaffolding. Index arithmetic cannot express this
        // pairing (K-P corners and facets adjoin fillets differently).
        let group_key = if route_face_edges_to_coplanar_facets {
            by_support_face
                .keys()
                .copied()
                .find(|support_face| {
                    face_edge_endpoints_fit_face_plane(
                        brep,
                        *support_face,
                        face_edge.support_half_edge,
                    )
                })
                .unwrap_or(face_edge.support_face)
        } else {
            face_edge.support_face
        };
        if let Some(half_edges) = by_support_face.get_mut(&group_key)
            && !half_edge_is_degenerate(brep, face_edge.support_half_edge)
        {
            if group_key != face_edge.support_face {
                splice_side_out_of_donor_loop(brep, face_edge.support_half_edge);
                if let Some(face) = brep.faces.get_mut(face_edge.support_face) {
                    face.role = TopologyRole::ConstructionResidual;
                    face.surface = None;
                }
            }
            half_edges.push(face_edge.support_half_edge);
            support_faces_with_existing_boundary.insert(group_key);
        }
    }

    let mut support_groups = by_support_face.into_iter().collect::<Vec<_>>();
    support_groups.sort_by_key(|(support_face, _)| {
        let generated_count = generated_trim_count_by_support_face
            .get(support_face)
            .copied()
            .unwrap_or_default();
        (
            std::cmp::Reverse(generated_count),
            support_faces_with_existing_boundary.contains(support_face),
            *support_face,
        )
    });

    for (support_face, mut half_edges) in support_groups {
        trace_topology_counts(brep, "support-face-before");
        half_edges.retain(|half_edge| !half_edge_is_degenerate(brep, *half_edge));
        if let Some(spliced) = assemble_extended_facet_cycle(brep, support_face, &half_edges)? {
            let loop_handle = brep
                .faces
                .get(support_face)
                .map(|face| face.outer_loop)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition cannot replace stale support face                          {support_face:?}"
                    ),
                })?;
            let displaced_old: Vec<HalfEdgeHandle> = collect_loop_half_edges(brep, loop_handle)
                .unwrap_or_default()
                .into_iter()
                .filter(|side| !spliced.contains(side))
                .collect();
            wire_existing_loop(brep, loop_handle, support_face, &spliced)?;
            for side in &displaced_old {
                let Some(data) = brep.half_edges.get(*side) else {
                    continue;
                };
                let side_start = data.origin;
                let Some(side_end) = brep.half_edges.get(data.twin).map(|twin| twin.origin) else {
                    continue;
                };
                let mut holding = vec![*side];
                if side_end != side_start {
                    let (seam_edge, seam_side, _) =
                        insert_unwired_edge_pair(brep, side_end, side_start);
                    brep.edges[seam_edge].role = TopologyRole::ConstructionResidual;
                    holding.push(seam_side);
                }
                let remainder_face = insert_face_with_outer_loop(brep, support_face, &holding)?;
                brep.faces[remainder_face].role = TopologyRole::ConstructionResidual;
                brep.faces[remainder_face].surface = None;
            }
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/support-loop-replace] support {support_face:?}                      explicit-assembled ({} sides; {} displaced retired)",
                    spliced.len(),
                    displaced_old.len()
                );
            }
            trace_topology_counts(brep, "support-face-after");
            continue;
        }
        let support_split_plan = materialize_support_split_cycles(brep, support_face, &half_edges)?;
        authored_boundary_spans.extend(support_split_plan.authored_boundary_spans.iter().copied());
        if !support_split_plan.is_complete() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot materialize support ownership on \
                     {support_face:?}: generated trim sides {:?} do not form closed support \
                     cycle(s) yet. Interior support subdivision must add the missing support \
                     boundary spans before these cylinder trims can be exposed. {}; {}",
                    half_edges,
                    support_cycle_open_chain_diagnostic(brep, &half_edges),
                    support_split_plan.diagnostic()
                ),
            });
        }
        if support_split_plan.cycles.len() != 1 {
            let counts = topology_count_summary(brep);
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot materialize support ownership on \
                     {support_face:?}: got {} closed support cycles. Multi-cycle support \
                     subdivision requires an Euler-balanced residual face/ring operator; direct \
                     inner-loop wiring or duplicate-face insertion would violate shell topology. \
                     Current topology deficit before this support face is {}; one face-loop per \
                     closed support cycle would add only {}.",
                    support_split_plan.cycles.len(),
                    counts.deficit_to_genus_zero(),
                    support_split_plan.cycles.len()
                ),
            });
        }
        let cycle = &support_split_plan.cycles[0];
        validate_half_edge_cycle(brep, cycle, "materialized support trim cycle", support_face)?;
        let loop_handle = brep
            .faces
            .get(support_face)
            .map(|face| face.outer_loop)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot replace stale support face \
                     {support_face:?}"
                ),
            })?;
        wire_existing_loop(brep, loop_handle, support_face, cycle)?;
        trace_topology_counts(brep, "support-face-after");
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/support-loop-replace] support {support_face:?} -> \
                 rewired outer loop {:?} using {:?}",
                loop_handle, cycle
            );
        }
    }
    retire_consumed_support_source_loops(brep, &source_loop_transfers)?;
    normalize_spilled_fillet_corners(brep, trim_sides, &source_loop_transfers)?;
    materialize_support_boundary_residual_cycles(brep, &authored_boundary_spans)?;
    Ok(())
}

fn retire_consumed_support_source_loops(
    brep: &mut BRep,
    transfers: &[SupportSourceLoopTransfer],
) -> Result<(), ApexBuildError> {
    let loops = transfers
        .iter()
        .map(|transfer| transfer.loop_)
        .collect::<std::collections::BTreeSet<_>>();
    for loop_handle in loops {
        let Some(loop_data) = brep.loops.get(loop_handle).cloned() else {
            continue;
        };
        if brep
            .half_edges
            .iter()
            .any(|(_, half_edge)| half_edge.loop_ == loop_handle)
        {
            continue;
        }
        let Some(face_data) = brep.faces.get(loop_data.face) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot retire consumed support source loop \
                     {loop_handle:?}: owner face {:?} is stale",
                    loop_data.face
                ),
            });
        };
        if face_data.outer_loop != loop_handle || !face_data.inner_loops.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot retire consumed support source loop \
                     {loop_handle:?}: owner face {:?} has non-simple loop ownership",
                    loop_data.face
                ),
            });
        }
        let shell = face_data.shell;
        brep.loops.remove(loop_handle);
        brep.faces.remove(loop_data.face);
        if brep
            .shells
            .get(shell)
            .is_some_and(|shell_data| shell_data.face == loop_data.face)
            && let Some(replacement_face) = brep
                .faces
                .iter()
                .find_map(|(face, data)| (data.shell == shell).then_some(face))
        {
            brep.shells[shell].face = replacement_face;
        }
    }
    Ok(())
}

fn trace_topology_counts(brep: &BRep, label: &'static str) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    eprintln!("{}", topology_count_summary(brep).trace_label(label));
}
