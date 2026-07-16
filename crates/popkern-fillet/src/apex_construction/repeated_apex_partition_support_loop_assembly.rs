//! Explicit six-side extended-facet cycle assembly for the repeated-apex
//! support-loop rewrite (PK-092 finale): [genA-twin, top, genB-twin,
//! span_b, base-contact, span_a], where corner spans are the twins of
//! the base loop's tangent-span sides.

use super::super::*;

/// Explicit extended-facet cycle: [genA-twin, top, genB-twin,
/// span(anchor->old-tangent, the base side's displaced twin),
/// kept base-contact, span] — every side identified semantically; no
/// chain search. Returns Ok(None) for groups that are not the
/// two-generator facet shape.
pub(super) fn assemble_extended_facet_cycle(
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

pub(super) fn trace_assembly_decline(support_face: FaceHandle, reason: &str) {
    if crate::direct_apex_trace::enabled() {
        eprintln!("[direct-apex-trace §14.9/assembly] {support_face:?} declined: {reason}");
    }
}
