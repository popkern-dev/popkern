use super::*;

#[derive(Clone)]
pub(super) struct SupportOwnedTrimSide {
    pub(super) cylinder_face: FaceHandle,
    pub(super) support_face: FaceHandle,
    pub(super) fillet_index: usize,
    pub(super) half_edge: HalfEdgeHandle,
    pub(super) support_side_owned: bool,
    pub(super) support_source: Option<SupportOwnedTrimSource>,
    pub(super) requires_support_ownership: bool,
}

#[derive(Clone, Copy)]
pub(super) struct SupportOwnedTrimSource {
    pub(super) loop_: LoopHandle,
    pub(super) prev: HalfEdgeHandle,
    pub(super) next: HalfEdgeHandle,
}

#[derive(Clone)]
pub(super) struct PreMaterializedTrimSide {
    pub(super) cylinder_face: FaceHandle,
    pub(super) fillet_index: usize,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) trim_sides: Vec<SupportOwnedTrimSide>,
}

pub(super) fn pre_materialized_trim_sides(
    pre_materialized: &[PreMaterializedTrimSide],
    cylinder_face: FaceHandle,
    fillet_index: usize,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<SupportOwnedTrimSide>> {
    pre_materialized
        .iter()
        .find(|entry| {
            entry.cylinder_face == cylinder_face
                && entry.fillet_index == fillet_index
                && entry.start == start
                && entry.end == end
        })
        .map(|entry| entry.trim_sides.clone())
}

pub(super) fn pre_materialize_far_side_trim_sides(
    brep: &mut BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
    include_selected_nonmaterialized: bool,
) -> Result<Vec<PreMaterializedTrimSide>, ApexBuildError> {
    let mut materialized = Vec::new();
    for (rebuild_plan, topology_plan) in rebuild_plans.iter().zip(topology_plans) {
        let Some(fragment) =
            rebuild_plan.fragment_for_strategy(topology_plan.loop_candidate.far_side_strategy)
        else {
            continue;
        };
        if !fragment.materialized_anchor_replacement && !include_selected_nonmaterialized {
            continue;
        }
        let allowed_support_faces =
            target_cylinder_support_faces(brep, rebuild_plan.cylinder_face, cap_face, face_edges);
        for segment in &topology_plan.loop_candidate.segments {
            if !segment.label.contains("far-side") {
                continue;
            }
            let trim_sides = materialized_far_side_trim_sides(
                brep,
                rebuild_plan.cylinder_face,
                segment.start,
                segment.end,
                segment.label,
                rebuild_plan.fillet_index,
                &allowed_support_faces,
            )?;
            materialized.push(PreMaterializedTrimSide {
                cylinder_face: rebuild_plan.cylinder_face,
                fillet_index: rebuild_plan.fillet_index,
                start: segment.start,
                end: segment.end,
                trim_sides,
            });
        }
    }
    Ok(materialized)
}
