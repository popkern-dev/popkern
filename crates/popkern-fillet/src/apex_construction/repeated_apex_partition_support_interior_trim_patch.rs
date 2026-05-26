//! Interior support trim patch authoring for repeated-apex cylinder rewrites.

use super::*;

pub(super) struct InteriorSupportTrimPatchRequest {
    pub(super) support_face: FaceHandle,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) curve: CurveOnEdge,
    pub(super) label: &'static str,
    pub(super) fillet_index: usize,
}

pub(super) fn author_interior_support_trim_patch_for_cylinder_trim(
    brep: &mut BRep,
    request: InteriorSupportTrimPatchRequest,
) -> Result<Option<MaterializedGeneratorTrimSide>, ApexBuildError> {
    let InteriorSupportTrimPatchRequest {
        support_face,
        start,
        end,
        curve,
        label,
        fillet_index,
    } = request;
    let Some(connector_curve) = support_trim_patch_connector_curve(brep, start, end) else {
        return Ok(None);
    };
    let rollback = brep.clone();
    let result = (|| {
        let (edge, cylinder_half_edge, support_half_edge) =
            insert_unwired_edge_pair(brep, start, end);
        brep.edges[edge].curve = Some(curve);
        let (connector_edge, connector_half_edge, _connector_twin) =
            insert_unwired_edge_pair(brep, start, end);
        brep.edges[connector_edge].curve = Some(connector_curve);
        brep.edges[connector_edge].role = TopologyRole::ConstructionResidual;
        let support_patch_face = insert_face_with_outer_loop(
            brep,
            support_face,
            &[support_half_edge, connector_half_edge],
        )?;
        brep.faces[support_patch_face].role = TopologyRole::ConstructionResidual;
        brep.faces[support_patch_face].surface = brep
            .faces
            .get(support_face)
            .and_then(|face| face.surface.clone());
        for half_edge in [support_half_edge, connector_half_edge] {
            let edge = brep.half_edges[half_edge].edge;
            brep.edges[edge].role = TopologyRole::ConstructionResidual;
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] authored {label} \
                 {start:?}->{end:?} as interior support trim patch {support_patch_face:?}; \
                 cylinder side {cylinder_half_edge:?}, support side {support_half_edge:?}"
            );
        }
        Ok(MaterializedGeneratorTrimSide {
            half_edge: cylinder_half_edge,
            support_face: support_patch_face,
            support_side_owned: true,
            support_source: None,
            authored_boundary_spans: vec![connector_half_edge],
        })
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result.map(Some)
}

fn support_trim_patch_connector_curve(
    brep: &BRep,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<CurveOnEdge> {
    let start = brep.vertices.get(start)?.position;
    let end = brep.vertices.get(end)?.position;
    CurveOnEdge::line_segment(start, end, brep.tolerance.linear())
}
