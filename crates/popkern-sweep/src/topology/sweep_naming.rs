//! Stable names for sweep, loft, and revolve results: every face,
//! edge, and vertex is `Generated` under the operation (sweeps create
//! entirely new geometry rather than modifying existing entities).

use popkern_core::history::TrackingContext;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Assign stable [`EntityName`]s to all entities in a sweep BRep.
///
/// All faces, edges, and vertices receive `Generated` names since sweeps
/// create entirely new geometry rather than modifying existing entities.
pub(crate) fn assign_sweep_names(brep: &mut BRep, ctx: &mut TrackingContext) {
    use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
    use popkern_core::naming_criteria::GeometricCriteria;
    use popkern_core::naming_criteria::{capture_edge, capture_face, capture_vertex};

    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    for (i, &fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[fh].name = Some(name_id);
        ctx.record_generated(name);
    }

    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for (i, &eh) in edge_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.edges[eh].name = Some(name_id);
        ctx.record_generated(name);
    }

    let vertex_handles: Vec<VertexHandle> = brep.vertices.keys().collect();
    for (i, &vh) in vertex_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.vertices[vh].name = Some(name_id);
        ctx.record_generated(name);
    }
}
