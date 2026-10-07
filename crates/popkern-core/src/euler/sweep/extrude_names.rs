//! Stable names for everything an extrude creates (split from
//! `extrude.rs`).

use crate::history::TrackingContext;
use crate::naming::{EntityName, EntityRole, NameId, RoleCategory};
use crate::naming_criteria::{GeometricCriteria, capture_edge, capture_face, capture_vertex};
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Assign stable topological names to all entities produced by an extrude.
///
/// The extrude creates:
/// - N lateral (side) faces, one per profile edge
/// - 1 cap (top) face
/// - The bottom face (the pre-existing face opposite the cap) is marked Modified
/// - N vertical edges (connecting profile vertices to their extruded copies)
/// - N top edges (connecting adjacent extruded vertices)
/// - N new vertices (extruded copies of the profile vertices)
///
/// Note: the input `face` handle becomes one of the side faces after extrude
/// because each `mef` carves off a portion of the original face. The bottom
/// face is the OTHER pre-existing face (typically the mvfs face).
#[allow(clippy::too_many_arguments)]
pub(super) fn assign_extrude_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    bottom_face: Option<FaceHandle>,
    input_face_name: Option<NameId>,
    boundary_edge_names: &[Option<NameId>],
    boundary_vertex_names: &[Option<NameId>],
    cap_face: FaceHandle,
    side_faces: &[FaceHandle],
    new_edges: &[EdgeHandle],
    new_vertices: &[VertexHandle],
) {
    let n = side_faces.len();
    let input_ancestors: Vec<NameId> = input_face_name.into_iter().collect();

    // --- Bottom face (modified) ---
    if let Some(bottom_fh) = bottom_face {
        let bottom_name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Modified,
                index: 0,
            },
            ancestors: input_ancestors.clone(),
        };
        let criteria = capture_face(brep, bottom_fh).map(GeometricCriteria::Face);
        let bottom_name_id = brep.name_pool.intern(bottom_name.clone(), criteria);
        brep.faces[bottom_fh].name = Some(bottom_name_id);
        if let Some(old_name) = input_face_name.and_then(|id| brep.name_pool.get(id).cloned()) {
            ctx.record_modified(old_name, bottom_name);
        } else {
            ctx.record_generated(bottom_name);
        }
    }

    // --- Cap (top) face ---
    let cap_name = EntityName {
        operation: ctx.operation,
        role: EntityRole {
            category: RoleCategory::Cap,
            index: 0,
        },
        ancestors: input_ancestors,
    };
    let criteria = capture_face(brep, cap_face).map(GeometricCriteria::Face);
    let cap_name_id = brep.name_pool.intern(cap_name.clone(), criteria);
    brep.faces[cap_face].name = Some(cap_name_id);
    ctx.record_generated(cap_name);

    // --- Lateral (side) faces ---
    for (i, &side_fh) in side_faces.iter().enumerate() {
        let ancestors: Vec<NameId> = boundary_edge_names
            .get(i)
            .copied()
            .flatten()
            .into_iter()
            .collect();
        let lateral_name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Lateral,
                index: i as u32,
            },
            ancestors,
        };
        let criteria = capture_face(brep, side_fh).map(GeometricCriteria::Face);
        let name_id = brep.name_pool.intern(lateral_name.clone(), criteria);
        brep.faces[side_fh].name = Some(name_id);
        ctx.record_generated(lateral_name);
    }

    // --- New vertices (extruded copies) ---
    // new_vertices[i] is the extruded copy of profile vertex i.
    for (i, &vh) in new_vertices.iter().enumerate() {
        let ancestors: Vec<NameId> = boundary_vertex_names
            .get(i)
            .copied()
            .flatten()
            .into_iter()
            .collect();
        let vertex_name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors,
        };
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let name_id = brep.name_pool.intern(vertex_name.clone(), criteria);
        brep.vertices[vh].name = Some(name_id);
        ctx.record_generated(vertex_name);
    }

    // --- New edges ---
    // new_edges layout: first N are vertical edges (mev), next N are top edges (mef).
    // Vertical edges: new_edges[0..n] — edge connecting profile vertex i to its
    // extruded copy.
    for i in 0..n {
        if let Some(&eh) = new_edges.get(i) {
            let ancestors: Vec<NameId> = boundary_vertex_names
                .get(i)
                .copied()
                .flatten()
                .into_iter()
                .collect();
            let edge_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Lateral,
                    index: i as u32,
                },
                ancestors,
            };
            let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
            let name_id = brep.name_pool.intern(edge_name.clone(), criteria);
            brep.edges[eh].name = Some(name_id);
            ctx.record_generated(edge_name);
        }
    }

    // Top edges: new_edges[n..2n] — edge connecting adjacent extruded vertices.
    for i in 0..n {
        if let Some(&eh) = new_edges.get(n + i) {
            let ancestors: Vec<NameId> = boundary_edge_names
                .get(i)
                .copied()
                .flatten()
                .into_iter()
                .collect();
            let edge_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Cap,
                    index: i as u32,
                },
                ancestors,
            };
            let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
            let name_id = brep.name_pool.intern(edge_name.clone(), criteria);
            brep.edges[eh].name = Some(name_id);
            ctx.record_generated(edge_name);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// Build a square face and return the face handle.
    fn build_square_face(brep: &mut BRep) -> crate::topology::FaceHandle {
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        mef.face
    }

    /// Build a triangle face.
    fn build_triangle_face(brep: &mut BRep) -> crate::topology::FaceHandle {
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(0.5, 1.0, 0.0))
            .unwrap();
        let mef = brep.mef(mev2.half_edge_twin, mev1.half_edge).unwrap();
        mef.face
    }

    #[test]
    fn extrude_square_creates_cube_topology() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_square_face(&mut brep);

        let _result = brep
            .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
            .unwrap();

        // Cube: V=8, E=12, F=6
        assert_eq!(brep.vertices.len(), 8, "Expected 8 vertices");
        assert_eq!(brep.edges.len(), 12, "Expected 12 edges");
        assert_eq!(brep.faces.len(), 6, "Expected 6 faces");
    }

    #[test]
    fn extrude_triangle_creates_prism_topology() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_triangle_face(&mut brep);

        let _result = brep
            .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
            .unwrap();

        // Triangular prism: V=6, E=9, F=5
        assert_eq!(brep.vertices.len(), 6, "Expected 6 vertices");
        assert_eq!(brep.edges.len(), 9, "Expected 9 edges");
        assert_eq!(brep.faces.len(), 5, "Expected 5 faces");
    }

    #[test]
    fn extrude_invalid_face_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.extrude(Default::default(), Vec3::new(0.0, 0.0, 1.0), 1.0, None);
        assert!(result.is_err());
    }

    #[test]
    fn extrude_face_without_edges_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // mvfs face has no edges
        let result = brep.extrude(mvfs.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None);
        assert!(result.is_err());
    }

    #[test]
    fn extruded_vertices_have_correct_positions() {
        let mut brep = BRep::new(Tolerance::default());
        let face = build_square_face(&mut brep);

        let result = brep
            .extrude(face, Vec3::new(0.0, 0.0, 2.0), 1.0, None)
            .unwrap();

        // All new vertices should have z = 2.0
        for &vh in &result.new_vertices {
            let v = &brep.vertices[vh];
            assert!(
                (v.position.z - 2.0).abs() < 1e-10,
                "New vertex z should be 2.0"
            );
        }
    }
}
