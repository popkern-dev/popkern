//! Extrude sweep — a higher-level operation built from Euler operators (mev + mef).

use crate::error::CoreError;
use crate::geometry::Vec3;
use crate::history::TrackingContext;
use crate::naming::{EntityName, EntityRole, NameId, RoleCategory};
use crate::naming_criteria::{GeometricCriteria, capture_edge, capture_face, capture_vertex};
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Result of the extrude sweep operation.
#[derive(Debug, Clone)]
pub struct ExtrudeResult {
    /// The cap (top) face created by extrusion.
    pub cap_face: FaceHandle,
    /// The side faces created during extrusion.
    pub side_faces: Vec<FaceHandle>,
    /// All new edges created during extrusion.
    pub new_edges: Vec<EdgeHandle>,
    /// All new vertices created during extrusion.
    pub new_vertices: Vec<VertexHandle>,
}

impl BRep {
    /// Extrude a face along a direction vector by a given distance.
    ///
    /// Creates side faces and a cap face by extruding each vertex of the face's
    /// outer loop along `direction * distance`, using only Euler operators (mev + mef).
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the face handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the face has no edges in its outer loop.
    pub fn extrude(
        &mut self,
        face: FaceHandle,
        direction: Vec3,
        distance: f64,
        tracking: Option<&mut TrackingContext>,
    ) -> Result<ExtrudeResult, CoreError> {
        let face_data = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("face handle is invalid".into()))?;

        let outer_loop = face_data.outer_loop;
        let loop_data = self
            .loops
            .get(outer_loop)
            .ok_or_else(|| CoreError::InvalidHandle("outer loop handle is invalid".into()))?;

        let start_he = loop_data.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return Err(CoreError::PreconditionFailed(
                "face has no edges in its outer loop".into(),
            ));
        }

        // Collect pre-extrude state for naming.
        let pre_extrude_faces: Vec<FaceHandle> = self.faces.keys().collect();
        let input_face_name: Option<NameId> = self.faces.get(face).and_then(|f| f.name);
        let boundary_edge_names: Vec<Option<NameId>> = {
            let mut names = Vec::new();
            let mut cursor = start_he;
            loop {
                let he = &self.half_edges[cursor];
                names.push(self.edges[he.edge].name);
                cursor = he.next;
                if cursor == start_he {
                    break;
                }
            }
            names
        };
        let boundary_vertex_names: Vec<Option<NameId>> = {
            let mut names = Vec::new();
            let mut cursor = start_he;
            loop {
                let he = &self.half_edges[cursor];
                names.push(self.vertices[he.origin].name);
                cursor = he.next;
                if cursor == start_he {
                    break;
                }
            }
            names
        };

        // Step 1: Walk the outer loop to collect half-edges and vertices.
        let mut loop_half_edges: Vec<HalfEdgeHandle> = Vec::new();
        let mut loop_vertices: Vec<VertexHandle> = Vec::new();

        let mut cursor = start_he;
        loop {
            let he_data = &self.half_edges[cursor];
            loop_half_edges.push(cursor);
            loop_vertices.push(he_data.origin);
            cursor = he_data.next;
            if cursor == start_he {
                break;
            }
        }

        let n = loop_vertices.len();
        let displacement = direction * distance;

        let mut new_vertices = Vec::with_capacity(n);
        let mut new_edges = Vec::with_capacity(3 * n);
        let mut side_faces = Vec::with_capacity(n);

        // Step 2: For each vertex, create an extruded vertex using mev.
        // mev(he_i) creates a spur edge before he_i: out(vi->vi') and in(vi'->vi).
        // We store the "in" half-edge (vi' -> vi) for each vertex.
        let mut in_half_edges: Vec<HalfEdgeHandle> = Vec::with_capacity(n);

        for i in 0..n {
            let old_vertex = loop_vertices[i];
            let old_pos = self.vertices[old_vertex].position;
            let new_pos = old_pos + displacement;

            let mev_result = self.mev(loop_half_edges[i], new_pos)?;

            new_vertices.push(mev_result.vertex);
            new_edges.push(mev_result.edge);
            in_half_edges.push(mev_result.half_edge_twin);
        }

        // Step 3: Close side faces using mef.
        //
        // After all mev calls, the loop has spurs at each vertex. We close side
        // faces by connecting adjacent extruded vertices.
        //
        // For i = 0..N-2: mef(in[i+1], in[i]) creates side face i.
        //   - The old face (kept with new_he) becomes the side face.
        //   - in[i+1] remains in the new (active) face loop.
        //
        // For i = N-1 (final): mef(top0_he, in[N-1]) creates the last side face
        // and the cap. top0_he is the twin from the first mef (origin=v0', in
        // the active loop).

        let mut top0_he: HalfEdgeHandle = HalfEdgeHandle::default();
        let mut cap_face = FaceHandle::default();

        for i in 0..n {
            let (he_a, he_b) = if i < n - 1 {
                (in_half_edges[i + 1], in_half_edges[i])
            } else {
                (top0_he, in_half_edges[n - 1])
            };

            let mef_result = self.mef(he_a, he_b)?;
            new_edges.push(mef_result.edge);

            if i == 0 {
                // The twin from the first mef has origin v0' and lives in
                // the active (new) face loop. Save for the final closure.
                top0_he = mef_result.half_edge_twin;
            }

            // The old face (containing new_he) is a side face.
            let side_loop = self.half_edges[mef_result.half_edge].loop_;
            side_faces.push(self.loops[side_loop].face);

            if i == n - 1 {
                // The new face from the final mef is the cap.
                cap_face = mef_result.face;
            }
        }

        // Assign stable names when tracking is active.
        if let Some(ctx) = tracking {
            // Find the bottom face: it's the pre-existing face that is NOT
            // in side_faces and NOT the cap. The input `face` handle becomes
            // a side face after extrude.
            let side_and_cap: std::collections::HashSet<FaceHandle> = {
                let mut s: std::collections::HashSet<FaceHandle> =
                    side_faces.iter().copied().collect();
                s.insert(cap_face);
                s
            };
            let bottom_face = pre_extrude_faces
                .iter()
                .copied()
                .find(|fh| !side_and_cap.contains(fh));

            assign_extrude_names(
                self,
                ctx,
                bottom_face,
                input_face_name,
                &boundary_edge_names,
                &boundary_vertex_names,
                cap_face,
                &side_faces,
                &new_edges,
                &new_vertices,
            );
        }

        let result = ExtrudeResult {
            cap_face,
            side_faces,
            new_edges,
            new_vertices,
        };
        // PK-316: author the geometry the sweep created — the Euler
        // topology alone leaves mef-inherited base surfaces on every
        // new face and no curves on the new edges (mechanism and
        // measured failure in `extrude_geometry`).
        super::extrude_geometry::author_extrude_geometry(self, &result, face, displacement)?;
        Ok(result)
    }
}

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
fn assign_extrude_names(
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
