//! Extrude sweep — a higher-level operation built from Euler operators (mev + mef).

use crate::error::CoreError;
use crate::geometry::Vec3;
use crate::history::TrackingContext;
use crate::naming::NameId;
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
    /// Each base-loop edge paired with the top (cap-rim) edge that
    /// `mef` created above it, in loop order. Recorded at Euler time
    /// so geometry authoring never re-derives the pairing from vertex
    /// pairs — ambiguous for a bigon lamina (a split circle's two
    /// half-arcs, a line-plus-arc "D"), where every top edge received
    /// the SAME base curve and one wall tessellated to nothing.
    pub top_edge_of_base: Vec<(EdgeHandle, EdgeHandle)>,
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
    /// - [`CoreError::PreconditionFailed`] if the face has no edges in its outer loop,
    ///   or when the finished topology fails the Euler postconditions (the
    ///   body is then unchanged).
    ///
    /// The extrude is one transaction: its n `mev` and n `mef` steps run
    /// on a working copy, which must show exactly their combined Euler
    /// delta and pass the structural checks (half-edge pairing, cycles,
    /// vertex anchors, manifold edges) before it replaces the body.
    /// Verifying once at this boundary instead of after every step keeps
    /// the extrude O(n) in checks (PK-709: per-step whole-model checks made
    /// a 1000-edge profile take 52 s), and a failure part-way leaves the
    /// body untouched instead of half extruded.
    pub fn extrude(
        &mut self,
        face: FaceHandle,
        direction: Vec3,
        distance: f64,
        tracking: Option<&mut TrackingContext>,
    ) -> Result<ExtrudeResult, CoreError> {
        use crate::euler::contract::{EulerDelta, EulerTopologyCounts};
        let before = EulerTopologyCounts::from_brep(self);
        let mut work = self.clone();
        let result = work.extrude_steps(face, direction, distance, tracking)?;
        let n = result.side_faces.len() as isize;
        EulerDelta::verify(
            "extrude",
            before,
            EulerTopologyCounts::from_brep(&work),
            EulerDelta::mev().times(n).then(EulerDelta::mef().times(n)),
        )?;
        work.validate_half_edge_pairing().map_err(CoreError::from)?;
        work.validate_half_edge_cycles().map_err(CoreError::from)?;
        work.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        work.validate_manifold_edges().map_err(CoreError::from)?;
        *self = work;
        Ok(result)
    }

    /// The extrude's Euler steps and geometry, unverified (see
    /// [`BRep::extrude`]).
    fn extrude_steps(
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

            let mev_result = self.apply_mev(loop_half_edges[i], new_pos)?;

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
        let mut top_edge_of_base: Vec<(EdgeHandle, EdgeHandle)> = Vec::with_capacity(n);

        for i in 0..n {
            let (he_a, he_b) = if i < n - 1 {
                (in_half_edges[i + 1], in_half_edges[i])
            } else {
                (top0_he, in_half_edges[n - 1])
            };

            let mef_result = self.apply_mef(he_a, he_b)?;
            new_edges.push(mef_result.edge);
            top_edge_of_base.push((self.half_edges[loop_half_edges[i]].edge, mef_result.edge));

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

            super::extrude_names::assign_extrude_names(
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
            top_edge_of_base,
        };
        // PK-316: author the geometry the sweep created — the Euler
        // topology alone leaves mef-inherited base surfaces on every
        // new face and no curves on the new edges (mechanism and
        // measured failure in `extrude_geometry`).
        super::extrude_geometry::author_extrude_geometry(self, &result, face, displacement)?;
        Ok(result)
    }
}
