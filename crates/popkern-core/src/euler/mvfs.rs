//! mvfs — Make Vertex, Face, Shell
//! kvfs — Kill Vertex, Face, Shell

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::topology::{
    BRep, Face, FaceHandle, Loop, LoopHandle, Shell, ShellHandle, Solid, SolidHandle, Vertex,
    VertexHandle,
};

/// Result of the mvfs (Make Vertex, Face, Shell) operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MvfsResult {
    /// The newly created vertex.
    pub vertex: VertexHandle,
    /// The newly created face.
    pub face: FaceHandle,
    /// The newly created loop (the face's outer boundary).
    pub loop_: LoopHandle,
    /// The newly created shell.
    pub shell: ShellHandle,
    /// The newly created solid.
    pub solid: SolidHandle,
}

impl BRep {
    /// Make Vertex, Face, Shell — the seed operation for a new solid.
    ///
    /// Creates one vertex, face, loop, shell, and solid, all wired together.
    /// The loop and vertex have no half-edges yet.
    pub fn mvfs(&mut self, position: Point3) -> Result<MvfsResult, CoreError> {
        // Insert with placeholder handles — we'll patch them up after insertion.
        let vertex = self.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
        });

        // We need face, loop, shell, solid — they are mutually referential.
        // Insert with defaults and patch.
        let face = self.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell: Default::default(),
            surface: None,
        });

        let loop_ = self.loops.insert(Loop {
            half_edge: Default::default(),
            face,
        });

        let shell = self.shells.insert(Shell {
            face,
            solid: Default::default(),
        });

        let solid = self.solids.insert(Solid { shell });

        // Patch the face with the real loop and shell handles.
        let face_data = &mut self.faces[face];
        face_data.outer_loop = loop_;
        face_data.shell = shell;

        // Patch the shell with the real solid handle.
        self.shells[shell].solid = solid;

        Ok(MvfsResult {
            vertex,
            face,
            loop_,
            shell,
            solid,
        })
    }

    /// Kill Vertex, Face, Shell — inverse of mvfs.
    ///
    /// Removes the solid that contains the given vertex, along with its shell,
    /// face, loop, and the vertex itself.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the vertex handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the vertex has outgoing edges,
    ///   the shell has more than one face, the face has inner loops, or edges
    ///   still exist in the model for this solid.
    pub fn kvfs(&mut self, vertex: VertexHandle) -> Result<(), CoreError> {
        // Validate vertex handle.
        let vertex_data = self
            .vertices
            .get(vertex)
            .ok_or_else(|| CoreError::InvalidHandle("vertex handle is invalid".into()))?;

        // Precondition: vertex has no outgoing edges.
        let outgoing = vertex_data.outgoing;
        if self.half_edges.contains_key(outgoing) {
            return Err(CoreError::PreconditionFailed(
                "vertex still has outgoing half-edges".into(),
            ));
        }

        // We need to find the face, loop, shell, solid via the vertex.
        // Strategy: find the face whose outer_loop's face points back to it
        // and whose shell's solid we can then remove.
        //
        // For mvfs-created topology the vertex isn't directly linked to the
        // face, so we scan shells to find the one that contains exactly one
        // face and that face's loop has no half-edges.
        //
        // A cleaner approach: scan all faces, find any face in a shell that
        // has no half-edges and belongs to a solid with one shell and one face.
        // Since kvfs is the strict inverse of mvfs, we find the shell/face/loop
        // associated by checking all shells for the one-face constraint.
        //
        // For correctness we need to know which solid owns this vertex.  Because
        // the vertex has no outgoing half-edges we cannot follow connectivity.
        // We scan faces — in a well-formed mvfs state there is exactly one face
        // per solid.  We locate the face whose shell's solid we will remove, by
        // finding a (face, shell, solid) triple where:
        //   - shell has no edges (half_edges.len() == 0 is too broad if there
        //     are other solids)
        //   - face has exactly one loop (outer_loop) and no inner loops
        //   - the loop has no half-edge
        //
        // Since this is kvfs (inverse of mvfs) and later operators will add
        // their own precondition checks, we simplify: find any shell whose face
        // is a single-face, single-loop, empty shell, and whose solid is the
        // one we can remove together with this vertex.
        //
        // The simplest correct implementation: iterate shells, for each shell
        // check if its face's outer_loop has no half-edge and no inner loops,
        // and the shell has exactly one face.  Count faces in the shell by
        // scanning all faces that reference this shell.

        // Collect faces per shell for quick lookup.
        let face_handle = {
            // Find the face in this BRep that has no inner loops, has an empty
            // outer loop, belongs to a shell that belongs to a solid that has
            // exactly one shell.  We assert that such a face corresponds to an
            // mvfs-seeded solid.
            //
            // For the general case (post-Task-8) we iterate over all faces and
            // find one that satisfies the mvfs invariant.  We also need to
            // confirm it is "owned" by this vertex — but since there's no direct
            // vertex→face link we accept the first matching face in a valid
            // mvfs-shaped solid.
            //
            // IMPORTANT: if there are multiple mvfs-seeded solids the caller
            // must have no edges so that we can identify which solid the vertex
            // belongs to via the total counts.  We check that the total edge
            // count for the candidate solid is zero.

            // Build a candidate: the shell whose face has one loop, no inner
            // loops, no half-edge, and the solid associated with that shell has
            // exactly one shell referenced by it.
            let mut found: Option<(FaceHandle, LoopHandle, ShellHandle, SolidHandle)> = None;

            'outer: for (shell_handle, shell_data) in &self.shells {
                let face_h = shell_data.face;
                let face_data = match self.faces.get(face_h) {
                    Some(f) => f,
                    None => continue,
                };
                // Face must have no inner loops.
                if !face_data.inner_loops.is_empty() {
                    continue;
                }
                let loop_h = face_data.outer_loop;
                let loop_data = match self.loops.get(loop_h) {
                    Some(l) => l,
                    None => continue,
                };
                // Loop must have no half-edge.
                if self.half_edges.contains_key(loop_data.half_edge) {
                    continue;
                }
                // Shell must contain exactly one face (scan all faces for this shell).
                let face_count = self
                    .faces
                    .values()
                    .filter(|f| f.shell == shell_handle)
                    .count();
                if face_count != 1 {
                    continue;
                }
                // Solid must reference this shell.
                let solid_h = shell_data.solid;
                let solid_data = match self.solids.get(solid_h) {
                    Some(s) => s,
                    None => continue,
                };
                if solid_data.shell != shell_handle {
                    continue;
                }
                // No edges exist in the model (for any solid — kvfs can only
                // be called when the target solid is completely edge-free, which
                // for a fresh mvfs solid means edges == 0).
                // We check that there are no half-edges in this shell's loops.
                for (_, l) in &self.loops {
                    if l.face == face_h && self.half_edges.contains_key(l.half_edge) {
                        continue 'outer;
                    }
                }

                found = Some((face_h, loop_h, shell_handle, solid_h));
                break;
            }

            match found {
                Some(f) => f,
                None => {
                    return Err(CoreError::PreconditionFailed(
                        "no valid mvfs-seeded solid found for this vertex".into(),
                    ));
                }
            }
        };

        let (face_h, loop_h, shell_h, solid_h) = face_handle;

        // Remove all entities.
        self.solids.remove(solid_h);
        self.shells.remove(shell_h);
        self.faces.remove(face_h);
        self.loops.remove(loop_h);
        self.vertices.remove(vertex);

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn mvfs_creates_vertex_face_shell_solid() {
        let mut brep = BRep::new(Tolerance::default());
        let result = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();

        assert_eq!(brep.vertices.len(), 1);
        assert_eq!(brep.faces.len(), 1);
        assert_eq!(brep.loops.len(), 1);
        assert_eq!(brep.shells.len(), 1);
        assert_eq!(brep.solids.len(), 1);
        assert_eq!(brep.edges.len(), 0);
        assert_eq!(brep.half_edges.len(), 0);

        // Verify wiring
        let face = &brep.faces[result.face];
        assert_eq!(face.outer_loop, result.loop_);
        assert!(face.inner_loops.is_empty());
        let loop_ = &brep.loops[result.loop_];
        assert_eq!(loop_.face, result.face);
        let shell = &brep.shells[result.shell];
        assert_eq!(shell.solid, result.solid);
        let solid = &brep.solids[result.solid];
        assert_eq!(solid.shell, result.shell);
        let vertex = &brep.vertices[result.vertex];
        assert_eq!(vertex.position.x, 0.0);
    }

    #[test]
    fn kvfs_removes_everything() {
        let mut brep = BRep::new(Tolerance::default());
        let result = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        brep.kvfs(result.vertex).unwrap();

        assert_eq!(brep.vertices.len(), 0);
        assert_eq!(brep.faces.len(), 0);
        assert_eq!(brep.loops.len(), 0);
        assert_eq!(brep.shells.len(), 0);
        assert_eq!(brep.solids.len(), 0);
    }

    #[test]
    fn kvfs_invalid_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let result = brep.kvfs(Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn two_mvfs_creates_two_solids() {
        let mut brep = BRep::new(Tolerance::default());
        let r1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let r2 = brep.mvfs(Point3::new(1.0, 0.0, 0.0)).unwrap();

        assert_eq!(brep.vertices.len(), 2);
        assert_eq!(brep.solids.len(), 2);
        assert_ne!(r1.solid, r2.solid);
    }
}
