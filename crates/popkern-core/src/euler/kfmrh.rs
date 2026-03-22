//! kfmrh — Kill Face, Make Ring-Hole
//! mfkrh — Make Face, Kill Ring-Hole

use crate::error::CoreError;
use crate::topology::{BRep, Face, FaceHandle, LoopHandle, Shell, ShellHandle};

/// Result of the mfkrh (Make Face, Kill Ring-Hole) operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MfkrhResult {
    /// The newly created face.
    pub face: FaceHandle,
    /// The newly created shell.
    pub shell: ShellHandle,
}

impl BRep {
    /// Kill Face, Make Ring-Hole — merges two shells by removing a face from one
    /// shell and making its outer loop an inner loop of a face in another shell.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if either face handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the faces are in the same shell,
    ///   or the face to kill has inner loops.
    pub fn kfmrh(&mut self, face: FaceHandle, target_face: FaceHandle) -> Result<(), CoreError> {
        // Validate handles.
        let face_data = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("face handle is invalid".into()))?;
        let target_data = self
            .faces
            .get(target_face)
            .ok_or_else(|| CoreError::InvalidHandle("target_face handle is invalid".into()))?;

        // Precondition: faces must be in different shells.
        if face_data.shell == target_data.shell {
            return Err(CoreError::PreconditionFailed(
                "face and target_face must be in different shells".into(),
            ));
        }

        // Precondition: the face to kill must have no inner loops.
        if !face_data.inner_loops.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "face to kill must have no inner loops".into(),
            ));
        }

        let loop_to_move = face_data.outer_loop;
        let killed_shell = face_data.shell;
        let killed_solid = self.shells[killed_shell].solid;

        // Move the loop to become an inner loop of target_face.
        self.loops[loop_to_move].face = target_face;
        self.faces[target_face].inner_loops.push(loop_to_move);

        // Update all faces that belonged to the killed shell to point to the target shell.
        let target_shell = self.faces[target_face].shell;
        let faces_to_update: Vec<FaceHandle> = self
            .faces
            .iter()
            .filter(|(fh, fd)| fd.shell == killed_shell && *fh != face)
            .map(|(fh, _)| fh)
            .collect();
        for fh in faces_to_update {
            self.faces[fh].shell = target_shell;
        }

        // Remove the killed face.
        self.faces.remove(face);

        // Remove the killed shell.
        self.shells.remove(killed_shell);

        // Remove the killed solid.
        self.solids.remove(killed_solid);

        Ok(())
    }

    /// Make Face, Kill Ring-Hole — promotes an inner loop of a face to a new
    /// face in a new shell. Inverse of kfmrh.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the loop handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the loop is not an inner loop of any face.
    pub fn mfkrh(&mut self, inner_loop: LoopHandle) -> Result<MfkrhResult, CoreError> {
        // Validate handle.
        let loop_data = self
            .loops
            .get(inner_loop)
            .ok_or_else(|| CoreError::InvalidHandle("inner_loop handle is invalid".into()))?;

        let owner_face = loop_data.face;

        // Find and remove this loop from the owner face's inner_loops.
        let face_data = self
            .faces
            .get(owner_face)
            .ok_or_else(|| CoreError::InvalidHandle("owner face handle is invalid".into()))?;

        let pos = face_data
            .inner_loops
            .iter()
            .position(|&l| l == inner_loop)
            .ok_or_else(|| {
                CoreError::PreconditionFailed("loop is not an inner loop of any face".into())
            })?;

        self.faces[owner_face].inner_loops.remove(pos);

        // Create a new face with this loop as its outer_loop.
        let new_face = self.faces.insert(Face {
            outer_loop: inner_loop,
            inner_loops: vec![],
            shell: Default::default(), // patched below
            surface: None,
        });

        // Update the loop to point to the new face.
        self.loops[inner_loop].face = new_face;

        // Create a new shell containing this face.
        let new_shell = self.shells.insert(Shell {
            face: new_face,
            solid: Default::default(), // patched below
        });

        // Create a new solid for the new shell.
        let new_solid = self
            .solids
            .insert(crate::topology::Solid { shell: new_shell });

        // Patch references.
        self.faces[new_face].shell = new_shell;
        self.shells[new_shell].solid = new_solid;

        Ok(MfkrhResult {
            face: new_face,
            shell: new_shell,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn kfmrh_merges_shells() {
        let mut brep = BRep::new(Tolerance::default());

        // Build first solid: a square face
        let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let _mef1 = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

        // Build second solid: another square face (inner)
        let mvfs2 = brep.mvfs(Point3::new(0.2, 0.2, 0.0)).unwrap();
        let mev4 = brep
            .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(0.8, 0.2, 0.0))
            .unwrap();
        let mev5 = brep
            .mev(mev4.half_edge_twin, Point3::new(0.8, 0.8, 0.0))
            .unwrap();
        let mev6 = brep
            .mev(mev5.half_edge_twin, Point3::new(0.2, 0.8, 0.0))
            .unwrap();
        let _mef2 = brep.mef(mev6.half_edge_twin, mev4.half_edge).unwrap();

        assert_eq!(brep.shells.len(), 2);
        assert_eq!(brep.solids.len(), 2);

        // kfmrh: kill a face from shell 2, make it an inner loop of a face in shell 1
        brep.kfmrh(mvfs2.face, mvfs1.face).unwrap();

        // After kfmrh: target face should have an inner loop
        let target = &brep.faces[mvfs1.face];
        assert!(
            !target.inner_loops.is_empty(),
            "Target face should have inner loop"
        );
    }

    #[test]
    fn kfmrh_invalid_face_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.kfmrh(Default::default(), mvfs.face);
        assert!(result.is_err());
        let result = brep.kfmrh(mvfs.face, Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn kfmrh_same_shell_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // Same shell, same face — should fail
        let result = brep.kfmrh(mvfs1.face, mvfs1.face);
        assert!(result.is_err());
    }

    #[test]
    fn kfmrh_face_with_inner_loops_fails() {
        let mut brep = BRep::new(Tolerance::default());
        // Build two solids, one with inner loops
        let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let _mef1 = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

        // Second solid for inner
        let mvfs2 = brep.mvfs(Point3::new(0.2, 0.2, 0.0)).unwrap();
        let mev4 = brep
            .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(0.8, 0.2, 0.0))
            .unwrap();
        let mev5 = brep
            .mev(mev4.half_edge_twin, Point3::new(0.8, 0.8, 0.0))
            .unwrap();
        let mev6 = brep
            .mev(mev5.half_edge_twin, Point3::new(0.2, 0.8, 0.0))
            .unwrap();
        let _mef2 = brep.mef(mev6.half_edge_twin, mev4.half_edge).unwrap();

        // kfmrh to create inner loop in face1
        brep.kfmrh(mvfs2.face, mvfs1.face).unwrap();

        // Now try kfmrh on mvfs1.face which has inner loops — third solid needed
        let mvfs3 = brep.mvfs(Point3::new(5.0, 5.0, 0.0)).unwrap();
        // Try to kill face1 (which has inner loops) into mvfs3's face
        let result = brep.kfmrh(mvfs1.face, mvfs3.face);
        assert!(result.is_err());
    }

    #[test]
    fn mfkrh_invalid_loop_handle_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let result = brep.mfkrh(Default::default());
        assert!(result.is_err());
    }

    #[test]
    fn mfkrh_outer_loop_fails() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // The outer loop is not an inner loop — should fail
        let result = brep.mfkrh(mvfs.loop_);
        assert!(result.is_err());
    }

    #[test]
    fn mfkrh_separates_shells() {
        let mut brep = BRep::new(Tolerance::default());

        // Same setup, then kfmrh, then mfkrh to undo
        let mvfs1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs1.vertex, mvfs1.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let _mef1 = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();

        let mvfs2 = brep.mvfs(Point3::new(0.2, 0.2, 0.0)).unwrap();
        let mev4 = brep
            .mev_at_vertex(mvfs2.vertex, mvfs2.loop_, Point3::new(0.8, 0.2, 0.0))
            .unwrap();
        let mev5 = brep
            .mev(mev4.half_edge_twin, Point3::new(0.8, 0.8, 0.0))
            .unwrap();
        let mev6 = brep
            .mev(mev5.half_edge_twin, Point3::new(0.2, 0.8, 0.0))
            .unwrap();
        let _mef2 = brep.mef(mev6.half_edge_twin, mev4.half_edge).unwrap();

        brep.kfmrh(mvfs2.face, mvfs1.face).unwrap();

        // mfkrh: separate the inner loop back out
        let inner_loop = brep.faces[mvfs1.face].inner_loops[0];
        let result = brep.mfkrh(inner_loop).unwrap();

        assert!(brep.faces.contains_key(result.face));
        assert!(brep.shells.contains_key(result.shell));
        assert!(brep.faces[mvfs1.face].inner_loops.is_empty());
    }
}
