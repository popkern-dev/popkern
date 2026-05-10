//! kfmrh — Kill Face, Make Ring-Hole
//! mfkrh — Make Face, Kill Ring-Hole

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts};
use crate::topology::{BRep, Face, FaceHandle, LoopHandle, Shell, ShellHandle, Solid};

#[derive(Debug, Clone)]
struct KfmrhApplied {
    loop_to_move: LoopHandle,
    killed_shell: ShellHandle,
    killed_solid: crate::topology::SolidHandle,
    target_shell: ShellHandle,
    reassigned_faces: Vec<FaceHandle>,
}

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
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let applied = candidate.apply_kfmrh(face, target_face)?;
        candidate.verify_kfmrh_postconditions(before, face, target_face, &applied)?;
        *self = candidate;
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
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let owner_face = candidate.loop_owner_face(inner_loop)?;
        let result = candidate.apply_mfkrh(inner_loop)?;
        candidate.verify_mfkrh_postconditions(before, inner_loop, owner_face, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_kfmrh(
        &mut self,
        face: FaceHandle,
        target_face: FaceHandle,
    ) -> Result<KfmrhApplied, CoreError> {
        let face_data = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("face handle is invalid".into()))?
            .clone();
        let target_data = self
            .faces
            .get(target_face)
            .ok_or_else(|| CoreError::InvalidHandle("target_face handle is invalid".into()))?
            .clone();

        if face_data.shell == target_data.shell {
            return Err(CoreError::PreconditionFailed(
                "face and target_face must be in different shells".into(),
            ));
        }

        if !face_data.inner_loops.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "face to kill must have no inner loops".into(),
            ));
        }

        let loop_to_move = face_data.outer_loop;
        let killed_shell = face_data.shell;
        let killed_shell_data = self
            .shells
            .get(killed_shell)
            .ok_or_else(|| CoreError::InvalidHandle("face shell handle is invalid".into()))?;
        let target_shell_data = self.shells.get(target_data.shell).ok_or_else(|| {
            CoreError::InvalidHandle("target face shell handle is invalid".into())
        })?;
        let killed_solid = killed_shell_data.solid;
        let target_solid = target_shell_data.solid;
        if killed_solid == target_solid {
            return Err(CoreError::PreconditionFailed(
                "face and target_face must belong to different solids".into(),
            ));
        }
        if !self.solids.contains_key(killed_solid) {
            return Err(CoreError::InvalidHandle(
                "face shell solid handle is invalid".into(),
            ));
        }
        if !self.solids.contains_key(target_solid) {
            return Err(CoreError::InvalidHandle(
                "target face shell solid handle is invalid".into(),
            ));
        }
        let loop_data = self
            .loops
            .get(loop_to_move)
            .ok_or_else(|| CoreError::InvalidHandle("face outer loop handle is invalid".into()))?;
        if loop_data.face != face {
            return Err(CoreError::PreconditionFailed(
                "face outer loop must be owned by the face being killed".into(),
            ));
        }

        self.loops[loop_to_move].face = target_face;
        self.faces[target_face].inner_loops.push(loop_to_move);

        let target_shell = target_data.shell;
        let faces_to_update: Vec<FaceHandle> = self
            .faces
            .iter()
            .filter(|(fh, fd)| fd.shell == killed_shell && *fh != face)
            .map(|(fh, _)| fh)
            .collect();
        for fh in faces_to_update.iter().copied() {
            self.faces[fh].shell = target_shell;
        }

        self.faces.remove(face);
        self.shells.remove(killed_shell);
        self.solids.remove(killed_solid);

        Ok(KfmrhApplied {
            loop_to_move,
            killed_shell,
            killed_solid,
            target_shell,
            reassigned_faces: faces_to_update,
        })
    }

    fn apply_mfkrh(&mut self, inner_loop: LoopHandle) -> Result<MfkrhResult, CoreError> {
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

        let new_face = self.faces.insert(Face {
            outer_loop: inner_loop,
            inner_loops: vec![],
            shell: Default::default(), // patched below
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });

        self.loops[inner_loop].face = new_face;

        let new_shell = self.shells.insert(Shell {
            face: new_face,
            solid: Default::default(), // patched below
        });

        let new_solid = self.solids.insert(Solid { shell: new_shell });

        self.faces[new_face].shell = new_shell;
        self.shells[new_shell].solid = new_solid;

        Ok(MfkrhResult {
            face: new_face,
            shell: new_shell,
        })
    }

    fn verify_kfmrh_postconditions(
        &self,
        before: EulerTopologyCounts,
        face: FaceHandle,
        target_face: FaceHandle,
        applied: &KfmrhApplied,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "kfmrh",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::kfmrh(),
        )?;
        if self.faces.contains_key(face) {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: killed face remains live".into(),
            ));
        }
        if self.shells.contains_key(applied.killed_shell) {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: killed shell remains live".into(),
            ));
        }
        if self.solids.contains_key(applied.killed_solid) {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: killed solid remains live".into(),
            ));
        }
        let target = self
            .faces
            .get(target_face)
            .ok_or_else(|| CoreError::PreconditionFailed("target face is not live".into()))?;
        if target.shell != applied.target_shell {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: target face changed shells".into(),
            ));
        }
        if !target.inner_loops.contains(&applied.loop_to_move) {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: target face does not own moved ring loop".into(),
            ));
        }
        let moved_loop = self
            .loops
            .get(applied.loop_to_move)
            .ok_or_else(|| CoreError::PreconditionFailed("moved loop is not live".into()))?;
        if moved_loop.face != target_face {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: moved ring loop points at wrong face".into(),
            ));
        }
        for moved_face in &applied.reassigned_faces {
            let moved = self
                .faces
                .get(*moved_face)
                .ok_or_else(|| CoreError::PreconditionFailed("moved face is not live".into()))?;
            if moved.shell != applied.target_shell {
                return Err(CoreError::PreconditionFailed(
                    "kfmrh postcondition failed: moved face has wrong target shell".into(),
                ));
            }
        }
        self.verify_no_references_to_killed_face_shell_solid(
            face,
            applied.killed_shell,
            applied.killed_solid,
        )?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        Ok(())
    }

    fn verify_mfkrh_postconditions(
        &self,
        before: EulerTopologyCounts,
        inner_loop: LoopHandle,
        owner_face: FaceHandle,
        result: MfkrhResult,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "mfkrh",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::mfkrh(),
        )?;
        let old_owner = self
            .faces
            .get(owner_face)
            .ok_or_else(|| CoreError::PreconditionFailed("old owner face is not live".into()))?;
        if old_owner.inner_loops.contains(&inner_loop) {
            return Err(CoreError::PreconditionFailed(
                "mfkrh postcondition failed: old owner still lists promoted loop".into(),
            ));
        }
        let new_face = self
            .faces
            .get(result.face)
            .ok_or_else(|| CoreError::PreconditionFailed("new face is not live".into()))?;
        if new_face.outer_loop != inner_loop || !new_face.inner_loops.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "mfkrh postcondition failed: new face does not use promoted outer loop".into(),
            ));
        }
        if new_face.shell != result.shell {
            return Err(CoreError::PreconditionFailed(
                "mfkrh postcondition failed: new face has wrong shell".into(),
            ));
        }
        let loop_data = self
            .loops
            .get(inner_loop)
            .ok_or_else(|| CoreError::PreconditionFailed("promoted loop is not live".into()))?;
        if loop_data.face != result.face {
            return Err(CoreError::PreconditionFailed(
                "mfkrh postcondition failed: promoted loop points at wrong face".into(),
            ));
        }
        let shell = self
            .shells
            .get(result.shell)
            .ok_or_else(|| CoreError::PreconditionFailed("new shell is not live".into()))?;
        if shell.face != result.face {
            return Err(CoreError::PreconditionFailed(
                "mfkrh postcondition failed: new shell does not reference new face".into(),
            ));
        }
        let solid = self
            .solids
            .get(shell.solid)
            .ok_or_else(|| CoreError::PreconditionFailed("new solid is not live".into()))?;
        if solid.shell != result.shell {
            return Err(CoreError::PreconditionFailed(
                "mfkrh postcondition failed: new solid does not reference new shell".into(),
            ));
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        Ok(())
    }

    fn verify_no_references_to_killed_face_shell_solid(
        &self,
        face: FaceHandle,
        shell: ShellHandle,
        solid: crate::topology::SolidHandle,
    ) -> Result<(), CoreError> {
        if self
            .loops
            .iter()
            .any(|(_, loop_data)| loop_data.face == face)
        {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: a loop still references the killed face".into(),
            ));
        }
        if self
            .faces
            .iter()
            .any(|(_, face_data)| face_data.shell == shell)
        {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: a face still references the killed shell".into(),
            ));
        }
        if self
            .shells
            .iter()
            .any(|(_, shell_data)| shell_data.solid == solid)
        {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: a shell still references the killed solid".into(),
            ));
        }
        if self
            .solids
            .iter()
            .any(|(_, solid_data)| solid_data.shell == shell)
        {
            return Err(CoreError::PreconditionFailed(
                "kfmrh postcondition failed: a solid still references the killed shell".into(),
            ));
        }
        Ok(())
    }

    fn loop_owner_face(&self, loop_: LoopHandle) -> Result<FaceHandle, CoreError> {
        self.loops
            .get(loop_)
            .map(|loop_data| loop_data.face)
            .ok_or_else(|| CoreError::InvalidHandle("inner_loop handle is invalid".into()))
    }
}

#[cfg(test)]
#[path = "kfmrh_tests.rs"]
mod tests;
