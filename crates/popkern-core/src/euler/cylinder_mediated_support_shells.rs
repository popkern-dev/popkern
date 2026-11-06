//! Shell-membership preconditions of a cylinder-mediated support
//! continuation: the support, outside and cylinder faces and every
//! context owner face lie in one shell.

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle};

impl BRep {
    pub(super) fn require_same_shell(
        &self,
        support_face: FaceHandle,
        outside_face: FaceHandle,
        cylinder_face: FaceHandle,
    ) -> Result<(), CoreError> {
        let shell = self.faces[support_face].back;
        if self.faces[outside_face].back != shell || self.faces[cylinder_face].back != shell {
            return Err(CoreError::PreconditionFailed(format!(
                "cylinder-mediated support continuation requires one shell: support={:?} outside={:?} cylinder={:?}",
                self.faces[support_face].back,
                self.faces[outside_face].back,
                self.faces[cylinder_face].back
            )));
        }
        Ok(())
    }

    pub(super) fn require_context_owner_faces(
        &self,
        support_face: FaceHandle,
        outside_face: FaceHandle,
        cylinder_face: FaceHandle,
        context_owner_faces: &[FaceHandle],
    ) -> Result<(), CoreError> {
        let shell = self.faces[support_face].back;
        for &face in context_owner_faces {
            self.require_face(face, "context owner face")?;
            if face == support_face || face == outside_face || face == cylinder_face {
                continue;
            }
            if self.faces[face].back != shell {
                return Err(CoreError::PreconditionFailed(format!(
                    "context owner face {face:?} belongs to shell {:?}, expected {shell:?}",
                    self.faces[face].back
                )));
            }
        }
        Ok(())
    }
}
