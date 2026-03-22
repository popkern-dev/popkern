use crate::topology::{FaceHandle, SolidHandle};
use serde::{Deserialize, Serialize};

/// A connected set of faces forming a closed surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    /// Any face in this shell.
    pub face: FaceHandle,
    /// Parent solid.
    pub solid: SolidHandle,
}
