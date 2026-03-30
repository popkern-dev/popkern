use crate::geometry::SurfaceOnFace;
use crate::naming::NameId;
use crate::topology::{LoopHandle, ShellHandle};
use serde::{Deserialize, Serialize};

/// A topological face — bounded by loops.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face {
    /// Outer boundary loop.
    pub outer_loop: LoopHandle,
    /// Inner loops (holes).
    pub inner_loops: Vec<LoopHandle>,
    /// Parent shell.
    pub shell: ShellHandle,
    /// Geometric surface.
    pub surface: Option<SurfaceOnFace>,
    /// Stable topological name for this face.
    #[serde(default)]
    pub name: Option<NameId>,
}
