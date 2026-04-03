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
    /// Per-entity tolerance (accuracy of this face's surface geometry).
    /// `None` means use the BRep-level default tolerance.
    #[serde(default)]
    pub tolerance: Option<f64>,
}
