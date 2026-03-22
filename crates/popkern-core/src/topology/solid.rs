use crate::topology::ShellHandle;
use serde::{Deserialize, Serialize};

/// A topological solid — bounded by a shell.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Solid {
    /// Outer shell.
    pub shell: ShellHandle,
}
