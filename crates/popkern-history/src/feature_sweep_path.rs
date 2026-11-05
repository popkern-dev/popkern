//! A sweep path of several edges (brep-bench §sweep: mitred elbows,
//! picture frames, pipes round a box's edges): the path edge and more
//! edges chained end to end, and what the sweep does where they meet
//! at a corner.

use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};

/// The rest of a sweep's path and its corner treatment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SweepPath {
    /// Edges after the path edge, in any order and direction (chained
    /// from the path edge's end).
    pub more_edges: Vec<EntityName>,
    /// What happens at a sharp corner.
    #[serde(default)]
    pub corner: SweepCorner,
}

/// A sweep's treatment of a sharp corner in its path (where the
/// tangent turns; smooth joints need none).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SweepCorner {
    /// Each run is cut by, and they meet in, the plane bisecting the
    /// corner.
    #[default]
    Mitre,
    /// The section turns about the corner from one run to the next
    /// (a ball's quarter at a pipe's bend).
    Round,
}
