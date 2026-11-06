use crate::topology::{FaceHandle, RegionHandle};
use serde::{Deserialize, Serialize};

/// A shell: a connected set of face sides bounding one region. A
/// closed solid's outer boundary is two shells, the solid region's
/// (holding its faces' back sides) and the void's across it (their
/// front sides); a sheet is one shell holding both sides of its faces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    /// A face with a side in this shell (the shell's anchor).
    pub face: FaceHandle,
    /// The region this shell bounds.
    #[serde(alias = "solid")]
    pub region: RegionHandle,
}
