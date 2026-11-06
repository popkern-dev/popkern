use crate::topology::{EdgeHandle, FaceHandle, RegionHandle, VertexHandle};
use serde::{Deserialize, Serialize};

/// A shell: a connected set of face sides, wire edges or an isolated
/// vertex bounding one region. A closed solid's outer boundary is two
/// shells, the solid region's (holding its faces' back sides) and the
/// void's across it (their front sides); a sheet is one shell holding
/// both sides of its faces; a wire is one shell of the region it lies
/// in holding its edges; an acorn is one shell holding one vertex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    /// A face with a side in this shell (the shell's anchor); null for a
    /// shell of wire edges or an acorn.
    pub face: FaceHandle,
    /// The region this shell bounds.
    #[serde(alias = "solid")]
    pub region: RegionHandle,
    /// The wire edges (edges on no face) this shell holds, connected
    /// through their vertices or through the shell's face sides.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub wire_edges: Vec<EdgeHandle>,
    /// The isolated vertex this shell holds when it is an acorn shell.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acorn: Option<VertexHandle>,
}

impl Shell {
    /// A shell of `region` anchored at `face`, holding no wire edge or
    /// acorn.
    pub fn new(face: FaceHandle, region: RegionHandle) -> Self {
        Self {
            face,
            region,
            wire_edges: Vec::new(),
            acorn: None,
        }
    }
}
