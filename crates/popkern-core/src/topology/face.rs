use crate::geometry::SurfaceOnFace;
use crate::naming::NameId;
use crate::topology::{LoopHandle, ShellHandle, TopologyRole};
use serde::{Deserialize, Serialize};

/// A topological face — bounded by loops.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face {
    /// Outer boundary loop.
    pub outer_loop: LoopHandle,
    /// Inner loops (holes).
    pub inner_loops: Vec<LoopHandle>,
    /// The shell holding this face's back side: the side its effective
    /// normal points away from, a solid's material for a solid's face.
    #[serde(alias = "shell")]
    pub back: ShellHandle,
    /// The shell holding this face's front side: the void across a
    /// solid's face, the same shell as `back` for a sheet's face. Bodies
    /// written before face sides were stored load with it null and have
    /// it derived ([`crate::topology::BRep::derive_regions`]).
    #[serde(default)]
    pub front: ShellHandle,
    /// Geometric surface.
    pub surface: Option<SurfaceOnFace>,
    /// Stable topological name for this face.
    #[serde(default)]
    pub name: Option<NameId>,
    /// Per-entity tolerance (accuracy of this face's surface geometry).
    /// `None` means use the BRep-level default tolerance.
    #[serde(default)]
    pub tolerance: Option<f64>,
    /// Whether this face belongs to production topology or retained construction residue.
    #[serde(default)]
    pub role: TopologyRole,
}

impl Face {
    /// The shells this face's two sides lie in.
    pub fn sides(&self) -> crate::topology::FaceSides {
        crate::topology::FaceSides {
            back: self.back,
            front: self.front,
        }
    }
}
