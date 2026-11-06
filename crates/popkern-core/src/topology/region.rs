use crate::topology::ShellHandle;
use serde::{Deserialize, Serialize};

/// What fills a region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum RegionKind {
    /// Material: a body's solid interior.
    #[default]
    Solid,
    /// Empty space: the unbounded outside, or a cavity.
    Void,
}

/// A region: a connected 3-cell of space, solid or void, bounded by
/// shells. Every body has exactly one unbounded void region, its
/// outside ([`crate::topology::BRep::outside`]), which has no peripheral
/// shell; every other region is bounded and has one.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Region {
    /// What fills the region. Bodies written before regions were
    /// stored held solid regions only.
    #[serde(default)]
    pub kind: RegionKind,
    /// The peripheral (outer) shell: the boundary between this region
    /// and the region enclosing it. Null for the outside.
    pub shell: ShellHandle,
}
