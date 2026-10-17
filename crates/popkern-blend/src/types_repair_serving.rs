//! PK-250/PK-251 serving-boundary records: the exclusion-retry and
//! footprint-repair evidence the fillet pipeline attaches to its
//! repair report at the serving boundary.

use popkern_core::topology::VertexHandle;

/// One PK-250 serving-boundary retry record.
#[derive(Debug, Clone)]
pub struct ServingSelfIntersectionRetry {
    /// Apex vertex whose first-pass direct-apex serve was rejected:
    /// its region self-intersected, its construction recorded a
    /// serving-shape decline, or the served-body gate declined the
    /// body (PK-533).
    pub vertex: VertexHandle,
    /// Why the first pass was rejected: the validator's offending
    /// face pair, the serving-shape reason, or the gate's reason.
    pub detail: String,
}
