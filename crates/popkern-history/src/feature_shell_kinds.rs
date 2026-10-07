//! Serializable shell / offset option kinds for the feature tree
//! (split from `feature.rs` under the 500-line cap; re-exported there,
//! so `crate::feature::ShellDirectionKind` and friends keep their
//! paths and serialized form).

use serde::{Deserialize, Serialize};

/// Shell offset direction.
///
/// This mirrors [`popkern_shell::types::ShellDirection`] but adds
/// `Serialize`/`Deserialize` so the feature tree can be persisted.
/// Convert to/from `ShellDirection` at the replay boundary.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ShellDirectionKind {
    /// Offset inward (remove material from the inside).
    Inward,
    /// Offset outward (add material to the outside).
    Outward,
    /// Split the thickness equally: half inward, half outward.
    Midplane,
}

/// How gap edges of an offset skin are closed (PK-328).
///
/// Mirrors [`popkern_shell::types::JoinStyle`] with serde; converted
/// at the replay boundary.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub enum JoinStyleKind {
    /// Rolling-ball fills (Minkowski-faithful; the default).
    #[default]
    Arc,
    /// Extend and intersect the adjacent offsets; sharp result.
    Intersection,
}

/// How globally self-intersecting offset regions are handled (PK-328).
///
/// Mirrors [`popkern_shell::types::TrimPolicy`] with serde; converted
/// at the replay boundary.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub enum TrimPolicyKind {
    /// Trim the colliding regions and re-close (the default).
    #[default]
    Trim,
    /// Decline with diagnostics instead of trimming.
    Decline,
}
