use serde::{Deserialize, Serialize};

/// Semantic ownership role for topology entities.
///
/// Production topology participates in exported model validation and user-visible
/// diagnostics. Construction residual topology is retained as an internal record
/// of retired or replaced material and must not be promoted back into production
/// by late repair passes.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TopologyRole {
    /// User-visible topology that defines the modeled shape.
    #[default]
    Production,
    /// Internal construction residue retained after a validated replacement.
    ConstructionResidual,
}
