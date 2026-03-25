//! Public types for shell/offset operations.

use crate::error::ShellWarning;
use popkern_core::topology::BRep;

/// Direction for the shell offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellDirection {
    /// Offset the inner skin inward (material is added outward).
    Inward,
    /// Offset the inner skin outward (material is removed inward).
    Outward,
    /// Split the thickness equally: half inward, half outward.
    Midplane,
}

/// Options controlling shell operation behavior.
#[derive(Debug, Clone)]
pub struct ShellOptions {
    /// When true, faces that fail to offset are auto-removed with a warning
    /// instead of causing the entire operation to fail.
    pub allow_partial: bool,
}

impl Default for ShellOptions {
    fn default() -> Self {
        Self {
            allow_partial: true,
        }
    }
}

/// Result of a successful shell operation.
#[derive(Debug, Clone)]
pub struct ShellResult {
    /// The resulting BRep solid.
    pub brep: BRep,
    /// Any warnings generated during the operation.
    pub warnings: Vec<ShellWarning>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shell_direction_variants() {
        assert_eq!(ShellDirection::Inward, ShellDirection::Inward);
        assert_ne!(ShellDirection::Inward, ShellDirection::Outward);
        let _ = format!("{:?}", ShellDirection::Midplane);
    }

    #[test]
    fn shell_options_default_allows_partial() {
        let opts = ShellOptions::default();
        assert!(opts.allow_partial);
    }

    #[test]
    fn shell_result_debug() {
        use popkern_core::tolerance::Tolerance;
        let result = ShellResult {
            brep: BRep::new(Tolerance::default()),
            warnings: vec![],
        };
        let _ = format!("{:?}", result);
    }
}
