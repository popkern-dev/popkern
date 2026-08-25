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

/// How gap edges — edges whose adjacent offset surfaces separate over
/// their natural windows (convex w.r.t. the displacement direction) —
/// are closed (design §5 Stage 3, PK-325).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum JoinStyle {
    /// Rolling-ball fills: cylinder/torus/canal windows along gap edges,
    /// sphere windows at gap corners. The Minkowski-faithful result — the
    /// shape the SDF offset definition `dist(p, S) ≤ d` describes.
    #[default]
    Arc,
    /// Extend the adjacent offset surfaces and intersect; the sharp
    /// result. Declines per-edge with [`crate::error::ShellError::JoinFailed`]
    /// (naming the arc join as the working alternative) where the
    /// extensions do not meet.
    Intersection,
}

/// Options controlling shell operation behavior.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ShellOptions {
    /// When true, faces that fail to offset are auto-removed with a warning
    /// instead of causing the entire operation to fail.
    pub allow_partial: bool,
    /// How gap edges are closed. Defaults to [`JoinStyle::Arc`]
    /// (Minkowski-faithful). Skins whose gap corners reach a vertex
    /// class the arc machinery does not yet serve (mixed
    /// convex/concave corners, rim vertices — PK-326) downgrade to
    /// [`JoinStyle::Intersection`] with a
    /// [`crate::error::ShellWarning::JoinStyleDowngraded`] warning.
    pub join_style: JoinStyle,
}

impl Default for ShellOptions {
    fn default() -> Self {
        Self {
            allow_partial: true,
            join_style: JoinStyle::default(),
        }
    }
}

impl ShellOptions {
    /// Builder: set the gap-edge join style (`#[non_exhaustive]` blocks
    /// cross-crate struct literals; downstream callers chain from
    /// [`ShellOptions::default`]).
    #[must_use]
    pub fn with_join_style(mut self, join_style: JoinStyle) -> Self {
        self.join_style = join_style;
        self
    }

    /// Builder: set whether failing faces are auto-removed with a
    /// warning instead of failing the whole operation.
    #[must_use]
    pub fn with_allow_partial(mut self, allow_partial: bool) -> Self {
        self.allow_partial = allow_partial;
        self
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
