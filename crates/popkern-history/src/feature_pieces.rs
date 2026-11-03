//! Which pieces of a split body a feature keeps (design: operations on
//! the cell complex, 5.4).

use serde::{Deserialize, Serialize};

/// Which pieces of a [`crate::feature::Feature::Split`] the body keeps,
/// recorded as data so a replay with moved tools selects again.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub enum PieceRule {
    /// Every piece.
    #[default]
    All,
    /// The pieces whose membership in the split's tools (tool k as body
    /// k) the expression accepts.
    InsideTools(popkern_bool::SetExpr),
    /// The piece containing a point (model units).
    Containing([f64; 3]),
    /// The `n` largest pieces by volume.
    Largest(usize),
}
