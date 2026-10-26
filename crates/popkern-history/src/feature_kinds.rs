//! Operation kinds shared by several features: the boolean set
//! operation and whether a feature makes a new body or merges with its
//! host. Re-exported from [`crate::feature`].

use serde::{Deserialize, Serialize};

/// Boolean set-operation kind.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum BoolOpKind {
    /// Unite (add) two bodies.
    Union,
    /// Retain only the intersection of two bodies.
    Intersection,
    /// Subtract the tool body from the target.
    Subtraction,
}

/// Controls whether an operation creates new geometry or merges with the
/// active body.
///
/// `NewBody` is the default and creates an independent solid.  `Boss` and
/// `Cut` automatically union or subtract from the nearest host body,
/// mirroring the behaviour of direct-modelling "Boss" and "Cut" operations.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub enum OperationMode {
    /// Create a new, independent body (default).
    #[default]
    NewBody,
    /// Unite the result with the host body (additive).
    Boss,
    /// Subtract the result from the host body (subtractive).
    Cut,
}
