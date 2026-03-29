//! Error and warning types for sketch-to-BRep conversion.

/// Errors that can occur during sketch-to-BRep conversion.
#[derive(Debug, thiserror::Error)]
pub enum SketchError {
    /// The face surface is not planar.
    #[error("face surface is not planar")]
    NonPlanarFace,

    /// The three points used to define a plane are collinear.
    #[error("collinear points cannot define a plane")]
    CollinearPoints,

    /// The computed normal vector is degenerate (zero length).
    #[error("degenerate normal vector")]
    DegenerateNormal,

    /// A profile edge crosses itself at the given 2D position.
    #[error("self-intersecting profile at ({0}, {1})")]
    SelfIntersectingProfile(f64, f64),

    /// No closed loops could be found in the sketch.
    #[error("no closed loops found in sketch")]
    NoClosedLoops,

    /// A topological invariant was violated during BRep construction.
    #[error("invalid topology: {0}")]
    InvalidTopology(String),

    /// An error propagated from the constraint solver.
    #[error(transparent)]
    ConstraintError(#[from] popkern_constraint::ConstraintError),

    /// An error propagated from the BRep kernel.
    #[error(transparent)]
    CoreError(#[from] popkern_core::error::CoreError),
}

/// Non-fatal warnings produced during sketch-to-BRep conversion.
#[derive(Debug, Clone)]
pub enum SketchWarning {
    /// A wire path was not closed (endpoints did not meet).
    UnclosedWire,
    /// An edge was degenerate (zero length) and was skipped.
    DegenerateEdge,
    /// A geometric validation warning from the output BRep.
    ValidationWarning,
}
