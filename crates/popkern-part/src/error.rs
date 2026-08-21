//! Part-level errors: each kernel error is preserved intact and
//! wrapped with the operation context, so a GUI (or a script author)
//! sees "what you asked for, what the limit was" — not a flattened
//! string.

use thiserror::Error;

/// Errors from part modeling operations.
///
/// Kernel error payloads are carried whole; `Display` renders a
/// human-readable message per the UI-forward rules.
#[derive(Debug, Error)]
pub enum PartError {
    /// Building a smooth body from a control mesh failed.
    ///
    /// Carries the pipeline stage that declined so the caller can
    /// tell a bad cage from a weld that could not close.
    #[error("smooth mesh failed at {stage}: {detail}")]
    SmoothMesh {
        /// Which stage declined ("obj parse", "spline construction",
        /// "boundary weld", "closure").
        stage: &'static str,
        /// The underlying reason.
        detail: String,
    },
    /// Profile extrusion / sweeping failed in the sweep kernel.
    #[error("sweep failed while building '{feature}': {source}")]
    Sweep {
        /// The feature being built (e.g. "extrude rect 2.4×1.6").
        feature: String,
        /// The kernel error.
        source: popkern_sweep::SweepError,
    },
    /// A boolean combination failed.
    #[error("boolean {operation} failed: {source}")]
    Boolean {
        /// Which combination ("union", "subtract", "intersect").
        operation: &'static str,
        /// The kernel error.
        source: popkern_bool::BoolError,
    },
    /// An edge blend (fillet/chamfer) failed.
    #[error("fillet failed on selection '{selection}': {source}")]
    Blend {
        /// Human-readable description of the resolved selection.
        selection: String,
        /// The kernel error.
        source: popkern_fillet::BlendError,
    },
    /// A selector matched nothing. Silent no-op dress-up operations
    /// are how modeling scripts rot, so empty selections are errors
    /// unless the query opts out with
    /// [`EdgeQuery::allow_empty`](crate::select::EdgeQuery::allow_empty).
    #[error(
        "selection '{selection}' matched no {entity} (part has {available}); \
         call .allow_empty() if a no-op is intended"
    )]
    EmptySelection {
        /// The query description.
        selection: String,
        /// What was being selected ("edges", "faces").
        entity: &'static str,
        /// How many candidates existed before filtering.
        available: usize,
    },
    /// The result of an operation failed kernel validation; the
    /// input part is unchanged.
    #[error("'{feature}' produced an invalid body: {source}")]
    Validation {
        /// The feature that produced the body.
        feature: String,
        /// The kernel validation error.
        source: popkern_core::error::ValidationError,
    },
    /// A core kernel operation failed.
    #[error("kernel error in '{feature}': {source}")]
    Core {
        /// The feature being built.
        feature: String,
        /// The kernel error.
        source: popkern_core::error::CoreError,
    },
    /// Tessellation failed.
    #[error("tessellation failed: {source}")]
    Mesh {
        /// The kernel error.
        source: popkern_mesh::MeshError,
    },
    /// The watertightness gate failed: the tessellated part has open
    /// (single-use) mesh edges.
    #[error("part is not watertight: {open_edges} open mesh edges")]
    NotWatertight {
        /// Number of open (boundary) edges in the welded mesh.
        open_edges: usize,
    },
    /// Export or import failed.
    #[error("I/O failed for '{path}': {source}")]
    Io {
        /// The file involved.
        path: String,
        /// The kernel error.
        source: popkern_io::error::IoError,
    },
    /// The operation does not apply to this part's geometry (e.g.
    /// engraving a non-planar face). The message names the entity
    /// and the reason.
    #[error("operation not applicable: {0}")]
    NotApplicable(String),
    /// A cancellation token fired (PK-276, §17-style clean
    /// termination): the operation stopped between stages; nothing was
    /// produced and the caller's inputs are untouched.
    #[error("part operation cancelled by token")]
    Cancelled,
}
