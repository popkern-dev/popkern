//! The uniform failure a 2D corner operation (fillet, chamfer) reports.

use popkern_core::topology::{BRep, VertexHandle};

use crate::error::{BlendError, BlendFailure};

/// Build a uniform `PartialFailure` wrapping an `InvalidTopology`
/// blend failure for a 2D corner operation. Shared between
/// fillet and chamfer so both emit the same error shape.
///
/// The location is looked up from `brep` when `vertex` resolves;
/// callers that fail outside the context of a specific corner
/// (e.g. chain-walk failures in `edges_2d`) pass
/// `VertexHandle::default()`, which yields `location: None`.
pub(crate) fn topology_failure(brep: &BRep, vertex: VertexHandle, reason: String) -> BlendError {
    let location = brep.vertices.get(vertex).map(|v| v.position);
    BlendError::partial_failure(
        1,
        BRep::new(popkern_core::tolerance::Tolerance::default()),
        vec![(
            popkern_core::topology::EdgeHandle::default(),
            BlendFailure::InvalidTopology {
                reason: format!("2D corner at {vertex:?}: {reason}"),
                location,
            },
        )],
    )
}
