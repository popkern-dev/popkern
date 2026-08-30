//! [`Feature::RawBody`](crate::Feature::RawBody): a body recorded as a
//! geometry snapshot instead of a recipe.
//!
//! Imports, polyhedra, smooth-mesh solids, and engraved results have no
//! parametric replay in v1; the document embeds the body they produced.
//! Replay materialises the snapshot as a new body, validates it at the
//! serving boundary, and mints every face, edge, and vertex a stable
//! [`EntityName`] under the feature's own
//! [`OperationId`](popkern_core::naming::OperationId) in slotmap order —
//! deterministic across save/load, so downstream features (fillets,
//! shells, transforms, references) address snapshot geometry exactly as
//! they address recipe geometry. Any names the snapshot arrived with are
//! discarded: the snapshot feature is the sole owner of its entities.

use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, NamePool, RoleCategory};
use popkern_core::naming_criteria::{
    GeometricCriteria, capture_edge, capture_face, capture_vertex,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

use crate::error::ReplayError;

/// Where a [`Feature::RawBody`](crate::Feature::RawBody) snapshot came
/// from. Provenance for reports and GUIs; replay does not branch on it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SnapshotOrigin {
    /// Imported from an external body or file (`Part::from_brep`, STEP).
    Import {
        /// Human-readable source: a file name, an interop origin.
        source: String,
    },
    /// Built as a convex polyhedron (`Part::polyhedron`).
    Polyhedron,
    /// Built from a control net by the smooth-mesh family.
    SmoothMesh,
    /// An engraved result (pockets cut into polyhedron faces).
    Engrave,
    /// Any other construction, described.
    Other {
        /// What produced the body.
        summary: String,
    },
}

/// Replace `brep` with the snapshot: validate, then own every entity.
///
/// The snapshot geometry is copied verbatim; only its tolerance is
/// re-stamped with the replay tolerance so downstream operations see
/// one explicit tolerance per body.
///
/// # Errors
///
/// [`ReplayError::FeatureFailed`] when the snapshot has no faces or
/// fails structural/geometric validation — the body is left untouched
/// by the caller's transactional restore.
pub(crate) fn materialize_snapshot(
    brep: &mut BRep,
    snapshot: &BRep,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
) -> Result<(), ReplayError> {
    let failed = |error: String| ReplayError::FeatureFailed {
        label: "raw_body".into(),
        error,
    };
    if snapshot.faces.is_empty() {
        return Err(failed("snapshot has no faces".into()));
    }
    let mut body = snapshot.clone();
    body.tolerance = *tolerance;
    body.validate()
        .map_err(|e| failed(format!("snapshot failed validation: {e}")))?;
    own_all_entities(&mut body, tracking);
    *brep = body;
    Ok(())
}

/// Discard inherited names and mint fresh ones for every entity under
/// the tracking operation: faces `Generated/i`, edges `Lateral/i`,
/// vertices `Junction/i`, each `i` in slotmap key order.
fn own_all_entities(brep: &mut BRep, ctx: &mut TrackingContext) {
    brep.name_pool = NamePool::new();
    let operation = ctx.operation;
    let name = move |category: RoleCategory, index: usize| EntityName {
        operation,
        role: EntityRole {
            category,
            index: index as u32,
        },
        ancestors: Vec::new(),
    };

    let faces: Vec<_> = brep.faces.keys().collect();
    for (i, fh) in faces.into_iter().enumerate() {
        let name = name(RoleCategory::Generated, i);
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[fh].name = Some(id);
        ctx.record_generated(name);
    }
    let edges: Vec<_> = brep.edges.keys().collect();
    for (i, eh) in edges.into_iter().enumerate() {
        let name = name(RoleCategory::Lateral, i);
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.edges[eh].name = Some(id);
        ctx.record_generated(name);
    }
    let vertices: Vec<_> = brep.vertices.keys().collect();
    for (i, vh) in vertices.into_iter().enumerate() {
        let name = name(RoleCategory::Junction, i);
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.vertices[vh].name = Some(id);
        ctx.record_generated(name);
    }
}

#[cfg(test)]
#[path = "feature_snapshot_tests.rs"]
mod tests;
