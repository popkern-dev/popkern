//! Rebinding a feature's references by provenance before it replays
//! (PK-677, `popkern_core::naming_provenance`).
//!
//! A recorded reference is an [`EntityName`]. Its role index and
//! ancestor pool ids are only meaningful in the pool the body had when
//! the feature was recorded, and an upstream edit (a feature inserted
//! at the history playhead, PK-666) shifts both. Each node therefore
//! keeps a by-value snapshot of every reference it resolved in its own
//! body ([`ReferenceProvenance`]). Before the feature runs:
//! - when the recorded name still names the snapshot's entity (its
//!   value matches), it is used as recorded;
//! - otherwise the entity the snapshot describes is found by
//!   provenance, the reference is rewritten to its current name for
//!   this execution only, and the row records the fallback;
//! - several candidates are a typed decline
//!   ([`ReplayError::AmbiguousReference`]), never a pick by index;
//! - none leaves the name as recorded, and the ordinary resolution
//!   reports it.
//!
//! A reference into another body (an extrude's sketch face) has no
//! snapshot and is never rebound. A document saved before snapshots
//! existed gets them on its first replay.

use popkern_core::naming::EntityName;
use popkern_core::naming_provenance::{
    EntityKind, NameValue, ProvenanceMatch, Relaxation, kind_of, resolve_by_provenance,
};
use popkern_core::naming_resolve::{Resolution, resolve_name};
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

use crate::error::ReplayError;
use crate::feature::Feature;
use crate::feature_names::for_each_entity_name_mut;
use crate::replay_report::{FallbackRecord, RecordSource};

/// One reference's by-value snapshot.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReferenceProvenance {
    /// The name as the feature records it.
    pub name: EntityName,
    /// The kind of entity it named.
    pub kind: EntityKind,
    /// That entity's name, by value, when the snapshot was taken.
    pub value: NameValue,
}

/// `feature` with its references rebound against `brep`, and the
/// fallbacks taken.
///
/// # Errors
///
/// [`ReplayError::AmbiguousReference`] when a snapshot matches several
/// entities.
pub(crate) fn rebind(
    feature: &Feature,
    provenance: &[ReferenceProvenance],
    brep: &BRep,
) -> Result<(Feature, Vec<FallbackRecord>), ReplayError> {
    let mut out = feature.clone();
    if provenance.is_empty() {
        return Ok((out, Vec::new()));
    }
    let mut fallbacks = Vec::new();
    let mut error = None;
    for_each_entity_name_mut(&mut out, &mut |name| {
        if error.is_some() {
            return;
        }
        let Some(snapshot) = provenance.iter().find(|p| p.name == *name) else {
            return;
        };
        let current = brep.name_pool.find(name);
        if current.and_then(|id| brep.name_pool.value_of(id)).as_ref() == Some(&snapshot.value) {
            return;
        }
        match resolve_by_provenance(&snapshot.value, snapshot.kind, brep) {
            ProvenanceMatch::Unique {
                name: rebound,
                relaxation,
                ..
            } => {
                fallbacks.push(FallbackRecord {
                    source: RecordSource::Naming,
                    mechanism: "recorded name -> provenance rebind".into(),
                    reason: format!(
                        "{} renumbered upstream: op {} {:?} index {} now {} ({})",
                        if current.is_some() {
                            "the recorded name now denotes another entity:"
                        } else {
                            "the recorded name is gone:"
                        },
                        name.operation.0,
                        name.role.category,
                        name.role.index,
                        rebound.role.index,
                        match relaxation {
                            Relaxation::Exact => "pool ids only",
                            Relaxation::OwnIndex => "its role index",
                            Relaxation::AncestorIndices => "role indices in its ancestry",
                        }
                    ),
                });
                *name = rebound;
            }
            ProvenanceMatch::Ambiguous { candidates, .. } => {
                error = Some(ReplayError::AmbiguousReference {
                    name: name.clone(),
                    candidates,
                });
            }
            ProvenanceMatch::NotFound => {}
        }
    });
    match error {
        Some(error) => Err(error),
        None => Ok((out, fallbacks)),
    }
}

/// Keep `provenance` in step with `feature`: drop snapshots of names
/// the feature no longer records, and take one for every recorded name
/// that resolves in `brep` (the body the feature ran on) and has none.
pub(crate) fn refresh(provenance: &mut Vec<ReferenceProvenance>, feature: &Feature, brep: &BRep) {
    let mut names = Vec::new();
    let mut copy = feature.clone();
    for_each_entity_name_mut(&mut copy, &mut |name| names.push(name.clone()));
    provenance.retain(|p| names.contains(&p.name));
    for name in names {
        if provenance.iter().any(|p| p.name == name) {
            continue;
        }
        let kind = match resolve_name(&name, brep) {
            Resolution::Exact(entity) => kind_of(entity),
            Resolution::Split(pieces) => match pieces.first() {
                Some(entity) => kind_of(*entity),
                None => continue,
            },
            Resolution::NotFound => continue,
        };
        let Some(value) = brep
            .name_pool
            .find(&name)
            .and_then(|id| brep.name_pool.value_of(id))
        else {
            continue;
        };
        provenance.push(ReferenceProvenance { name, kind, value });
    }
}

#[cfg(test)]
#[path = "reference_rebind_tests.rs"]
mod tests;
