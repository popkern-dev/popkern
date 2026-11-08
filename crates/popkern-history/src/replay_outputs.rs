//! The bodies a feature serving several bodies writes (PK-1120).
//!
//! A split, or a boolean whose touching pieces are served under the
//! split policy, serves several bodies: no body holds lumps that touch
//! without sharing (PK-1115). Each piece is served as a body of its own,
//! the feature's own body one of them. The node records each body with
//! a key, a face name only its piece carried when it was assigned, and
//! a later replay serves the piece carrying that name as that body, so
//! a piece keeps its body when an upstream edit reorders the pieces or
//! changes how many there are. The assignment reads only the recorded
//! keys and the pieces, never an earlier replay's result: a reloaded
//! document assigns the same bodies.
//!
//! The feature's own body always takes a piece (its key's, else the
//! first piece no key claimed), so features after it keep a body to act
//! on. A recorded body whose key no piece carries (its piece's faces
//! renamed) takes the next unclaimed piece in order, recorded as a
//! fallback, so the bodies hold while the count holds; one left over
//! when pieces are gone is absent from this replay and keeps its key,
//! so it returns with its piece; a piece left over is served as a new
//! body, allocated and recorded. Keys prefer names with lineage, which
//! follow their piece, over generated names numbered by position.

use std::collections::{BTreeSet, HashSet};

use popkern_core::naming::{EntityName, NameId, RoleCategory};
use popkern_core::topology::BRep;

use crate::FeatureTree;
use crate::node::{BodyId, FeatureId, OutputBody};
use crate::replay_report::{FallbackRecord, FeatureEvidence, RecordSource};
use crate::working_bodies::WorkingBodies;

impl FeatureTree {
    /// Serve feature `id`'s `pieces` as its bodies ([`assign`]): every
    /// piece but its own body's into `bodies`, a body whose piece is
    /// gone out of it, each noted in `touched`, new bodies allocated
    /// and the node's record updated, a body matched by order rather
    /// than by its key recorded as a fallback. Returns its own body's piece
    /// (always one when there are pieces).
    pub(crate) fn serve_outputs(
        &mut self,
        id: FeatureId,
        pieces: Vec<BRep>,
        bodies: &mut WorkingBodies,
        touched: &mut HashSet<BodyId>,
        evidence: &mut FeatureEvidence,
    ) -> Option<BRep> {
        if pieces.is_empty() {
            return None;
        }
        let own = self.features[id].body;
        let recorded = self.features[id].output_bodies.clone();
        let manager = &mut self.body_manager;
        let mut allocate = || {
            let body = manager.create_body();
            manager.assign_feature(body, id);
            body
        };
        let assigned = assign(own, &recorded, pieces, &mut allocate);
        for gone in assigned.absent {
            bodies.remove(&gone);
            touched.insert(gone);
        }
        let mut own_piece = None;
        for (body, piece) in assigned.served {
            if body == own {
                own_piece = Some(piece);
            } else {
                bodies.insert(body, piece);
                touched.insert(body);
            }
        }
        for body in assigned.by_order {
            evidence.fallbacks.push(FallbackRecord {
                source: RecordSource::Naming,
                mechanism: "output body key name -> piece order".into(),
                reason: format!(
                    "no piece carries the face name body {} was keyed by; it takes the next \
                     unclaimed piece",
                    body.0
                ),
            });
        }
        self.features[id].output_bodies = assigned.outputs;
        own_piece
    }
}

/// Served pieces assigned to the feature's bodies.
pub(crate) struct Assigned {
    /// Each served piece with its body, the feature's own body first.
    pub served: Vec<(BodyId, BRep)>,
    /// The node's output bodies after this assignment (empty while the
    /// feature has only ever served one body).
    pub outputs: Vec<OutputBody>,
    /// Recorded bodies no piece is served as on this replay.
    pub absent: Vec<BodyId>,
    /// Recorded bodies whose key no piece carried, served the next
    /// unclaimed piece in order instead.
    pub by_order: Vec<BodyId>,
}

/// Assign `pieces` (in the order the operation served them) to the
/// bodies `recorded` for a feature whose own body is `primary`,
/// allocating a body for each piece no recorded key claims.
pub(crate) fn assign(
    primary: BodyId,
    recorded: &[OutputBody],
    pieces: Vec<BRep>,
    allocate: &mut dyn FnMut() -> BodyId,
) -> Assigned {
    let names: Vec<BTreeSet<u32>> = pieces.iter().map(face_names).collect();
    let carries = |i: usize, key: &Option<EntityName>| {
        key.as_ref()
            .and_then(|k| pieces[i].name_pool.find(k))
            .is_some_and(|id| names[i].contains(&id.0))
    };
    let mut owner: Vec<Option<usize>> = vec![None; pieces.len()];
    let mut fallbacks = Vec::new();
    let mut outputs: Vec<OutputBody> = recorded.to_vec();
    if outputs.is_empty() {
        outputs.push(OutputBody {
            body: primary,
            key: None,
        });
    }
    // The feature's own body first, then the others, each by its key.
    let order: Vec<usize> = (0..outputs.len())
        .filter(|&k| outputs[k].body == primary)
        .chain((0..outputs.len()).filter(|&k| outputs[k].body != primary))
        .collect();
    for &k in &order {
        if let Some(i) =
            (0..pieces.len()).find(|&i| owner[i].is_none() && carries(i, &outputs[k].key))
        {
            owner[i] = Some(k);
        }
    }
    let own = order.first().copied();
    if let Some(p) = own
        && !owner.contains(&Some(p))
    {
        // The own body's key is gone: the first unclaimed piece, else
        // the first piece (whose body is then absent).
        let i = (0..pieces.len())
            .find(|&i| owner[i].is_none())
            .or((!pieces.is_empty()).then_some(0));
        if let Some(i) = i {
            owner[i] = Some(p);
            if outputs[p].key.is_some() {
                fallbacks.push(outputs[p].body);
            }
        }
    }
    // A recorded body whose key no piece carries (its piece's faces
    // renamed: a name of no lineage is numbered by position) takes the
    // next unclaimed piece in order, so the bodies hold while the count
    // holds; only a piece left over is a new body.
    for &k in order.iter().skip(1) {
        if owner.contains(&Some(k)) {
            continue;
        }
        let Some(i) = (0..pieces.len()).find(|&i| owner[i].is_none()) else {
            break;
        };
        owner[i] = Some(k);
        if outputs[k].key.is_some() {
            fallbacks.push(outputs[k].body);
        }
    }
    for slot in owner.iter_mut().filter(|o| o.is_none()) {
        *slot = Some(outputs.len());
        outputs.push(OutputBody {
            body: allocate(),
            key: None,
        });
    }
    // Keys are fixed once recorded; a body assigned now without one
    // takes its piece's.
    for (i, k) in owner.iter().enumerate() {
        if let Some(k) = *k
            && outputs[k].key.is_none()
        {
            outputs[k].key = key_of(i, &pieces, &names);
        }
    }
    let mut served: Vec<(usize, BodyId, BRep)> = pieces
        .into_iter()
        .zip(&owner)
        .filter_map(|(piece, k)| k.map(|k| (k, outputs[k].body, piece)))
        .collect();
    served.sort_by_key(|(k, body, _)| (*body != primary, *k));
    let absent = outputs
        .iter()
        .enumerate()
        .filter(|(k, _)| !owner.contains(&Some(*k)))
        .map(|(_, o)| o.body)
        .collect();
    if outputs.len() == 1 && recorded.is_empty() {
        outputs.clear();
    }
    Assigned {
        served: served.into_iter().map(|(_, b, p)| (b, p)).collect(),
        outputs,
        absent,
        by_order: fallbacks,
    }
}

/// The names of `piece`'s faces, by id in its pool.
fn face_names(piece: &BRep) -> BTreeSet<u32> {
    piece
        .faces
        .values()
        .filter_map(|f| f.name.map(|n| n.0))
        .collect()
}

/// The key of piece `i`: its least face name no other piece carries,
/// else its least face name, names with lineage (an ancestor, or a role
/// other than generated) before generated names of none, which are
/// numbered by position and move between pieces when they reorder.
fn key_of(i: usize, pieces: &[BRep], names: &[BTreeSet<u32>]) -> Option<EntityName> {
    let pool = &pieces[i].name_pool;
    let positional = |id: u32| {
        pool.get(NameId(id))
            .is_none_or(|n| n.ancestors.is_empty() && n.role.category == RoleCategory::Generated)
    };
    let unique = |id: u32| {
        let name = pool.get(NameId(id));
        names.iter().enumerate().all(|(j, other)| {
            j == i
                || name
                    .and_then(|n| pieces[j].name_pool.find(n))
                    .is_none_or(|o| !other.contains(&o.0))
        })
    };
    let ranked = |id: &u32| (!unique(*id), positional(*id), *id);
    let id = names[i].iter().min_by_key(|id| ranked(id))?;
    pool.get(NameId(*id)).cloned()
}

#[cfg(test)]
#[path = "replay_outputs_tests.rs"]
mod tests;
