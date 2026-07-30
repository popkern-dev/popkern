//! Anchor transplantation for the chain-holder targeting fix
//! (PK-141): move a pair of valence-1 spur-tip anchors onto the
//! face that carries the fillet's own chain descendants.

use popkern_core::topology::BRep;

/// Transplant a pair of valence-1 spur-tip anchors onto `holder`:
/// kev each spur and re-insert the anchor at the same position with
/// a spur rooted on the holder's loops. All-or-nothing — on any
/// failure the candidate work is discarded by the caller re-cloning
/// (the pass runs on the apply candidate). Returns the new anchor
/// vertices `(start, end)`.
pub(super) fn transplant_spur_anchor_pair(
    brep: &mut BRep,
    start_v: popkern_core::topology::VertexHandle,
    end_v: popkern_core::topology::VertexHandle,
    holder: popkern_core::topology::FaceHandle,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<(
    popkern_core::topology::VertexHandle,
    popkern_core::topology::VertexHandle,
)> {
    // Transactional: probe the whole transplant on a candidate
    // clone and commit only when BOTH spurs re-land. The former
    // in-place sequence kev'd a spur before knowing whether the
    // holder accepts the re-insert (planar holders decline the
    // curved-support spur), silently destroying both anchors and
    // leaving the caller holding dangling vertex handles (measured:
    // SlotMap panic on the skewed-rib junction trim, PK-185).
    let mut candidate = brep.clone();
    let mut out = Vec::with_capacity(2);
    for v in [start_v, end_v] {
        let pos = candidate.vertices.get(v)?.position;
        // Valence-1 spur tip: exactly one outgoing half-edge.
        let outgoing: Vec<_> = candidate
            .half_edges
            .iter()
            .filter(|(_, he)| he.origin == v)
            .map(|(h, _)| h)
            .collect();
        let [he_out] = outgoing.as_slice() else {
            if std::env::var_os("PK_MULTICHORD_TRACE").is_some() {
                eprintln!(
                    "[mc-target]   transplant declined: {v:?} valence {}",
                    outgoing.len()
                );
            }
            return None;
        };
        let he_in = candidate.half_edges.get(*he_out)?.twin;
        candidate.kev(he_in).ok()?;
        let nv = super::helpers::helpers_anchor_insert::insert_anchor_with_spur(
            &mut candidate,
            pos,
            holder,
            tolerance,
        )
        .ok()?;
        out.push(nv);
    }
    *brep = candidate;
    Some((out[0], out[1]))
}
