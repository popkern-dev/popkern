//! Blend-tail dirty-region derivation and the whole-model stale
//! cross-check (PK-274).
//!
//! The pipeline tail's pcurve passes are priced by the blend
//! neighbourhood, not the whole model: the operation only invalidates
//! pcurves on the faces it created plus the retained faces it
//! rewired, and every rewired retained face borders a new face
//! (Euler surgery), so the vertex-sharing one-ring of the new faces
//! covers them. See docs/architecture/region-scoped-validation.md,
//! companion section.

use std::collections::HashSet;

use popkern_core::topology::{BRep, FaceHandle};

/// The tail region: faces created by the blend (absent from
/// `pre_faces`) dilated by the vertex-sharing one-ring. Derived
/// fresh at each tail pass — intermediate repair stages can delete
/// or replace faces.
pub(crate) fn blend_tail_region(brep: &BRep, pre_faces: &HashSet<FaceHandle>) -> Vec<FaceHandle> {
    let blend_dirty: Vec<FaceHandle> = brep
        .faces
        .keys()
        .filter(|f| !pre_faces.contains(f))
        .collect();
    brep.validation_region(&blend_dirty)
}

/// First tail pass: refill and de-stale the pcurves of the blend
/// region (recompute + deviation-gated repair, both scoped).
pub(crate) fn recompute_and_fix_tail(brep: &mut BRep, pre_faces: &HashSet<FaceHandle>) {
    let region = blend_tail_region(brep, pre_faces);
    let _ = brep.recompute_pcurves_for(&region);
    brep.fix_stale_pcurves_among(&region);
}

/// Final tail fidelity pass (PK-095) over a freshly derived region —
/// intermediate stages can delete or replace faces — followed by the
/// row-407 cross-check.
pub(crate) fn final_fix_and_cross_check(
    brep: &mut BRep,
    pre_faces: &HashSet<FaceHandle>,
    entry_stale: &[popkern_core::topology::HalfEdgeHandle],
) {
    let region = blend_tail_region(brep, pre_faces);
    brep.fix_stale_pcurves_among(&region);
    debug_assert_tail_region_complete(brep, &region, entry_stale);
}

/// Debug-only snapshot of the half-edges whose pcurves are ALREADY
/// stale at pipeline entry. Input bodies can carry staleness the
/// blend never caused (measured: open-box builder caps, closed-rim
/// modify/unfillet round-trip bodies — 14 debug fixtures); the old
/// whole-model tail pass silently healed those as a side effect, and
/// the cross-check must not blame the region derivation for them.
/// Returns empty in release builds (the whole-model scan is the cost
/// the scoping exists to remove).
pub(crate) fn entry_stale_snapshot(brep: &BRep) -> Vec<popkern_core::topology::HalfEdgeHandle> {
    #[cfg(debug_assertions)]
    {
        brep.stale_pcurve_handles()
    }
    #[cfg(not(debug_assertions))]
    {
        let _ = brep;
        Vec::new()
    }
}

/// PK-274 row-407 cross-check: a NEW stale pcurve OUTSIDE the scoped
/// tail region means a dirty-region derivation (the pipeline's or a
/// consolidation pass's) lost a face — that class ships as a silent
/// mesher/export defect, so debug builds fail loudly here instead.
/// Half-edges already stale at entry ([`entry_stale_snapshot`]) are
/// exempt: their staleness predates the operation. Compiled out of
/// release builds.
#[cfg(debug_assertions)]
pub(crate) fn debug_assert_tail_region_complete(
    brep: &BRep,
    tail_region: &[FaceHandle],
    entry_stale: &[popkern_core::topology::HalfEdgeHandle],
) {
    let region: HashSet<FaceHandle> = tail_region.iter().copied().collect();
    let entry: HashSet<popkern_core::topology::HalfEdgeHandle> =
        entry_stale.iter().copied().collect();
    for heh in brep.stale_pcurve_handles() {
        if entry.contains(&heh) {
            continue;
        }
        let face = brep
            .half_edges
            .get(heh)
            .and_then(|he| brep.loops.get(he.loop_))
            .map(|l| l.face);
        debug_assert!(
            face.is_some_and(|f| region.contains(&f)),
            "PK-274: NEW stale pcurve on {heh:?} (face {face:?}) escaped the scoped tail \
             region — the dirty-region derivation missed a face"
        );
    }
}

#[cfg(not(debug_assertions))]
pub(crate) fn debug_assert_tail_region_complete(
    _brep: &BRep,
    _tail_region: &[FaceHandle],
    _entry_stale: &[popkern_core::topology::HalfEdgeHandle],
) {
}
