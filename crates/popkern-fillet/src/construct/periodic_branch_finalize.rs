//! Final periodic pcurve-branch reconciliation (PK-141).
//!
//! Late pipeline stages (strip cleanup, junction settles, seam
//! welds) can re-author pcurves on a PERIODIC support whose loop
//! crosses the seam (the boss-on-wall hole ring), leaving mixed
//! period branches: the UV wire self-intersects even though the 3D
//! body is watertight (measured: OCC SelfIntersectingWire plus
//! UnorientableShape on the boss wall cylinder). Reconcile every
//! periodic-surfaced face once, after the last pcurve author.

use popkern_core::topology::BRep;

/// Reconcile period-branch choices on every periodic-surfaced face.
/// Returns how many faces changed.
pub(crate) fn reconcile_all_periodic_pcurve_branches(brep: &mut BRep) -> usize {
    let periodic_faces: Vec<_> = brep
        .faces
        .iter()
        .filter(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| s.geometry.periodicity().0.is_some())
        })
        .map(|(h, _)| h)
        .collect();
    let mut changed = 0usize;
    for face in periodic_faces {
        if brep.reconcile_periodic_pcurve_branches(face) {
            changed += 1;
            if std::env::var_os("PK_RECONCILE_TRACE").is_some() {
                eprintln!("[periodic-final] reconciled period branches on {face:?}");
            }
        }
    }
    changed
}
