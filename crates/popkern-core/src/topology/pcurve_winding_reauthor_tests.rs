//! Winding-honesty pins for full pcurve re-authoring (PK-322).

use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Re-authoring EVERY pcurve from scratch must reproduce winding-honest
/// loops (PK-322). The closed-edge canonical direction bit cannot
/// orient a face-loop traversal by itself — the loop-winding repair
/// postcondition (pcurve_loop_winding_repair) fixes measurable
/// dishonest loops by flipping closed members. Before it, a stripped
/// cylinder re-authored every cap and lateral loop against its
/// orientation flag (the shell/offset assembly path, which copies
/// topology without pcurves, served WindingDishonest hollow bodies).
#[test]
fn full_pcurve_reauthoring_stays_winding_honest() {
    for (name, mut brep) in [
        (
            "cylinder",
            BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap(),
        ),
        ("cone", BRep::cone(1.0, 2.0, Tolerance::default()).unwrap()),
        ("cube", BRep::unit_cube(Tolerance::default()).unwrap()),
    ] {
        let he_handles: Vec<_> = brep.half_edges.keys().collect();
        for heh in he_handles {
            brep.half_edges[heh].pcurve = None;
        }
        brep.recompute_all_pcurves_winding_repaired().unwrap();
        let winding_errors: Vec<_> = brep
            .validate_all()
            .into_iter()
            .filter(|e| matches!(e, crate::error::ValidationError::WindingDishonest { .. }))
            .collect();
        assert!(
            winding_errors.is_empty(),
            "{name}: re-authored pcurves wind against orientation flags: \
             {winding_errors:?}"
        );
    }
}

/// Reversed bodies must ALSO re-author winding-honest: the shell/offset
/// inner-skin path is clone → reverse_body_orientation → assembly →
/// full pcurve re-authoring. Reversal re-purposes every half-edge to
/// traverse the opposite 3D direction, so it must swap each edge's
/// canonical `half_edge` to the twin — otherwise the closed-edge
/// direction rule authors every reversed skin's rim pcurves backward
/// and the loops wind against their (flipped) orientation flags.
#[test]
fn reversed_body_reauthoring_stays_winding_honest() {
    for (name, mut brep) in [
        (
            "cylinder",
            BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap(),
        ),
        ("cone", BRep::cone(1.0, 2.0, Tolerance::default()).unwrap()),
        ("cube", BRep::unit_cube(Tolerance::default()).unwrap()),
    ] {
        brep.reverse_body_orientation();
        let he_handles: Vec<_> = brep.half_edges.keys().collect();
        for heh in he_handles {
            brep.half_edges[heh].pcurve = None;
        }
        brep.recompute_all_pcurves_winding_repaired().unwrap();
        let winding_errors: Vec<_> = brep
            .validate_all()
            .into_iter()
            .filter(|e| matches!(e, crate::error::ValidationError::WindingDishonest { .. }))
            .collect();
        assert!(
            winding_errors.is_empty(),
            "{name} (reversed): re-authored pcurves wind against orientation \
             flags: {winding_errors:?}"
        );
    }
}
