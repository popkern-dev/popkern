//! Winding-honesty pins for full pcurve re-authoring (PK-322).

use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Re-authoring EVERY pcurve from scratch must reproduce winding-honest
/// loops (PK-322) without moving a correct flag or running a pcurve
/// against its edge (PK-878): the shell/offset assembly path copies
/// topology without pcurves and re-authors all of them.
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
        let flags: Vec<_> = brep
            .faces
            .values()
            .map(|f| f.surface.as_ref().map(|s| s.orientation))
            .collect();
        brep.recompute_all_pcurves_winding_repaired().unwrap();
        // The flags were right before the strip: honest re-authoring
        // must reproduce windings that agree with them, not move them
        // (PK-878: a repair that rewrites either side to agree hides
        // which one lied).
        let after: Vec<_> = brep
            .faces
            .values()
            .map(|f| f.surface.as_ref().map(|s| s.orientation))
            .collect();
        assert_eq!(flags, after, "{name}: re-authoring moved orientation flags");
        let miss = brep
            .half_edges
            .keys()
            .filter_map(|h| brep.measure_half_edge_pcurve_deviation(h))
            .map(|r| r.max_deviation)
            .fold(0.0_f64, f64::max);
        assert!(
            miss < 1e-6,
            "{name}: a pcurve runs against its edge ({miss:e})"
        );
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
        let flags: Vec<_> = brep
            .faces
            .values()
            .map(|f| f.surface.as_ref().map(|s| s.orientation))
            .collect();
        brep.recompute_all_pcurves_winding_repaired().unwrap();
        // The flags were right before the strip: honest re-authoring
        // must reproduce windings that agree with them, not move them
        // (PK-878: a repair that rewrites either side to agree hides
        // which one lied).
        let after: Vec<_> = brep
            .faces
            .values()
            .map(|f| f.surface.as_ref().map(|s| s.orientation))
            .collect();
        assert_eq!(flags, after, "{name}: re-authoring moved orientation flags");
        let miss = brep
            .half_edges
            .keys()
            .filter_map(|h| brep.measure_half_edge_pcurve_deviation(h))
            .map(|r| r.max_deviation)
            .fold(0.0_f64, f64::max);
        assert!(
            miss < 1e-6,
            "{name}: a pcurve runs against its edge ({miss:e})"
        );
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
