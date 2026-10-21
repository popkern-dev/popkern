//! PK-906: the pocketed tray's fillet end arcs declare the side their
//! edge window runs on. Edge reconciliation picked the complement of
//! the trig dispatcher's arc and re-windowed the edge onto it, while the
//! curve kept the dispatcher's span; the window check flags any edge
//! left that way.

use popkern_core::error::ValidationError;
use popkern_core::geometry::CurveGeometry;

#[path = "fillet_mileage_corpus_support/parts.rs"]
#[allow(dead_code)]
mod parts;

#[test]
fn the_pocketed_tray_arcs_declare_their_windows_and_a_complement_window_is_flagged() {
    let part = parts::build_pocketed_tray().expect("pocketed tray");
    let mut brep = part.brep().clone();
    let arcs: Vec<_> = brep
        .edges
        .iter()
        .filter(|(_, e)| {
            e.curve
                .as_ref()
                .is_some_and(|c| matches!(c.geometry, CurveGeometry::Arc { .. }))
        })
        .map(|(h, _)| h)
        .collect();
    assert!(arcs.len() >= 8, "{} arcs", arcs.len());
    assert!(
        brep.validate_curve_window_all().is_empty(),
        "{:?}",
        brep.validate_curve_window_all()
    );
    // Negative control: move one arc's window onto its complement.
    let curve = brep.edges[arcs[0]].curve.as_mut().expect("curve");
    let (t0, t1) = (curve.t_start, curve.t_end);
    (curve.t_start, curve.t_end) = (t1, t0 + std::f64::consts::TAU);
    let flagged = brep.validate_curve_window_all();
    assert!(
        matches!(
            flagged[..],
            [ValidationError::CurveWindowOutsideSpan { edge, .. }] if edge == arcs[0]
        ),
        "{flagged:?}"
    );
    assert!(
        brep.validate_all()
            .iter()
            .any(|e| matches!(e, ValidationError::CurveWindowOutsideSpan { .. }))
    );
}
