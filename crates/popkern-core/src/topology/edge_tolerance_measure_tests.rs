//! Measured edge tolerance pins (PK-826).

use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn exact_bodies_keep_their_default_tolerances() {
    for mut b in [
        BRep::unit_cube(Tolerance::default()).unwrap(),
        BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap(),
        BRep::cone(1.0, 2.0, Tolerance::default()).unwrap(),
    ] {
        assert!(b.record_measured_edge_tolerances().is_empty());
        assert!(b.edges.values().all(|e| e.tolerance.is_none()));
    }
}

#[test]
fn a_misfit_pcurve_raises_its_edge_to_the_measured_gap() {
    let mut b = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // Nudge one rim use's pcurve 2e-5 in v (under the model's 1e-4
    // max_imprecision): its image leaves the rim by exactly that much
    // along the cylinder's axis.
    let (he, edge) = b
        .half_edges
        .iter()
        .find(|(_, h)| h.origin == b.half_edges[h.twin].origin)
        .map(|(he, h)| (he, h.edge))
        .unwrap();
    let pc = b.half_edges[he].pcurve.clone().unwrap();
    b.half_edges[he].pcurve = Some(pc.translated_uv(0.0, 2.0e-5));
    let raised = b.record_measured_edge_tolerances();
    assert_eq!(raised, vec![edge]);
    let tol = b.edges[edge].tolerance.unwrap();
    assert!((2.0e-5..4.0e-5).contains(&tol), "recorded {tol}");
    assert!(b.validate_pcurve_consistency_all().is_empty());
}

#[test]
fn a_misfit_past_the_ceiling_is_recorded_as_measured() {
    // A 3e-4 nudge, past the 1e-4 max_imprecision: recorded as 3e-4,
    // not hidden at the ceiling (PK-885 removed that containment).
    let mut b = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let (he, edge) = b
        .half_edges
        .iter()
        .find(|(_, h)| h.origin == b.half_edges[h.twin].origin)
        .map(|(he, h)| (he, h.edge))
        .unwrap();
    let pc = b.half_edges[he].pcurve.clone().unwrap();
    b.half_edges[he].pcurve = Some(pc.translated_uv(0.0, 3.0e-4));
    b.record_measured_edge_tolerances();
    let tol = b.edges[edge].tolerance.unwrap();
    assert!((3.0e-4..3.1e-4).contains(&tol), "recorded {tol}");
    assert!(tol > b.max_imprecision);
}
