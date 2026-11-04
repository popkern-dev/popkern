//! The table closes each closed record at a vertex the body has anyway.

use popkern_core::geometry::Transform3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::intersect::IntersectionRecord;

/// A ball drilled across, both turned off their seams, and the
/// boolean's records before the table decides them.
fn drilled_ball() -> (BRep, BRep, crate::intersect::IntersectionResult) {
    let tol = Tolerance::default();
    let turn = Transform3::rotation_z(0.7);
    let ball = BRep::sphere(0.5, tol)
        .and_then(|s| s.transform(&turn))
        .expect("sphere");
    let bore = BRep::cylinder(0.2, 2.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.0)))
        .and_then(|c| c.transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)))
        .and_then(|c| c.transform(&turn))
        .expect("bore");
    let ix = crate::intersect::find_intersections(&ball, &bore, &tol).expect("intersections");
    (ball, bore, ix)
}

fn closure(r: &IntersectionRecord) -> popkern_core::geometry::Point3 {
    let c = &r.curve.curve_3d;
    c.point_at(c.parameter_range().0)
}

/// Every closed rim closes where an interfering hit lies (the seams
/// crossing it), and the table over the moved records moves nothing
/// more: no closure is left where only its record's own ends are.
#[test]
fn every_closed_rim_closes_at_a_hit_on_it() {
    let (ball, bore, ix) = drilled_ball();
    let tol = Tolerance::default();
    let (records, identity) = super::decide((&ball, &bore), &ix, ix.records.clone(), &tol, false);
    assert!(!identity.graph_identity());
    let closed: Vec<&IntersectionRecord> = records.iter().filter(|r| r.curve.is_closed).collect();
    assert!(!closed.is_empty());
    let moved = records
        .iter()
        .zip(&ix.records)
        .filter(|(r, o)| closure(r).distance_to(&closure(o)) > 1e-9)
        .count();
    assert!(moved > 0, "the tracer's closures already sat at the hits");
    for r in &closed {
        let p = closure(r);
        assert!(
            identity.hits.iter().any(|h| h.point.distance_to(&p) < 1e-9),
            "a rim closes at {p:?}, at no hit"
        );
    }
    let hits = crate::tolerance_zone::candidates::interfering_hits(
        (&ball, &bore),
        &ix,
        &records,
        tol.linear(),
    );
    let table = crate::tolerance_zone::table::build((&ball, &bore), &hits, &records, tol.linear());
    assert!(table.closure_moves(&records).iter().all(Option::is_none));
}
