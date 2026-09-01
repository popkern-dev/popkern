//! PK-394 — union of a tool whose cap is coplanar with a host face that
//! already carries an inner loop serves a valid body in either order.
//!
//! The boss's bottom rim lies ON the plate top; its closed record was
//! dropped by the UV clip (all-outside verdicts on the tool lateral
//! whose on-boundary midpoints failed a sample-sparse, unfolded
//! clearance test), so the boss's cap survived as an interior face
//! (pillar-first) or the assembly violated Euler (boss-first). The
//! flush rim is now kept as the imprint record.

use popkern_part::prelude::*;

fn parts() -> (Part, Part, Part) {
    let plate = Profile::on(Plane::XY)
        .rect(2.4, 1.6)
        .extrude(0.4)
        .expect("plate");
    let pillar = Profile::on(Plane::xy_at(0.4))
        .rect(0.4, 0.4)
        .at(0.9, 0.5)
        .extrude(1.0)
        .expect("pillar");
    let boss = Profile::on(Plane::xy_at(0.4))
        .circle(0.35)
        .extrude(1.0)
        .expect("boss");
    (plate, pillar, boss)
}

fn check(label: &str, part: &Part) {
    assert_eq!(
        part.face_count(),
        13,
        "{label}: plate 6 + pillar 5 + boss wall + boss top"
    );
    let issues = part.brep().validate_all();
    assert!(issues.is_empty(), "{label}: {issues:?}");
    assert_eq!(
        part.open_edge_count().expect("mesh"),
        0,
        "{label}: watertight"
    );
}

#[test]
fn three_solid_coplanar_union_serves_in_either_order() {
    let (plate, pillar, boss) = parts();
    let pillar_first = plate
        .union(&pillar)
        .expect("plate ∪ pillar")
        .union(&boss)
        .expect("∪ boss");
    check("pillar-first", &pillar_first);
    let boss_first = plate
        .union(&boss)
        .expect("plate ∪ boss")
        .union(&pillar)
        .expect("∪ pillar");
    check("boss-first", &boss_first);
}
