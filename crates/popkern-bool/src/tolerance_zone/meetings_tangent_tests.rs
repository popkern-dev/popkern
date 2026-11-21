//! A cylinder and a ball tangent at one point with relative curvature of
//! both signs (a V configuration): their intersection is a figure-eight
//! crossing itself at the tangency.

use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The cylinder of radius 1 on the z axis and the ball of radius 1.5
/// centred at (-0.5, 0, 0): both pass (1, 0, 0) with normal +x; across z
/// the cylinder is flat and the ball bends away, around the axis the
/// cylinder bends more: the faces cross both ways there.
#[test]
fn a_tangent_crossing_is_one_vertex_at_the_contact() {
    let tol = Tolerance::default();
    let cylinder = BRep::cylinder(1.0, 4.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -2.0)))
        .expect("cylinder");
    let ball = BRep::sphere(1.5, tol)
        .and_then(|s| s.transform(&Transform3::rotation_z(0.9)))
        .and_then(|s| s.transform(&Transform3::translation(-0.5, 0.0, 0.0)))
        .expect("ball");
    let ix = crate::intersect::find_intersections(&cylinder, &ball, &tol).expect("intersect");
    let hits = super::super::candidates::interfering_hits(
        (&cylinder, &ball),
        &ix,
        &ix.records,
        tol.linear(),
    );
    let table = super::super::table::build(
        (&cylinder, &ball),
        (&hits, &ix.contacts),
        &ix.records,
        tol.linear(),
    );
    let p = Point3::new(1.0, 0.0, 0.0);
    // Every candidate at the crossing (the records' ends there, or their
    // meetings refined onto it) is one vertex, standing at the point.
    let res = table.resolution.as_ref().expect("resolved");
    let near: Vec<usize> = (0..table.cands.len())
        .filter(|&i| table.cands[i].position.distance_to(&p) < 1e-3)
        .collect();
    assert!(!near.is_empty());
    let group = res.group_of[near[0]];
    assert!(near.iter().all(|&i| res.group_of[i] == group), "{near:?}");
    let anchor = table.cands[res.groups[group].anchor].position;
    assert!(anchor.distance_to(&p) < 1e-9, "{anchor:?}");
}
