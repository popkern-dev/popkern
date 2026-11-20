//! A ball resting on a box's top face where the operands otherwise
//! cross (PK-1155's one-shell touch): the face-pair pass finds the touch,
//! the table resolves it on no record, and the result's contact is that
//! candidate's (tangency design §4.8, W3b).

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{BoolError, BoolOperation, boolean};

/// The box [0, 10]³ and a C prism (a plate under the box's foot, a
/// column, a roof) holding a ball of radius 3 at (6, 5, 13) resting on
/// the box's top at (6, 5, 10); `turned` turns the ball about its
/// centre, so no vertex of it lies near the touch.
fn box_and_holder(turned: bool) -> (BRep, BRep) {
    let tol = Tolerance::default();
    let rect: Vec<Point3> = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
        .iter()
        .map(|&(x, y)| Point3::new(x, y, 0.0))
        .collect();
    let block = BRep::prism(&rect, Vec3::new(0.0, 0.0, 1.0), 10.0, tol).expect("box");
    let profile: Vec<Point3> = [
        (8.0, 0.0),
        (20.0, 0.0),
        (20.0, 20.0),
        (2.0, 20.0),
        (2.0, 15.0),
        (12.0, 15.0),
        (12.0, 2.0),
        (8.0, 2.0),
    ]
    .iter()
    .map(|&(x, z)| Point3::new(x, 0.0, z))
    .collect();
    let c = BRep::prism(&profile, Vec3::new(0.0, 1.0, 0.0), 10.0, tol).expect("prism");
    let centre = Point3::new(6.0, 5.0, 13.0);
    let mut ball = BRep::sphere(3.0, tol)
        .and_then(|s| s.transform(&Transform3::translation(6.0, 5.0, 13.0)))
        .expect("ball");
    if turned {
        let turn =
            Transform3::rotation_around_axis(centre, Vec3::new(1.0, 2.0, 3.0), 0.7).expect("axis");
        ball = ball.transform(&turn).expect("turned");
    }
    let holder = boolean(BoolOperation::Union, c, ball, tol, None).expect("the roof holds it");
    (block, holder)
}

const TOUCH: Point3 = Point3 {
    x: 6.0,
    y: 5.0,
    z: 10.0,
};

#[test]
fn a_touch_on_no_record_is_a_contact_vertex_of_the_result() {
    for turned in [false, true] {
        let (block, holder) = box_and_holder(turned);
        let tol = Tolerance::default();
        let ix = crate::intersect::find_intersections(&block, &holder, &tol).expect("intersect");
        assert_eq!(ix.contacts.len(), 1, "turned {turned}: {:?}", ix.contacts);
        assert!(ix.contacts[0].point.distance_to(&TOUCH) < 1e-6);
        let hits = super::super::candidates::interfering_hits(
            (&block, &holder),
            &ix,
            &ix.records,
            tol.linear(),
        );
        let table = super::super::table::build(
            (&block, &holder),
            (&hits, &ix.contacts),
            &ix.records,
            tol.linear(),
        );
        let isolated = table.isolated_contacts((&block, &holder), tol.linear());
        assert_eq!(isolated.len(), 1, "turned {turned}");
        // The decline names the candidate's vertex.
        match boolean(BoolOperation::Union, block, holder, tol, None) {
            Err(BoolError::NonManifoldResult { contacts }) => {
                assert!(
                    contacts
                        .iter()
                        .any(|c| c.location.distance_to(&isolated[0]) <= tol.linear_reach()),
                    "turned {turned}: {contacts:?} miss {:?}",
                    isolated[0]
                );
            }
            other => panic!("turned {turned}: {other:?}"),
        }
    }
}
