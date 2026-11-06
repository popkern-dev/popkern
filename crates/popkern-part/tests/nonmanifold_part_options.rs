//! Non-manifold policy on the Part surface (PK-964): parts touching
//! along an edge decline to union by default and combine into
//! manifold pieces with recorded contacts on request; the symmetric
//! difference of overlapping parts serves its two differences.

use popkern_bool::{BoolError, BoolOperation, BoolOptions, NonManifoldPolicy};
use popkern_part::prelude::*;

fn cube() -> Part {
    Profile::on(Plane::xy_at(0.0))
        .rect(1.0, 1.0)
        .extrude(1.0)
        .expect("unit cube")
}

const SPLIT: BoolOptions = BoolOptions {
    non_manifold: NonManifoldPolicy::SplitToManifold,
};

#[test]
fn edge_touching_parts_decline_then_split_on_request() {
    let other = cube().translate(1.0, 1.0, 0.0).expect("translate");
    let declined = cube().union(&other);
    assert!(
        matches!(
            declined,
            Err(PartError::Boolean {
                source: BoolError::NonManifoldResult { .. },
                ..
            })
        ),
        "{declined:?}"
    );
    let split = cube()
        .boolean_with(&other, BoolOperation::Union, &SPLIT)
        .expect("splits");
    assert!(!split.brep().topology_contacts.is_empty());
    assert!(split.brep().validate_all().is_empty());
}

#[test]
fn symmetric_difference_of_overlapping_parts_serves_split() {
    let other = cube().translate(0.5, 0.5, 0.5).expect("translate");
    let xor = cube()
        .symmetric_difference(&other, &SPLIT)
        .expect("serves split");
    assert!(xor.brep().validate_all().is_empty());
    assert_eq!(xor.brep().back_shell_count(), 2);
    assert!(
        cube()
            .symmetric_difference(&other, &BoolOptions::default())
            .is_err()
    );
}
