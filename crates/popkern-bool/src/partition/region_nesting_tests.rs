//! Regions bounded by several shells: bodies floating inside others.
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::super::{BodySet, SetExpr, partition};

/// An axis-aligned cube of side `side` with its low corner at `at`.
fn cube(side: f64, at: f64) -> BRep {
    BRep::unit_cube(Tolerance::default())
        .and_then(|c| c.transform(&Transform3::scale_uniform(side)))
        .and_then(|c| c.transform(&Transform3::translation(at, at, at)))
        .expect("cube")
}

/// (faces, shells) of each region solid, sorted.
fn solids(bodies: &[&BRep]) -> Vec<(usize, usize, BodySet)> {
    let p = partition(bodies, Tolerance::default()).expect("partition");
    let cells = p.region_solids(&SetExpr::Any).expect("region solids");
    let mut out: Vec<_> = cells
        .solids
        .iter()
        .map(|s| {
            let findings = s.body.validate_all();
            assert!(findings.is_empty(), "{findings:?}");
            (s.body.faces.len(), s.body.back_shell_count(), s.bodies)
        })
        .collect();
    out.sort_unstable_by_key(|s| (s.0, s.1));
    out
}

#[test]
fn a_cube_floating_in_a_cube_leaves_it_hollow() {
    let (big, small) = (cube(1.0, 0.0), cube(0.5, 0.25));
    assert_eq!(
        solids(&[&big, &small]),
        vec![(6, 1, BodySet::of(&[0, 1])), (12, 2, BodySet::of(&[0]))]
    );
    let cut = crate::split(&big, &[&small], Tolerance::default()).expect("split");
    let mut faces: Vec<_> = cut.pieces.iter().map(|p| p.body.faces.len()).collect();
    faces.sort_unstable();
    assert_eq!(faces, vec![6, 12]);
}

#[test]
fn three_nested_cubes_give_two_hollow_shells_and_a_core() {
    let (a, b, c) = (cube(1.0, 0.0), cube(0.6, 0.2), cube(0.2, 0.4));
    assert_eq!(
        solids(&[&a, &b, &c]),
        vec![
            (6, 1, BodySet::of(&[0, 1, 2])),
            (12, 2, BodySet::of(&[0])),
            (12, 2, BodySet::of(&[0, 1])),
        ]
    );
}

#[test]
fn a_hollow_operand_is_one_region() {
    let tol = Tolerance::default();
    let hollow = crate::boolean(
        crate::BoolOperation::Subtraction,
        cube(1.0, 0.0),
        cube(0.5, 0.25),
        tol,
        None,
    )
    .expect("hollow");
    assert_eq!(hollow.back_shell_count(), 2);
    assert_eq!(solids(&[&hollow]), vec![(12, 2, BodySet::of(&[0]))]);
    // A cube in the void is a region of its own, the void unkept.
    let core = cube(0.2, 0.4);
    assert_eq!(
        solids(&[&hollow, &core]),
        vec![(6, 1, BodySet::of(&[1])), (12, 2, BodySet::of(&[0]))]
    );
}

#[test]
fn apart_cubes_stay_two_regions() {
    let (a, b) = (cube(1.0, 0.0), cube(1.0, 3.0));
    assert_eq!(
        solids(&[&a, &b]),
        vec![(6, 1, BodySet::of(&[0])), (6, 1, BodySet::of(&[1]))]
    );
}
