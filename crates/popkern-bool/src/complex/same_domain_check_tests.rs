//! A unit cube lifted alone stands for a join whose partner copies
//! were dropped: its x = 1 face, coincident with a second cube's x = 0
//! face, is carried by one owner where the complex should hold both.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::super::ids::FaceId;
use super::super::model::FaceOwner;
use super::super::tests::complex_of_brep;
use super::same_domain_mismatches;

/// The face of `brep` whose vertices all have x = `x`.
fn face_at_x(brep: &BRep, x: f64) -> FaceHandle {
    brep.faces
        .keys()
        .find(|&f| {
            brep.face_outer_loop_positions(f)
                .iter()
                .all(|p| (p.x - x).abs() < 1e-9)
        })
        .expect("a face at x")
}

fn cubes(dy: f64) -> (BRep, BRep, FaceHandle, FaceHandle) {
    let tol = Tolerance::default();
    let a = BRep::unit_cube(tol).expect("cube");
    let b = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&Transform3::translation(1.0, dy, 0.0))
        .expect("placed");
    let (fa, fb) = (face_at_x(&a, 1.0), face_at_x(&b, 1.0));
    (a, b, fa, fb)
}

fn lifted_face(c: &super::super::model::CellComplex, face: FaceHandle) -> usize {
    c.faces
        .iter()
        .position(|f| f.owners[0].face == face)
        .expect("lifted")
}

#[test]
fn a_coincident_face_carried_by_one_operand_is_a_mismatch() {
    let (a, b, fa, fb) = cubes(0.0);
    let c = complex_of_brep(&a);
    let found = same_domain_mismatches(&c, &[&a, &b], &[(fa, fb)], 1e-6);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].face, FaceId::from_index(lifted_face(&c, fa)));
    assert_eq!(found[0].partner, fb);
}

#[test]
fn a_face_carrying_both_copies_is_not_a_mismatch() {
    let (a, b, fa, fb) = cubes(0.0);
    let mut c = complex_of_brep(&a);
    let k = lifted_face(&c, fa);
    c.faces[k].owners.push(FaceOwner {
        operand: 1,
        face: fb,
        same_sense: false,
    });
    assert!(same_domain_mismatches(&c, &[&a, &b], &[(fa, fb)], 1e-6).is_empty());
}

#[test]
fn a_face_outside_its_partner_region_is_not_a_mismatch() {
    // The partner's plane holds the face, its region does not: the
    // part of a partially overlapping face beyond the patch.
    let (a, b, fa, fb) = cubes(2.0);
    let c = complex_of_brep(&a);
    assert!(same_domain_mismatches(&c, &[&a, &b], &[(fa, fb)], 1e-6).is_empty());
}

#[test]
fn no_overlaps_no_mismatches() {
    let (a, b, ..) = cubes(0.0);
    let c = complex_of_brep(&a);
    assert!(same_domain_mismatches(&c, &[&a, &b], &[], 1e-6).is_empty());
}

/// A box 7.5e-4 thick in y beside a unit block, their tops coplanar:
/// the thin top's inner point, stepped 1e-2 in from a long side, would
/// cross onto the block's top. Every inner point stays on its face,
/// and the thin top (not the block's) is no mismatch.
#[test]
fn a_sliver_face_s_inner_point_stays_on_the_sliver() {
    let tol = Tolerance::default();
    let thin = 7.5e-4;
    let slab = |y0: f64, y1: f64| {
        let corners = [(0.0, y0), (1.0, y0), (1.0, y1), (0.0, y1)];
        let base = corners.map(|(x, y)| Point3::new(x, y, 0.0));
        BRep::prism(&base, Vec3::new(0.0, 0.0, 1.0), 1.0, tol).expect("slab")
    };
    let (a, b) = (slab(0.0, thin), slab(thin, 1.0 + thin));
    let top = |brep: &BRep| {
        brep.faces
            .keys()
            .find(|&f| {
                brep.face_outer_loop_positions(f)
                    .iter()
                    .all(|p| (p.z - 1.0).abs() < 1e-9)
            })
            .expect("top")
    };
    let c = complex_of_brep(&a);
    for k in 0..c.faces.len() {
        let at =
            super::inner_point(&c, &[&a, &b], 1e-6, FaceId::from_index(k)).expect("an inner point");
        let p = at.point;
        assert!(
            (-1e-9..=thin + 1e-9).contains(&p.y)
                && (-1e-9..=1.0 + 1e-9).contains(&p.x)
                && (-1e-9..=1.0 + 1e-9).contains(&p.z),
            "face {k}: {p:?}"
        );
    }
    let found = same_domain_mismatches(&c, &[&a, &b], &[(top(&a), top(&b))], 1e-6);
    assert!(found.is_empty(), "{found:?}");
}
