use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::*;
use crate::error::BoolOperation;

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(Tolerance::default())
        .expect("cube builds")
        .transform(&Transform3::translation(x, y, z))
        .expect("translates")
}

#[test]
fn separated_and_face_touching_cubes_have_disjoint_interiors() {
    assert!(verify_disjoint_interiors(
        &cube_at(0.0, 0.0, 0.0),
        &cube_at(3.0, 0.0, 0.0)
    ));
    assert!(verify_disjoint_interiors(
        &cube_at(0.0, 0.0, 0.0),
        &cube_at(1.0, 0.0, 0.0)
    ));
    assert!(verify_disjoint_interiors(
        &cube_at(0.0, 0.0, 0.0),
        &cube_at(1.0, 1.0, 1.0)
    ));
}

/// Negative controls: overlapping, nested and identical operands are
/// never certified empty.
#[test]
fn overlapping_nested_and_identical_operands_are_not_disjoint() {
    let unit = cube_at(0.0, 0.0, 0.0);
    assert!(!verify_disjoint_interiors(&unit, &cube_at(0.5, 0.5, 0.5)));
    let small = BRep::unit_cube(Tolerance::default())
        .expect("cube builds")
        .transform(&Transform3::scale_uniform(0.25))
        .expect("scales")
        .transform(&Transform3::translation(0.4, 0.4, 0.4))
        .expect("translates");
    assert!(!verify_disjoint_interiors(&unit, &small));
    assert!(!verify_disjoint_interiors(&small, &unit));
    assert!(!verify_disjoint_interiors(&unit, &cube_at(0.0, 0.0, 0.0)));
}

/// A point on a contact line is on the boundary even where the ray
/// vote would call it inside.
#[test]
fn a_line_contact_point_is_on_the_rod_boundary() {
    let tol = Tolerance::default();
    let rod = BRep::cylinder(0.5, 1.0, tol)
        .expect("cylinder builds")
        .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
        .expect("rotates")
        .transform(&Transform3::translation(0.0, 0.5, 1.5))
        .expect("translates");
    let aabbs = crate::intersect::compute_face_aabbs(&rod);
    let on = popkern_core::geometry::Point3::new(0.5, 0.5, 1.0);
    let off = popkern_core::geometry::Point3::new(0.5, 0.5, 0.99);
    use crate::postprocess::zero_thickness::on_body_boundary;
    assert!(on_body_boundary(&on, &rod, &aabbs, &tol));
    assert!(!on_body_boundary(&off, &rod, &aabbs, &tol));
}

#[test]
fn an_empty_union_is_never_certified() {
    let r = certify_empty_result(
        BoolOperation::Union,
        &cube_at(0.0, 0.0, 0.0),
        &cube_at(3.0, 0.0, 0.0),
    );
    assert!(r.is_err());
}

#[test]
fn an_empty_subtraction_needs_a_inside_b() {
    let unit = cube_at(0.0, 0.0, 0.0);
    let big = BRep::unit_cube(Tolerance::default())
        .expect("cube builds")
        .transform(&Transform3::scale_uniform(3.0))
        .expect("scales")
        .transform(&Transform3::translation(-1.0, -1.0, -1.0))
        .expect("translates");
    assert!(certify_empty_result(BoolOperation::Subtraction, &unit, &big).is_ok());
    assert!(certify_empty_result(BoolOperation::Subtraction, &big, &unit).is_err());
}

/// A body coincident with the other lies within it: its corners sit
/// on the other's boundary, where the ray vote has no verdict.
#[test]
fn a_coincident_body_is_certified_inside() {
    let unit = cube_at(0.0, 0.0, 0.0);
    assert!(verify_whole_body_containment(
        &unit,
        &cube_at(0.0, 0.0, 0.0)
    ));
}

/// Negative control: every vertex and boundary sample of the cube lies
/// on the notched cube's boundary, but its top face runs over the notch.
#[test]
fn a_face_over_a_notch_is_not_certified_inside() {
    let tol = Tolerance::default();
    let notch = BRep::unit_cube(tol)
        .expect("cube builds")
        .transform(&Transform3::scale(0.4, 0.4, 0.4))
        .expect("scales")
        .transform(&Transform3::translation(0.3, 0.3, 0.8))
        .expect("translates");
    let notched = crate::boolean(
        BoolOperation::Subtraction,
        cube_at(0.0, 0.0, 0.0),
        notch,
        tol,
        None,
    )
    .expect("notch cut");
    assert!(!verify_whole_body_containment(
        &cube_at(0.0, 0.0, 0.0),
        &notched
    ));
}

/// Cubes sharing a face share area; cubes meeting along an edge or at
/// a corner, and a ball resting on a cube, do not (PK-1124).
#[test]
fn only_face_sharing_operands_share_area() {
    let tol = Tolerance::default();
    let unit = cube_at(0.0, 0.0, 0.0);
    assert!(shares_area(&unit, &cube_at(1.0, 0.0, 0.0), &tol));
    assert!(shares_area(&unit, &cube_at(1.0, 0.5, 0.0), &tol));
    assert!(!shares_area(&unit, &cube_at(1.0, 1.0, 0.0), &tol));
    assert!(!shares_area(&unit, &cube_at(1.0, 1.0, 1.0), &tol));
    let ball = BRep::sphere(0.5, tol)
        .expect("sphere builds")
        .transform(&Transform3::translation(0.5, 0.5, 1.5))
        .expect("translates");
    assert!(!shares_area(&unit, &ball, &tol));
    assert!(!shares_area(&ball, &unit, &tol));
}
