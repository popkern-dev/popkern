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
