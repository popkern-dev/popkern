use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::shell_contacts;
use crate::operation::helpers::combine_breps;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn ball(r: f64, c: (f64, f64, f64)) -> BRep {
    BRep::sphere(r, tol())
        .expect("sphere builds")
        .transform(&Transform3::translation(c.0, c.1, c.2))
        .expect("translates")
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .expect("cube builds")
        .transform(&Transform3::translation(x, y, z))
        .expect("translates")
}

#[test]
fn a_ball_resting_on_a_face_is_a_contact() {
    let body = combine_breps(cube_at(0.0, 0.0, 0.0), ball(0.4, (0.5, 0.5, 1.4))).expect("combines");
    let contacts = shell_contacts(&body, &tol());
    assert_eq!(contacts.len(), 1, "{contacts:?}");
    let at = contacts[0].location;
    assert!(
        at.distance_to(&popkern_core::geometry::Point3::new(0.5, 0.5, 1.0)) < 1e-5,
        "{at:?}"
    );
}

/// Negative controls: a ball one hundred tolerances clear of the face,
/// a ball beside the cube with overlapping boxes but no contact, and a
/// cavity shell with a wall.
#[test]
fn clear_shells_are_not_contacts() {
    let lifted =
        combine_breps(cube_at(0.0, 0.0, 0.0), ball(0.4, (0.5, 0.5, 1.4001))).expect("combines");
    assert!(shell_contacts(&lifted, &tol()).is_empty());
    // The ball's box overlaps the cube's corner region; the sphere
    // clears the corner (distance from (1,1,1) to the centre is 0.52).
    let beside =
        combine_breps(cube_at(0.0, 0.0, 0.0), ball(0.4, (1.3, 1.3, 1.3))).expect("combines");
    assert!(shell_contacts(&beside, &tol()).is_empty());
    let cavity = crate::operation::helpers::copy_brep_into(
        BRep::unit_cube(tol())
            .expect("cube")
            .transform(&Transform3::scale_uniform(3.0))
            .expect("scales"),
        ball(0.5, (1.5, 1.5, 1.5)),
        true,
    )
    .expect("combines");
    assert!(shell_contacts(&cavity, &tol()).is_empty());
}

#[test]
fn a_cavity_touching_the_outer_wall_is_a_contact() {
    let cavity = crate::operation::helpers::copy_brep_into(
        BRep::unit_cube(tol())
            .expect("cube")
            .transform(&Transform3::scale_uniform(3.0))
            .expect("scales"),
        ball(0.5, (1.5, 1.5, 2.5)),
        true,
    )
    .expect("combines");
    assert_eq!(shell_contacts(&cavity, &tol()).len(), 1);
}

#[test]
fn one_shell_has_no_shell_contacts() {
    assert!(shell_contacts(&cube_at(0.0, 0.0, 0.0), &tol()).is_empty());
}

#[test]
fn contact_extent_is_point_line_or_area() {
    use super::{ContactExtent, shell_findings};
    let extents = |body: &BRep| -> Vec<ContactExtent> {
        shell_findings(body, &tol())
            .into_iter()
            .map(|f| f.extent)
            .collect()
    };
    let resting =
        combine_breps(cube_at(0.0, 0.0, 0.0), ball(0.4, (0.5, 0.5, 1.4))).expect("combines");
    assert_eq!(extents(&resting), vec![ContactExtent::Point]);
    let rod = BRep::cylinder(0.5, 1.0, tol())
        .expect("cylinder")
        .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
        .expect("rotates")
        .transform(&Transform3::translation(0.0, 0.5, 1.5))
        .expect("translates");
    let lying = combine_breps(cube_at(0.0, 0.0, 0.0), rod).expect("combines");
    assert!(
        extents(&lying).contains(&ContactExtent::Line),
        "{:?}",
        extents(&lying)
    );
    assert!(!extents(&lying).contains(&ContactExtent::Area));
    let flush = combine_breps(cube_at(0.0, 0.0, 0.0), cube_at(1.0, 0.0, 0.0)).expect("combines");
    assert!(
        extents(&flush).contains(&ContactExtent::Area),
        "{:?}",
        extents(&flush)
    );
}

/// An area touch is never served split: the bodies should have merged.
#[test]
fn an_area_touch_is_not_accepted_as_a_split() {
    let flush = combine_breps(cube_at(0.0, 0.0, 0.0), cube_at(1.0, 0.0, 0.0)).expect("combines");
    assert!(
        super::super::split_pieces::accept_split(
            flush,
            &tol(),
            crate::NonManifoldPolicy::SplitToManifold
        )
        .is_err()
    );
}

/// A ball resting on a cube's edge touches its two faces only at their
/// common boundary: the surfaces pass through each other beyond it, so
/// the touch is an edge against a face (R5.8). Lifted ten tolerances
/// off the edge it is clear.
#[test]
fn a_ball_resting_on_an_edge_is_a_contact() {
    let r = 0.4;
    let at = |gap: f64| {
        let d = (r + gap) / std::f64::consts::SQRT_2;
        combine_breps(cube_at(0.0, 0.0, 0.0), ball(r, (1.0 + d, 0.5, 1.0 + d))).expect("combines")
    };
    let contacts = shell_contacts(&at(0.0), &tol());
    assert_eq!(contacts.len(), 1, "{contacts:?}");
    let p = contacts[0].location;
    assert!(
        p.distance_to(&popkern_core::geometry::Point3::new(1.0, 0.5, 1.0)) < 1e-5,
        "{p:?}"
    );
    assert!(shell_contacts(&at(1e-5), &tol()).is_empty());
}
