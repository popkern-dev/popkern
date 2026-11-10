use super::first_self_contact;
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// A box [-10, 10]² × [0, 10] with a through hole of radius 3 along x
/// at height 5: the hole's lowest line runs 2 above the floor and its
/// highest 2 below the top.
fn holed_box(tol: Tolerance) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(-10.0, -10.0, z),
            Point3::new(10.0, -10.0, z),
            Point3::new(10.0, 10.0, z),
            Point3::new(-10.0, 10.0, z),
        ]
    };
    let bx = BRep::polygon_frustum(&rect(0.0), &rect(10.0), tol).expect("box");
    let rod = BRep::cylinder(3.0, 30.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -15.0)))
        .and_then(|c| c.transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)))
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 5.0)))
        .expect("rod");
    popkern_bool::boolean(popkern_bool::BoolOperation::Subtraction, bx, rod, tol, None)
        .expect("hole")
}

fn shelled(t: f64, tol: Tolerance) -> BRep {
    crate::shell_impl_unchecked(
        &holed_box(tol),
        t,
        crate::ShellDirection::Inward,
        &[],
        &[],
        &crate::ShellOptions::default(),
        tol,
        &mut popkern_core::progress::StageTracker::inert(),
        None,
    )
    .expect("shell")
    .brep
}

#[test]
fn a_body_that_touches_nowhere_has_no_self_contact() {
    let tol = Tolerance::default();
    assert_eq!(first_self_contact(&holed_box(tol), &tol), None);
    assert_eq!(first_self_contact(&shelled(0.9, tol), &tol), None);
}

#[test]
fn skins_meeting_where_a_wall_closes_are_a_self_contact() {
    // Shelled by 1 the hole's offset cylinder (radius 4) is tangent to
    // the floor's offset plane z = 1 along y = 0 (and to the top's,
    // z = 9).
    let tol = Tolerance::default();
    let contact = first_self_contact(&shelled(1.0, tol), &tol).expect("contact");
    let p = contact.location;
    assert!(p.y.abs() < 1e-3, "{p:?}");
    assert!(
        (p.z - 1.0).abs() < 1e-5 || (p.z - 9.0).abs() < 1e-5,
        "{p:?}"
    );
}
