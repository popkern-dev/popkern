use super::*;
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

/// `body`'s curved faces' boundary-integral areas, summed.
fn curved_area(body: &BRep) -> f64 {
    body.faces
        .keys()
        .filter_map(|f| boundary_integral_area(body, f))
        .sum()
}

#[test]
fn primitive_walls_match_their_closed_forms() {
    let tol = Tolerance::default();
    let tilt = Transform3::rotation_x(1.1).compose(&Transform3::rotation_z(0.7));
    let cases = [
        (
            "cylinder",
            BRep::cylinder(0.7, 1.3, tol),
            2.0 * PI * 0.7 * 1.3,
        ),
        ("cone", BRep::cone(1.0, 2.0, tol), PI * 5.0_f64.sqrt()),
        ("sphere", BRep::sphere(0.5, tol), PI),
    ];
    for (name, body, want) in cases {
        let body = body.expect("primitive");
        for (frame, b) in [
            ("home", body.clone()),
            ("tilted", body.transform(&tilt).expect("tilt")),
        ] {
            let got = curved_area(&b);
            assert!(
                (got - want).abs() < 1e-12 * want,
                "{name} {frame}: {got} vs {want}"
            );
        }
    }
}

#[test]
fn a_reversed_face_flag_reads_the_same_area() {
    let tol = Tolerance::default();
    let mut ball = BRep::sphere(0.5, tol).expect("sphere");
    // Turning the body inside out flips every flag and loop together.
    ball.reverse_body_orientation();
    let got = curved_area(&ball);
    assert!((got - PI).abs() < 1e-12 * PI, "{got}");
}

#[test]
fn planes_have_no_boundary_integral_here() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    assert!(
        cube.faces
            .keys()
            .all(|f| boundary_integral_area(&cube, f).is_none())
    );
}
