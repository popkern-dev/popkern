//! face_area of trimmed cylinder, sphere and cone faces.
//!
//! PK-934: a face without inner loops read the closed form of its
//! pcurves' uv BOUNDING rectangle: the notched cone wall of a cone minus
//! a corner box read 8.781 for a lateral area of at most 7.025, the
//! octant-cut sphere half of its true 7/8. The closed form now applies
//! only when the outer loop bounds the rectangle.
//!
//! PK-935: any other trim, holes included, integrates its area exactly
//! over its boundary (it read the measuring mesh, 1.4e-4 relative off),
//! in every orientation.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::face_area;

const PI: f64 = std::f64::consts::PI;

fn curved_area(body: &BRep, tol: &Tolerance, pick: fn(&SurfaceGeometry) -> bool) -> f64 {
    body.faces
        .iter()
        .filter(|(_, f)| f.surface.as_ref().is_some_and(|s| pick(&s.geometry)))
        .map(|(h, _)| face_area(body, h, tol).expect("area"))
        .sum()
}

fn is_cone(g: &SurfaceGeometry) -> bool {
    matches!(g, SurfaceGeometry::Cone { .. })
}
fn is_sphere(g: &SurfaceGeometry) -> bool {
    matches!(g, SurfaceGeometry::Sphere { .. })
}
fn is_cylinder(g: &SurfaceGeometry) -> bool {
    matches!(g, SurfaceGeometry::Cylinder { .. })
}

fn rotations() -> Vec<(&'static str, Transform3)> {
    vec![
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
        (
            "mixed",
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
    ]
}

fn op(kind: BoolOperation, a: BRep, b: BRep, r: &Transform3, tol: Tolerance) -> BRep {
    let a = a.transform(r).expect("rotate");
    let b = b.transform(r).expect("rotate");
    boolean(kind, a, b, tol, None).expect("boolean")
}

fn assert_close(what: &str, got: f64, want: f64) {
    assert!((got - want).abs() < 1e-10 * want, "{what}: {got} vs {want}");
}

#[test]
fn notched_and_octant_faces_measure_their_region_exactly() {
    let tol = Tolerance::default();
    for (name, r) in rotations() {
        // Cone (radius 1, height 2) less the corner box [0,1]^2 x
        // [0.8,1.8]: the lateral area less a quarter of the frustum band
        // between radii 0.6 and 0.1 (slant √1.25).
        let cube = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 0.8)))
            .expect("box");
        let cone = BRep::cone(1.0, 2.0, tol).expect("cone");
        let notched = op(BoolOperation::Subtraction, cone, cube, &r, tol);
        let want = PI * 5.0_f64.sqrt() - PI * 0.7 * 1.25_f64.sqrt() / 4.0;
        assert_close(
            &format!("{name}: cone wall"),
            curved_area(&notched, &tol, is_cone),
            want,
        );

        // A box and a sphere on its corner: 7/8 of the sphere shows on
        // the union, the inward octant on the difference.
        let ball = || {
            BRep::sphere(0.5, tol)
                .and_then(|s| s.transform(&Transform3::translation(1.0, 1.0, 1.0)))
                .expect("sphere")
        };
        let cube = || BRep::unit_cube(tol).expect("box");
        let united = op(BoolOperation::Union, cube(), ball(), &r, tol);
        assert_close(
            &format!("{name}: 7/8 sphere"),
            curved_area(&united, &tol, is_sphere),
            7.0 * PI / 8.0,
        );
        let cut = op(BoolOperation::Subtraction, cube(), ball(), &r, tol);
        assert_close(
            &format!("{name}: inward octant"),
            curved_area(&cut, &tol, is_sphere),
            PI / 8.0,
        );
    }
}

#[test]
fn holed_cylinder_and_sphere_faces_measure_exactly() {
    let tol = Tolerance::default();
    for (name, r) in rotations() {
        // A box pocketed into a cylinder wall (radius 0.5) at u = π/2:
        // the hole spans |x| < 0.1, i.e. 2·asin(0.2) of azimuth, over
        // height 0.2.
        let cylinder = BRep::cylinder(0.5, 2.0, tol).expect("cylinder");
        let pocket = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::scale_uniform(0.2)))
            .and_then(|c| c.transform(&Transform3::translation(-0.1, 0.4, 0.5)))
            .expect("pocket");
        let holed = op(BoolOperation::Subtraction, cylinder, pocket, &r, tol);
        let want = 2.0 * PI * 0.5 * 2.0 - 0.5 * 2.0 * 0.2_f64.asin() * 0.2;
        assert_close(
            &format!("{name}: holed wall"),
            curved_area(&holed, &tol, is_cylinder),
            want,
        );

        // A sphere (radius 0.5) with a cap cut off by the plane x = 0.3:
        // its hole is a circle off the poles, and the cap has height 0.2.
        let ball = BRep::sphere(0.5, tol).expect("sphere");
        let slab = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::translation(0.3, -0.5, -0.5)))
            .expect("slab");
        let capped = op(BoolOperation::Subtraction, ball, slab, &r, tol);
        let want = PI - 2.0 * PI * 0.5 * 0.2;
        assert_close(
            &format!("{name}: capped sphere"),
            curved_area(&capped, &tol, is_sphere),
            want,
        );
    }
}

/// A sphere (radius 0.5) drilled along x by a radius-0.2 bore: two caps
/// of height 0.5 − √0.21 leave it. The boolean's hole arcs sit up to
/// 7e-8 off both surfaces (marched, then promoted to arcs; PK-937), so
/// the area is exact to that trim, 2e-7 relative.
#[test]
fn a_drilled_sphere_measures_to_its_trim() {
    let tol = Tolerance::default();
    let ball = BRep::sphere(0.5, tol).expect("sphere");
    let bore = BRep::cylinder(0.2, 2.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.0)))
        .and_then(|c| c.transform(&Transform3::rotation_y(PI / 2.0)))
        .expect("bore");
    let drilled = boolean(BoolOperation::Subtraction, ball, bore, tol, None).expect("drill");
    let want = PI - 2.0 * 2.0 * PI * 0.5 * (0.5 - 0.21_f64.sqrt());
    let got = curved_area(&drilled, &tol, is_sphere);
    assert!(
        (got - want).abs() < 1e-6 * want,
        "drilled sphere {got} vs {want}"
    );
}

#[test]
fn a_rectangle_bounded_face_keeps_its_closed_form() {
    let tol = Tolerance::default();
    let cylinder = BRep::cylinder(0.5, 2.0, tol).expect("cylinder");
    let got = curved_area(&cylinder, &tol, is_cylinder);
    let want = 2.0 * PI * 0.5 * 2.0;
    assert!(
        (got - want).abs() < 1e-12 * want,
        "cylinder wall {got} vs {want}"
    );
}
