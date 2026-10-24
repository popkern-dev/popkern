//! PK-934: face_area of a cylinder, sphere, cone or torus face used the
//! closed form of its pcurves' uv BOUNDING rectangle for any face
//! without inner loops. A face whose region is not that rectangle read
//! the rectangle: the notched cone wall of a cone minus a corner box
//! read 8.781 for a lateral area of at most 7.025, the octant-cut
//! sphere half of its true 7/8. The closed form now applies only when
//! the outer loop bounds the rectangle (modulo the chart period);
//! anything else is tessellated.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::face_area;

fn curved_area(body: &BRep, tol: &Tolerance, pick: fn(&SurfaceGeometry) -> bool) -> f64 {
    body.faces
        .iter()
        .filter(|(_, f)| f.surface.as_ref().is_some_and(|s| pick(&s.geometry)))
        .map(|(h, _)| face_area(body, h, tol).expect("area"))
        .sum()
}

#[test]
fn notched_and_octant_faces_measure_their_region() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;

    // Cone (radius 1, height 2) less the corner box [0,1]^2 x [0.8,1.8]:
    // the full lateral area less a quarter of the frustum band between
    // radii 0.6 and 0.1 (slant 1.118).
    let cone = BRep::cone(1.0, 2.0, tol).expect("cone");
    let cube = BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 0.8)))
        .expect("box");
    let notched = boolean(BoolOperation::Subtraction, cone, cube, tol, None).expect("cut");
    let slant = (1.0_f64 + 0.25).sqrt();
    let want = pi * 5.0_f64.sqrt() - pi * (0.6 + 0.1) * slant / 4.0;
    let got = curved_area(&notched, &tol, |g| {
        matches!(g, SurfaceGeometry::Cone { .. })
    });
    assert!(
        (got - want).abs() < 1e-3 * want,
        "cone wall {got} vs {want}"
    );

    // A box united with a sphere on its corner: 7/8 of the sphere shows.
    let cube = BRep::unit_cube(tol).expect("box");
    let sphere = BRep::sphere(0.5, tol)
        .and_then(|s| s.transform(&Transform3::translation(1.0, 1.0, 1.0)))
        .expect("sphere");
    let united = boolean(BoolOperation::Union, cube, sphere, tol, None).expect("union");
    let want = 7.0 / 8.0 * 4.0 * pi * 0.25;
    let got = curved_area(&united, &tol, |g| {
        matches!(g, SurfaceGeometry::Sphere { .. })
    });
    assert!(
        (got - want).abs() < 1e-3 * want,
        "sphere face {got} vs {want}"
    );

    // A rectangle-bounded face keeps the exact closed form: a plain
    // cylinder's wall.
    let cylinder = BRep::cylinder(0.5, 2.0, tol).expect("cylinder");
    let got = curved_area(&cylinder, &tol, |g| {
        matches!(g, SurfaceGeometry::Cylinder { .. })
    });
    let want = 2.0 * pi * 0.5 * 2.0;
    assert!(
        (got - want).abs() < 1e-12 * want,
        "cylinder wall {got} vs {want}"
    );
}
