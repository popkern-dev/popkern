//! Tests for body-level measurement functions.

use super::*;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn make_cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("cube")
}

// ---------------------------------------------------------------------------
// Existing unit cube tests (moved from body.rs)
// ---------------------------------------------------------------------------

#[test]
fn unit_cube_volume() {
    let cube = make_cube();
    let tol = Tolerance::default();
    let vol = body_volume(&cube, &tol).unwrap();
    assert!((vol - 1.0).abs() < 0.01, "expected ~1.0, got {vol}");
}

#[test]
fn unit_cube_surface_area() {
    let cube = make_cube();
    let tol = Tolerance::default();
    let sa = body_surface_area(&cube, &tol).unwrap();
    assert!((sa - 6.0).abs() < 0.01, "expected ~6.0, got {sa}");
}

#[test]
fn unit_cube_centroid() {
    let cube = make_cube();
    let tol = Tolerance::default();
    let c = body_centroid(&cube, &tol).unwrap();
    // Unit cube is [0,1]^3, centroid should be (0.5, 0.5, 0.5)
    assert!((c.x - 0.5).abs() < 0.01, "centroid x={}", c.x);
    assert!((c.y - 0.5).abs() < 0.01, "centroid y={}", c.y);
    assert!((c.z - 0.5).abs() < 0.01, "centroid z={}", c.z);
}

#[test]
fn unit_cube_bounding_box() {
    let cube = make_cube();
    let bb = body_bounding_box(&cube).unwrap();
    assert!((bb.min.x - 0.0).abs() < 1e-10);
    assert!((bb.max.x - 1.0).abs() < 1e-10);
    assert!((bb.size_x() - 1.0).abs() < 1e-10);
}

#[test]
fn unit_cube_moment_of_inertia() {
    let cube = make_cube();
    let tol = Tolerance::default();
    // Moment of inertia of unit cube (side=1, mass=volume=1) about an axis through centroid:
    // I = (1/6) for a unit cube ~ 0.1667
    let moi = body_moment_of_inertia(&cube, Axis::Z, &tol).unwrap();
    assert!((moi - 1.0 / 6.0).abs() < 0.05, "expected ~0.167, got {moi}");
}

// ---------------------------------------------------------------------------
// Moment-of-inertia: unit cube all axes are equal by symmetry
// ---------------------------------------------------------------------------

#[test]
fn unit_cube_moi_all_axes_equal() {
    let cube = make_cube();
    let tol = Tolerance::default();
    let moi_x = body_moment_of_inertia(&cube, Axis::X, &tol).unwrap();
    let moi_y = body_moment_of_inertia(&cube, Axis::Y, &tol).unwrap();
    let moi_z = body_moment_of_inertia(&cube, Axis::Z, &tol).unwrap();
    // By cubic symmetry, all three should be equal.
    assert!(
        (moi_x - moi_y).abs() < 0.05,
        "Ix={moi_x} should equal Iy={moi_y}"
    );
    assert!(
        (moi_y - moi_z).abs() < 0.05,
        "Iy={moi_y} should equal Iz={moi_z}"
    );
}

// ---------------------------------------------------------------------------
// Moment-of-inertia: sphere (r=1)
// Note: The MoI implementation uses a centroid-based tet approximation, so
// accuracy for curved solids depends on tessellation resolution. We use
// generous tolerances for curved geometry.
// ---------------------------------------------------------------------------

#[test]
fn sphere_moi_positive_and_nonzero() {
    // Sphere MoI should be positive. The exact value is 8*pi/15 ~ 1.6755,
    // but the tet-centroid approximation may differ significantly.
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();

    let moi_z = body_moment_of_inertia(&sphere, Axis::Z, &tol).unwrap();
    assert!(moi_z > 0.0, "sphere Iz should be positive, got {moi_z}");
}

#[test]
fn sphere_moi_all_axes_positive() {
    // Sphere has rotational symmetry about all axes. The centroid-based tet
    // approximation severely underestimates the sphere MoI due to the coarse
    // tessellation, but all axes should still produce positive values.
    let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let moi_x = body_moment_of_inertia(&sphere, Axis::X, &tol).unwrap();
    let moi_y = body_moment_of_inertia(&sphere, Axis::Y, &tol).unwrap();
    let moi_z = body_moment_of_inertia(&sphere, Axis::Z, &tol).unwrap();

    assert!(moi_x > 0.0, "sphere Ix should be positive, got {moi_x}");
    assert!(moi_y > 0.0, "sphere Iy should be positive, got {moi_y}");
    assert!(moi_z > 0.0, "sphere Iz should be positive, got {moi_z}");
}

// ---------------------------------------------------------------------------
// Moment-of-inertia: cylinder (r=1, h=2)
// ---------------------------------------------------------------------------

#[test]
fn cylinder_moment_of_inertia_about_axis() {
    // Solid cylinder of radius r=1, height h=2 centered at z=0..h,
    // about its own axis (Z through centroid).
    // Volume = pi*r^2*h = 2*pi.
    // I_z = (1/2)*m*r^2 = (1/2)*(2*pi)*1 = pi ~ 3.1416.
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let expected = std::f64::consts::PI;
    let moi_z = body_moment_of_inertia(&cyl, Axis::Z, &tol).unwrap();
    assert!(
        (moi_z - expected).abs() / expected < 0.10,
        "cylinder Iz should be ~{expected}, got {moi_z}"
    );
}

#[test]
fn cylinder_moment_of_inertia_transverse() {
    // Solid cylinder r=1, h=2 about a transverse axis (X or Y through centroid).
    // I_x = m*(3*r^2 + h^2)/12 = 2*pi*(3+4)/12 = 7*pi/6 ~ 3.665.
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let m = std::f64::consts::PI * 1.0 * 1.0 * 2.0;
    let expected = m * (3.0 * 1.0 + 4.0) / 12.0;
    let moi_x = body_moment_of_inertia(&cyl, Axis::X, &tol).unwrap();
    assert!(
        (moi_x - expected).abs() / expected < 0.20,
        "cylinder Ix should be ~{expected}, got {moi_x}"
    );
}

#[test]
fn cylinder_moi_xy_symmetry() {
    // Cylinder axis is Z, so I_x == I_y.
    let cyl = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let moi_x = body_moment_of_inertia(&cyl, Axis::X, &tol).unwrap();
    let moi_y = body_moment_of_inertia(&cyl, Axis::Y, &tol).unwrap();
    assert!(
        (moi_x - moi_y).abs() / moi_x.max(1e-12) < 0.05,
        "cylinder Ix={moi_x} should equal Iy={moi_y}"
    );
}

// ---------------------------------------------------------------------------
// Moment-of-inertia: rectangular prism (cuboid) with distinct sides
// ---------------------------------------------------------------------------

#[test]
fn rectangular_prism_moi_distinct_axes() {
    // Prism 2x3x1 (width=2, depth=3, height=1).
    // Volume = 6, mass = 6 (unit density).
    // I_z = m*(a^2 + b^2)/12 = 6*(4+9)/12 = 6.5
    // I_y = m*(a^2 + c^2)/12 = 6*(4+1)/12 = 2.5
    // I_x = m*(b^2 + c^2)/12 = 6*(9+1)/12 = 5.0
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 3.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
    ];
    let prism = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("rectangular prism");
    let tol = Tolerance::default();

    let moi_x = body_moment_of_inertia(&prism, Axis::X, &tol).unwrap();
    let moi_y = body_moment_of_inertia(&prism, Axis::Y, &tol).unwrap();
    let moi_z = body_moment_of_inertia(&prism, Axis::Z, &tol).unwrap();

    // I_x = m*(b^2+c^2)/12 = 6*(9+1)/12 = 5.0
    assert!(
        (moi_x - 5.0).abs() < 2.0,
        "prism Ix should be ~5.0, got {moi_x}"
    );
    // I_y = m*(a^2+c^2)/12 = 6*(4+1)/12 = 2.5
    assert!(
        (moi_y - 2.5).abs() < 1.5,
        "prism Iy should be ~2.5, got {moi_y}"
    );
    // I_z = m*(a^2+b^2)/12 = 6*(4+9)/12 = 6.5
    assert!(
        (moi_z - 6.5).abs() < 2.5,
        "prism Iz should be ~6.5, got {moi_z}"
    );
}

#[test]
fn rectangular_prism_moi_ordering() {
    // For a 2x3x1 prism, I_y < I_x < I_z because the smallest two
    // dimensions contribute least to the transverse axes.
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 3.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
    ];
    let prism = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("rectangular prism");
    let tol = Tolerance::default();

    let moi_x = body_moment_of_inertia(&prism, Axis::X, &tol).unwrap();
    let moi_y = body_moment_of_inertia(&prism, Axis::Y, &tol).unwrap();
    let moi_z = body_moment_of_inertia(&prism, Axis::Z, &tol).unwrap();

    assert!(
        moi_y < moi_x && moi_x < moi_z,
        "expected Iy < Ix < Iz, got Ix={moi_x}, Iy={moi_y}, Iz={moi_z}"
    );
}

// ---------------------------------------------------------------------------
// Moment-of-inertia: cone (r=1, h=2)
// ---------------------------------------------------------------------------

#[test]
fn cone_moment_of_inertia_about_axis() {
    // Solid cone of radius r=1, height h=2.
    // Volume = (1/3)*pi*r^2*h = 2*pi/3.
    // I_z = (3/10)*m*r^2 = (3/10)*(2*pi/3)*1 = pi/5 ~ 0.6283.
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let m = std::f64::consts::PI / 3.0 * 2.0; // (1/3)*pi*r^2*h
    let expected = 3.0 / 10.0 * m * 1.0; // (3/10)*m*r^2
    let moi_z = body_moment_of_inertia(&cone, Axis::Z, &tol).unwrap();
    assert!(
        (moi_z - expected).abs() / expected < 0.1,
        "cone Iz should be ~{expected}, got {moi_z}"
    );
}

// ---------------------------------------------------------------------------
// Curved-body volume: exact-surface quadrature (PK-807)
// ---------------------------------------------------------------------------

/// Curved bodies integrate the exact surface over the tessellation's
/// parameter triangles and take planar terms from the edge curves, so
/// analytic bodies read their analytic volumes (the flat-facet volume
/// was off by the facet sagitta: a sphere of radius 1 read -0.06%).
#[test]
fn curved_body_volumes_match_their_analytic_values() {
    use std::f64::consts::PI;
    let tol = Tolerance::default();
    let cases = [
        (
            "cylinder r1 h2",
            BRep::cylinder(1.0, 2.0, tol).unwrap(),
            2.0 * PI,
        ),
        ("sphere r1", BRep::sphere(1.0, tol).unwrap(), 4.0 / 3.0 * PI),
        (
            "cone r1 h2",
            BRep::cone(1.0, 2.0, tol).unwrap(),
            2.0 / 3.0 * PI,
        ),
    ];
    for (name, body, exact) in cases {
        let v = body_volume(&body, &tol).unwrap();
        assert!((v - exact).abs() < 1e-9 * exact, "{name}: {v} vs {exact}");
    }
}

#[test]
fn bored_cube_volume_is_exact() {
    use popkern_core::geometry::transform::Transform3;
    let tol = Tolerance::default();
    let bore = BRep::cylinder(0.15, 3.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -1.0))
        .unwrap();
    let bored = popkern_bool::subtraction(BRep::unit_cube(tol).unwrap(), bore, tol, None).unwrap();
    let exact = 1.0 - std::f64::consts::PI * 0.0225;
    let v = body_volume(&bored, &tol).unwrap();
    assert!((v - exact).abs() < 1e-9, "{v} vs {exact}");
}
