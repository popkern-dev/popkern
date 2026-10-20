// This file is only compiled under #[cfg(test)] (see cylinder_cone.rs).
//! Parallel-offset cylinder x cone pairs (PK-876): a tapped hole's
//! entry chamfer cone beside a flange rim.
use super::*;

fn rim(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

/// A 45° chamfer cone on the axis x = 35, apex at z = 12 opening down
/// (-z): radius v at depth v below the apex.
fn chamfer() -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(35.0, 0.0, 12.0),
        axis: Vec3::new(0.0, 0.0, -1.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    }
}

fn solve(cyl: &SurfaceGeometry, v_cone: (f64, f64)) -> Option<usize> {
    try_cylinder_cone(
        cyl,
        (0.0, std::f64::consts::TAU),
        (-5.0, 20.0),
        &chamfer(),
        (0.0, std::f64::consts::TAU),
        v_cone,
        &Tolerance::default(),
    )
    .map(|c| c.len())
}

#[test]
fn a_chamfer_frustum_inside_the_rim_is_exactly_empty() {
    // Radii 1..3 at 35 off the axis: 38 < 39.65.
    assert_eq!(solve(&rim(39.65), (1.0, 3.0)), Some(0));
    // Wholly outside a smaller cylinder: 35 - 3 > 20.
    assert_eq!(solve(&rim(20.0), (1.0, 3.0)), Some(0));
}

#[test]
fn a_crossing_frustum_is_left_to_the_solver() {
    // Radius up to 6 reaches 41 > 39.65: the rim crosses the cone.
    assert_eq!(solve(&rim(39.65), (1.0, 6.0)), None);
    let found = crate::geometry::intersection::surface_surface::intersect_surfaces(
        &rim(39.65),
        (0.0, std::f64::consts::TAU),
        (-5.0, 20.0),
        &chamfer(),
        (0.0, std::f64::consts::TAU),
        (1.0, 6.0),
        &Tolerance::default(),
    )
    .expect("ssi");
    assert!(
        !found.intersections.is_empty(),
        "the crossing is still traced"
    );
}

#[test]
fn the_window_through_the_apex_has_a_zero_radius() {
    // A window through the apex: radii 0..3. A large cylinder nearly
    // on the cone's axis still holds it; a thin one cannot be wholly
    // inside, outside or around it.
    let hub = SurfaceGeometry::Cylinder {
        origin: Point3::new(35.5, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 10.0,
    };
    // Offset 0.5 + r_max 3 < 10: inside, empty.
    assert_eq!(solve(&hub, (-1.0, 3.0)), Some(0));
    assert!(!parallel_offset_disjoint(
        0.1,
        0.5,
        std::f64::consts::FRAC_PI_4,
        (-1.0, 3.0),
        &Tolerance::default()
    ));
}
