// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use super::*;
use popkern_core::{
    geometry::Point3, geometry::curve::CurveEval, geometry::nurbs_curve::NurbsCurve,
    tolerance::Tolerance, topology::BRep,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Build a small closed circle in UV space centered at (cx, cy) with radius r.
/// Uses 9 rational NURBS control points (degree 2) for a full circle.
fn make_uv_circle(cx: f64, cy: f64, r: f64) -> NurbsCurve {
    // Full circle via rational quadratic NURBS (9 control points, degree 2).
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let control_points = vec![
        Point3::new(cx + r, cy, 0.0),
        Point3::new(cx + r, cy + r, 0.0),
        Point3::new(cx, cy + r, 0.0),
        Point3::new(cx - r, cy + r, 0.0),
        Point3::new(cx - r, cy, 0.0),
        Point3::new(cx - r, cy - r, 0.0),
        Point3::new(cx, cy - r, 0.0),
        Point3::new(cx + r, cy - r, 0.0),
        Point3::new(cx + r, cy, 0.0), // closes the circle
    ];
    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    NurbsCurve {
        degree: 2,
        knots,
        control_points,
        weights,
    }
}

// -- ray_cast_contains tests --

#[test]
fn ray_cast_unit_square_center_inside() {
    let polygon = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    assert!(ray_cast_contains(&polygon, 0.5, 0.5));
}

#[test]
fn ray_cast_unit_square_outside() {
    let polygon = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    assert!(!ray_cast_contains(&polygon, 2.0, 0.5));
    assert!(!ray_cast_contains(&polygon, 0.5, 2.0));
    assert!(!ray_cast_contains(&polygon, -1.0, 0.5));
}

#[test]
fn ray_cast_degenerate_polygon() {
    // Fewer than 3 points -> always false.
    let polygon = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    assert!(!ray_cast_contains(&polygon, 0.5, 0.0));
}

#[test]
fn ray_cast_triangle() {
    let polygon = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(1.0, 2.0, 0.0),
    ];
    assert!(ray_cast_contains(&polygon, 1.0, 0.5));
    assert!(!ray_cast_contains(&polygon, 3.0, 0.5));
}

// -- point_inside_face_uv tests --

#[test]
fn point_inside_cube_face() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    // The cube face UV coordinates should form a [0,1]^2 square.
    // Pick any face and test a point at (0.5, 0.5) in UV space.
    let face = cube.faces.keys().next().expect("has faces");
    let inside = Point3::new(0.5, 0.5, 0.0);
    assert!(point_inside_face_uv(&cube, face, inside, &tol()));
}

#[test]
fn point_outside_cube_face() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("has faces");
    // Way outside the face UV boundary.
    let outside = Point3::new(10.0, 10.0, 0.0);
    assert!(!point_inside_face_uv(&cube, face, outside, &tol()));
}

// -- inner loop integration tests --

#[test]
fn closed_curve_inside_face_creates_inner_loop() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("has faces");

    // Build a small circle in 3D space (z=0 plane for the bottom face).
    // The exact position doesn't matter for topology — we just need it closed.
    let curve_3d = make_uv_circle(0.5, 0.5, 0.1);
    let curve_uv = make_uv_circle(0.5, 0.5, 0.1);

    let inner_loops_before = cube.faces[face].inner_loops.len();
    let result = cube
        .add_inner_loop(face, &curve_3d, &curve_uv, 12, &tol())
        .expect("add_inner_loop");
    let inner_loops_after = cube.faces[face].inner_loops.len();

    assert_eq!(inner_loops_after, inner_loops_before + 1);
    assert_eq!(result.edges.len(), 12);
    assert_eq!(result.vertices.len(), 12);
}

#[test]
fn closed_curve_outside_face_is_skipped() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("has faces");

    // Circle centered far outside the face boundary.
    let curve_uv = make_uv_circle(10.0, 10.0, 0.1);
    let mid_t = {
        let (t0, t1) = curve_uv.parameter_range();
        (t0 + t1) / 2.0
    };
    let mid_uv = curve_uv.point_at(mid_t);

    // The point-in-face test should return false.
    assert!(!point_inside_face_uv(&cube, face, mid_uv, &tol()));
}
