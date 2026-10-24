//! PK-936: the trim-grid face mesher (`tessellate_face`, also the
//! fallback for faces the boundary-first mesher declines) on loops that
//! cover a whole period of a periodic chart by passing through its pole.
//!
//! A cone minus a box whose corner edge runs down the cone's axis
//! leaves a wall whose loop runs the rim, down the seam to the apex,
//! and back up the seam's other copy around the notch; the folded
//! convention writes both seam copies at the same u. The trim polygon
//! was folded about its circular mean into a self-crossing scribble and
//! the grid meshed 0.614 of the wall's 6.410 in every orientation. The
//! period now rides the apex row, and the grid widens to the lifted
//! trim.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, SurfaceGeometry, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate_face};

fn area(triangles: &[[usize; 3]], at: &dyn Fn(usize) -> Point3) -> f64 {
    triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(at);
            (b - a).cross(&(c - a)).length() / 2.0
        })
        .sum()
}

fn rotations() -> Vec<(&'static str, Transform3)> {
    vec![
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("z pi", Transform3::rotation_z(std::f64::consts::PI)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
    ]
}

/// The trim-grid area of `body`'s faces on surfaces `kind` picks.
fn grid_area(body: &BRep, kind: fn(&SurfaceGeometry) -> bool) -> f64 {
    let options = TessellationOptions::default();
    body.faces
        .iter()
        .filter(|(_, f)| f.surface.as_ref().is_some_and(|s| kind(&s.geometry)))
        .map(|(h, _)| {
            let m = tessellate_face(body, h, &options).expect("face mesh");
            area(&m.triangles, &|i| m.vertices[i].position)
        })
        .sum()
}

#[test]
fn a_cone_wall_folded_through_its_apex_meshes_its_whole_region() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    // Lateral area π·r·slant less the notch, a quarter of the band
    // between radii 0.1 and 0.6 (area π·slant·r² below radius r).
    let slant = 5.0_f64.sqrt();
    let want = pi * slant * (1.0 - (0.36 - 0.01) / 4.0);
    for (name, r) in rotations() {
        let cone = BRep::cone(1.0, 2.0, tol)
            .and_then(|c| c.transform(&r))
            .expect("cone");
        let cube = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 0.8)))
            .and_then(|c| c.transform(&r))
            .expect("box");
        let body = boolean(BoolOperation::Subtraction, cone, cube, tol, None).expect("sub");
        let got = grid_area(&body, |s| matches!(s, SurfaceGeometry::Cone { .. }));
        assert!(
            (got - want).abs() < 2e-3 * want,
            "{name}: wall area {got} vs {want}"
        );
    }
}

/// The 7/8 sphere face of a cube ∪ a sphere centred on its corner runs
/// a seam spine to the north pole and back, then enters and leaves its
/// octant hole along the south pole row; the grid meshed nothing. The
/// grid's own chord error on a whole sphere is 4.4e-3 here.
#[test]
fn a_sphere_face_through_both_poles_meshes_its_whole_region() {
    let tol = Tolerance::default();
    let want = 7.0 / 8.0 * std::f64::consts::PI;
    for (name, r) in rotations() {
        let cube = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&r))
            .expect("cube");
        let ball = BRep::sphere(0.5, tol)
            .and_then(|s| s.transform(&Transform3::translation(1.0, 1.0, 1.0)))
            .and_then(|s| s.transform(&r))
            .expect("ball");
        let body = boolean(BoolOperation::Union, cube, ball, tol, None).expect("union");
        let got = grid_area(&body, |s| matches!(s, SurfaceGeometry::Sphere { .. }));
        assert!(
            (got - want).abs() < 5e-3 * want,
            "{name}: sphere face area {got} vs {want}"
        );
    }
}
