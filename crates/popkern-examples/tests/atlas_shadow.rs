//! The mesher's atlas in its shadow phase (chart-service.md §6a, PK-1267):
//! faces reaching a sphere's poles meshed with `TessellationOptions::atlas`
//! (each piece in a chart with no singular row) must be watertight, and
//! every sphere face's facet area must lie within the inscribed band of
//! its exact area (under it by at most the chord deficit, never over: a
//! piece meshed twice or a hole meshed as material reads high).
//!
//! Each scene was a recorded failing state of the atlas as built: the
//! turned Viviani union's figure-eight lobes came first among the band's
//! rings and its hole meshed as material (1247.7 for 1028.3); the
//! midsphere's patches read their row through a boundary contact and
//! meshed the whole sphere (618.5 for 9.5); the cube-and-sphere union's
//! caps sat at the chord height and its volume read 0.27% low.

use std::f64::consts::{FRAC_PI_2, PI};

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::TessellationOptions;

fn atlas(display: bool) -> TessellationOptions {
    let base = if display {
        TessellationOptions::display()
    } else {
        TessellationOptions::default()
    };
    TessellationOptions {
        atlas: true,
        ..base
    }
}

/// Watertight, and each sphere face's facet area within the inscribed
/// band of its exact area: under it by at most the chord deficit (twice
/// the chord height over the radius, relative), over it by at most `over`
/// (1e-6; the midsphere's patches, which do not reach a pole and so do
/// not take the atlas, read 1.5e-3 high on both paths, PK-1272).
fn assert_meshes(name: &str, body: &BRep, tol: &Tolerance, over: f64) {
    for display in [false, true] {
        let mesh = popkern_mesh::tessellate(body, &atlas(display)).expect("mesh");
        assert_eq!(
            popkern_mesh::welded_open_edges(&mesh),
            0,
            "{name} (display {display}): watertight"
        );
        for fm in &mesh.face_meshes {
            let on_sphere = body
                .faces
                .get(fm.face)
                .and_then(|f| f.surface.as_ref())
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Sphere { .. }));
            if !on_sphere {
                continue;
            }
            let area: f64 = fm
                .triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| fm.vertices[i].position);
                    0.5 * (b - a).cross(&(c - a)).length()
                })
                .sum();
            let exact = popkern_measure::face_area(body, fm.face, tol).expect("area");
            let rel = (area - exact) / exact;
            let radius = match body.faces[fm.face].surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Sphere { radius, .. }) => *radius,
                _ => f64::INFINITY,
            };
            let chord = atlas(display).chord_height;
            assert!(
                (-2.0 * chord / radius..over).contains(&rel),
                "{name} (display {display}): face {:?} mesh area {area} for {exact}",
                fm.face
            );
        }
    }
}

fn viviani(tol: Tolerance) -> (BRep, BRep) {
    let ball = BRep::sphere(10.0, tol)
        .and_then(|b| b.transform(&Transform3::rotation_x(-FRAC_PI_2)))
        .and_then(|b| b.transform(&Transform3::rotation_y(PI)))
        .expect("ball");
    let rod = BRep::cylinder(5.0, 24.0, tol)
        .and_then(|b| b.transform(&Transform3::rotation_z(PI)))
        .and_then(|b| b.transform(&Transform3::translation(5.0, 0.0, -12.0)))
        .expect("rod");
    (ball, rod)
}

#[test]
fn turned_viviani_operands_mesh_on_the_atlas() {
    let tol = Tolerance::default();
    for (name, op) in [
        ("union", popkern_bool::BoolOperation::Union),
        ("common", popkern_bool::BoolOperation::Intersection),
    ] {
        let (ball, rod) = viviani(tol);
        let body = popkern_bool::boolean(op, ball, rod, tol, None).expect("boolean");
        assert_meshes(&format!("viviani {name}"), &body, &tol, 1e-6);
    }
}

#[test]
fn the_midsphere_patches_mesh_on_the_atlas() {
    let tol = Tolerance::default();
    let rect = |z: f64| {
        vec![
            Point3::new(-5.0, -5.0, z),
            Point3::new(5.0, -5.0, z),
            Point3::new(5.0, 5.0, z),
            Point3::new(-5.0, 5.0, z),
        ]
    };
    let cube = BRep::polygon_frustum(&rect(-5.0), &rect(5.0), tol).expect("cube");
    let ball = BRep::sphere(5.0 * 2f64.sqrt(), tol).expect("ball");
    let body = popkern_bool::boolean(
        popkern_bool::BoolOperation::Intersection,
        cube,
        ball,
        tol,
        None,
    )
    .expect("common");
    assert_meshes("midsphere", &body, &tol, 2e-3);
}

#[test]
fn a_cube_and_sphere_union_meshes_on_the_atlas() {
    let tol = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(tol).expect("cube");
    let sphere = BRep::sphere(0.8, tol)
        .and_then(|b| b.transform(&Transform3::translation(0.5, 0.5, 0.5)))
        .expect("sphere");
    let body = popkern_bool::union(cube, sphere, tol, None).expect("union");
    assert_meshes("cube and sphere", &body, &tol, 1e-6);
}

/// A ball of radius 10 less the half-space beyond the plane n·x = d
/// (n = (0.3, −0.5, 0.2) normalised, d = −6 scaled to it): the rest of
/// the sphere beyond a cap clear of both poles. Its band's ring wraps the
/// period and the band's slit grazed the cut circle 7e-4 away, 4 edges
/// open (the slit now keeps clear of the face's other rings).
#[test]
fn the_rest_of_a_sphere_beyond_a_cap_meshes_on_the_atlas() {
    let tol = Tolerance::default();
    let raw = [0.3_f64, -0.5, 0.2];
    let len = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
    let n = popkern_core::geometry::Vec3::new(raw[0] / len, raw[1] / len, raw[2] / len);
    let d = -6.0 / len;
    let e1 = n
        .cross(&popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0))
        .normalized()
        .expect("frame");
    let e2 = n.cross(&e1);
    let corner = |a: f64, b: f64, h: f64| {
        let p = n * (d + h) + e1 * a + e2 * b;
        Point3::new(p.x, p.y, p.z)
    };
    let rect = |h: f64| {
        vec![
            corner(-30.0, -30.0, h),
            corner(30.0, -30.0, h),
            corner(30.0, 30.0, h),
            corner(-30.0, 30.0, h),
        ]
    };
    let block = BRep::polygon_frustum(&rect(-30.0), &rect(0.0), tol).expect("block");
    let ball = BRep::sphere(10.0, tol).expect("ball");
    let body = popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        ball,
        block,
        tol,
        None,
    )
    .expect("difference");
    assert_meshes("sphere beyond a cap", &body, &tol, 1e-6);
}
