//! PK-930: a box and a sphere centred on one of its corners, so the
//! box's corner edge runs along the sphere's axis through its pole and
//! the octant cut meets the sphere's seam at that pole. Union,
//! subtraction and intersection, under rigid rotations of both
//! operands.
//!
//! Stacked defects, in pipeline order:
//! - the clip chorded the plane x sphere record's fitted polyline, so
//!   the meridian from the pole read as a diagonal across the chart;
//! - the sphere meridian's polar angle came from acos at the pole (half
//!   its digits gone), so an arc ending on the pole lost its exact
//!   pcurve to a chart fit;
//! - the mesher's pole crossing offered a seam pair over the pole only
//!   the short way, found no winding-0 ring, and its winding-1 fallback
//!   picked the far pole's direction arbitrarily.
//!
//! The union served 0.80 for 1.46, inside out, validator-clean.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::body_volume;
use popkern_mesh::{TessellationOptions, tessellate};

fn mesh_volume(brep: &BRep) -> f64 {
    let options = TessellationOptions {
        chord_height: 1.0e-4,
        ..Default::default()
    };
    let mesh = tessellate(brep, &options).expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.vertices[i].position);
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

#[test]
fn a_sphere_on_a_box_corner_combines_exactly_under_rotation() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    let octant = pi / 6.0 / 8.0; // an eighth of the radius-0.5 sphere
    let rotations = [
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("z pi", Transform3::rotation_z(pi)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
        (
            "mixed",
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
    ];
    let ops = [
        (BoolOperation::Union, 1.0 + 7.0 * octant, 7),
        (BoolOperation::Subtraction, 1.0 - octant, 7),
        (BoolOperation::Intersection, octant, 4),
    ];
    for (name, r) in rotations {
        let cube = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&r))
            .expect("box");
        let sphere = BRep::sphere(0.5, tol)
            .and_then(|s| s.transform(&Transform3::translation(1.0, 1.0, 1.0)))
            .and_then(|s| s.transform(&r))
            .expect("sphere");
        for (op, want, faces) in ops {
            let tag = format!("{op:?}, {name}");
            let body = boolean(op, cube.clone(), sphere.clone(), tol, None)
                .unwrap_or_else(|e| panic!("{tag}: {e}"));
            let errors = body.validate_all();
            assert!(errors.is_empty(), "{tag}: {errors:?}");
            assert_eq!(body.faces.len(), faces, "{tag}: face count");
            let exact = body_volume(&body, &tol).expect("volume");
            assert!(
                (exact - want).abs() < 1e-9 * want,
                "{tag}: volume {exact} vs {want}"
            );
            let meshed = mesh_volume(&body);
            assert!(
                (meshed - want).abs() < 2e-3 * want,
                "{tag}: mesh volume {meshed} vs {want}"
            );
        }
    }
}
