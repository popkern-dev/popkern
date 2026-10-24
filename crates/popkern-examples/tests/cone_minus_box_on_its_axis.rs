//! PK-931: a cone minus a box whose corner edge runs down the cone's
//! axis, under rigid rotations of both operands.
//!
//! The box's two faces through the axis cut the cone along generator
//! lines; one lies on the cone's own seam. Four defects stacked here:
//! - plane x cone through the apex fell to the marching solver, whose
//!   fitted generator straddled the chart seam once the cone's frame put
//!   u = 0 on it (exact lines now, `plane_cone_apex`);
//! - a record on the cone's seam line read inside or outside the cone
//!   face by rotation (the clip now decides it off the seam);
//! - the cut quadrant's classification sample sat on its own seam
//!   boundary, on the box face, and read Outside (a sample clear of the
//!   other body now decides);
//! - the cone wall's loop runs rim, seam up, apex, seam down around the
//!   notch: the mesher declined the structured revolution and fell to
//!   the uniform grid, meshing 0.18 of the wall's 2.05 volume share —
//!   every orientation served a validator-clean 0.1126 for 1.9818 (the
//!   ring now closes along the apex row).

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
fn a_cone_minus_a_box_on_its_axis_serves_its_volume_under_rotation() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    // The cone (radius 1, height 2) less the quarter frustum the box
    // [0,1]^2 x [0.8, 1.8] removes (radii 0.6 to 0.1 over height 1).
    let want = 2.0 * pi / 3.0 - pi * (0.36 + 0.06 + 0.01) / 12.0;
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
    for (name, r) in rotations {
        let cone = BRep::cone(1.0, 2.0, tol)
            .and_then(|c| c.transform(&r))
            .expect("cone");
        let cube = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 0.8)))
            .and_then(|c| c.transform(&r))
            .expect("box");
        let body = boolean(BoolOperation::Subtraction, cone, cube, tol, None)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let errors = body.validate_all();
        assert!(errors.is_empty(), "{name}: {errors:?}");
        assert_eq!(body.faces.len(), 6, "{name}: wall, base, four box faces");
        let exact = body_volume(&body, &tol).expect("volume");
        assert!(
            (exact - want).abs() < 1e-9 * want,
            "{name}: volume {exact} vs {want}"
        );
        let meshed = mesh_volume(&body);
        assert!(
            (meshed - want).abs() < 2e-3 * want,
            "{name}: mesh volume {meshed} vs {want}"
        );
    }
}
