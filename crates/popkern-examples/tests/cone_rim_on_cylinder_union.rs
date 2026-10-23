//! PK-925: a cone whose base rim lies on a coaxial cylinder wall of the
//! same radius unions with the cylinder in either order. With the cone
//! first, the record the cylinder wall traced along the rim and the
//! rim itself sat a period apart on the cone's chart, the duplicate
//! cover was never dropped, the cone wall fell back to the 3D-Newell
//! walk and served the rim as a face of its own: a body with a hole
//! wound like its outer loop (declined since PK-923, served before).
//!
//! PK-929: the same union under rigid rotations of both operands. About
//! z the cone's apex cap flipped its rim pcurve whenever the rim vertex
//! sat off u = 0 (a 1.0 trim miss, declined); about y the disc's clip
//! kept the base rim record one chord short and the wall served split
//! in two. A rigid rotation changes neither the faces nor the volume.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::body_volume;
use popkern_mesh::{TessellationOptions, tessellate};

fn signed_volume(brep: &BRep) -> f64 {
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
fn a_cone_on_its_rim_cylinder_unions_in_either_order_under_rotation() {
    let tol = Tolerance::default();
    // The cylinder, plus the cone's tip above its top (radius 0.5 at
    // z 0.5, height 0.5); the rest of the cone lies inside.
    let want = std::f64::consts::PI * (1.0 + 0.25 * 0.5 / 3.0);
    let rotations = [
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("z pi", Transform3::rotation_z(std::f64::consts::PI)),
        ("z -2", Transform3::rotation_z(-2.0)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
    ];
    for (turn, r) in rotations {
        let cone = BRep::cone(1.0, 1.0, tol)
            .and_then(|c| c.transform(&r))
            .expect("cone");
        let cylinder = BRep::cylinder(1.0, 1.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -0.5)))
            .and_then(|c| c.transform(&r))
            .expect("cylinder");
        for (order, a, b) in [
            ("cone ∪ cylinder", &cone, &cylinder),
            ("cylinder ∪ cone", &cylinder, &cone),
        ] {
            let name = format!("{order}, {turn}");
            assert_union_serves(&name, a, b, want, tol);
        }
    }
}

/// The union serves the four faces, valid, with the exact volume.
fn assert_union_serves(name: &str, a: &BRep, b: &BRep, want: f64, tol: Tolerance) {
    let body = boolean(BoolOperation::Union, a.clone(), b.clone(), tol, None)
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{name}: {errors:?}");
    assert_eq!(
        body.faces.len(),
        4,
        "{name}: bottom cap, wall, top annulus, cone tip"
    );
    let v = signed_volume(&body);
    assert!(
        (v - want).abs() < 2e-3 * want,
        "{name}: signed volume {v} vs {want}"
    );
    let exact = body_volume(&body, &tol).expect("volume");
    assert!(
        (exact - want).abs() < 1e-9 * want,
        "{name}: volume {exact} vs {want}"
    );
}
