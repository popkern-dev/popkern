//! PK-111 follow-up pin: the boolean kernel emits an OCC-exact
//! cube∪sphere (BRepCheck valid, zero free bounds, volume 2.14682 vs
//! the grid-integrated oracle 2.14679), but the tessellator cannot
//! yet mesh the FULL-SPHERE wall face — a periodic face whose outer
//! loop is the seam-pair spine (two meridian half-edges through both
//! poles) carrying eight triangular imprint holes. Pinned so a mesher
//! fix flips this test instead of hiding behind the boolean suites.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{TessellationOptions, tessellate};

fn union_mesh_open_pairs_and_volume() -> (usize, f64) {
    let t = Tolerance::default().with_linear(1e-2);
    let cube = BRep::unit_cube(t).expect("cube");
    let sphere = BRep::sphere(0.8, t)
        .expect("sphere")
        .transform(&Transform3::translation(0.5, 0.5, 0.5))
        .expect("translate");
    let out = popkern_bool::union(cube, sphere, t, None).expect("union (kernel-exact per OCC)");
    let mesh = tessellate(&out, &TessellationOptions::default()).expect("tessellate");
    let mut counts: std::collections::HashMap<(usize, usize), i64> =
        std::collections::HashMap::new();
    for tri in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += if a < b { 1 } else { -1 };
        }
    }
    let opens = counts.values().filter(|&&v| v != 0).count();
    let volume: f64 = mesh
        .triangles
        .iter()
        .map(|tri| {
            let a = mesh.vertices[tri[0]].position;
            let b = mesh.vertices[tri[1]].position;
            let c = mesh.vertices[tri[2]].position;
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum::<f64>()
        .abs();
    (opens, volume)
}

#[test]
fn sphere_union_wall_face_mesh_gap_is_pinned() {
    let (opens, volume) = union_mesh_open_pairs_and_volume();
    assert!(
        opens > 0,
        "mesher gap CLOSED (0 open pairs) — promote this pin to a closure \
         gate: assert opens == 0 and volume within 5e-3 of 2.14679 \
         (measured {volume:.5})"
    );
}
