//! PK-114 gate: an OPEN intersection record genuinely crossing a
//! reflex face's bounds four times must clip to honest per-window
//! records instead of a single crossing hull spanning the
//! out-of-face notch.
//!
//! Fixture: a U-shaped (reflex) plate minus a 45°-rotated slab whose
//! cutting plane meets the plate's top face in the line y + z = 47 —
//! at z = 20 that line enters and leaves BOTH arms of the U
//! (boundary crossings at x = 0, 20, 40, 60), with the span between
//! the arms lying outside the face. Decisive multi-window clipping
//! (transversal record, out-of-face span clear of the boundary)
//! keeps only the two real arm segments; the pre-PK-114 hull ran the
//! record THROUGH the notch. Acceptance is structural and metric:
//! watertight mesh and the analytic volume oracle.

use popkern_part::prelude::*;

fn mesh_volume(mesh: &popkern_mesh::TriMesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let a = mesh.vertices[t[0]].position;
            let b = mesh.vertices[t[1]].position;
            let c = mesh.vertices[t[2]].position;
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum::<f64>()
        .abs()
}

#[test]
fn reflex_face_four_crossing_open_record_clips_to_windows() {
    // U plate: outer 60 × 40 × 20 with a 20 × 25 notch opening
    // toward +y between the arms. Reflex corners at (20, 15) and
    // (40, 15).
    let plate = Profile::on(Plane::XY)
        .polygon(&[
            (0.0, 0.0),
            (60.0, 0.0),
            (60.0, 40.0),
            (40.0, 40.0),
            (40.0, 15.0),
            (20.0, 15.0),
            (20.0, 40.0),
            (0.0, 40.0),
        ])
        .extrude(20.0)
        .expect("U plate");

    // 45° slab: bottom plane z = y after rotate_x(-45°), translated
    // to z = 47 - y. Subtracting removes material above that plane —
    // a wedge across BOTH arms for y in [27, 40].
    let cutter = Profile::on(Plane::XY)
        .rect(400.0, 400.0)
        .at(30.0, 0.0)
        .extrude(400.0)
        .expect("slab")
        .rotate_x(-45.0_f64.to_radians())
        .expect("rotate")
        .translate(0.0, 0.0, 47.0)
        .expect("place");

    let cut = plate
        .subtract(&cutter)
        .expect("wedge cut across the U arms");

    let opens = cut.open_edge_count().expect("mesh");
    // Analytic oracle: base U volume 1900 * 20 = 38000; the wedge
    // removes a triangular prism of cross-section 0.5 * 13^2 over
    // both 20-wide arms: 84.5 * 40 = 3380.
    let expected = 38_000.0 - 3_380.0;
    let mesh = cut
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("tessellation");
    let volume = mesh_volume(&mesh);
    assert!(
        opens == 0 && (volume - expected).abs() <= expected * 1.0e-3,
        "reflex four-crossing open-record clip regressed: opens={opens} volume={volume:.3} \
         expected={expected:.3}"
    );
}
