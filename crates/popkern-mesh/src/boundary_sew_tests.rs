//! PK-458: the projection sew must never touch an internally TORN
//! face mesh. A tear turns interior vertices into "boundary", and the
//! generous projection band (20% of the face diagonal — larger than a
//! shallow face's own depth) then mass-snaps them onto the rim
//! polylines (measured: 2168 vertices of a broken dimple cap
//! collapsed onto its rim circle). The tear is structural — more
//! boundary components than BRep loops — and a torn face is skipped
//! so the open-edge census reports the genuine defect.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::sew_boundaries_to_edge_polylines;
use crate::TessellationOptions;

fn subdivided_cube_face(brep: &BRep) -> crate::mesh::FaceMesh {
    let face = brep.faces.keys().next().expect("cube face");
    let options = TessellationOptions {
        chord_height: 1.0e-3,
        min_subdivisions: 8,
        max_subdivisions: 8,
        ..Default::default()
    };
    crate::tessellate_face(brep, face, &options).expect("face mesh")
}

#[test]
fn torn_face_is_not_sewn_and_healthy_face_still_is() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let healthy = subdivided_cube_face(&brep);

    // Tear the interior: drop one triangle whose corners are all off
    // the face border (its vertices become mesh-boundary without
    // being anywhere near the BRep boundary polylines' true role).
    let mut torn = healthy.clone();
    let interior = torn
        .triangles
        .iter()
        .position(|t| {
            t.iter().all(|&i| {
                let (u, v) = torn.vertices[i].uv;
                u > 0.3 && u < 0.7 && v > 0.3 && v < 0.7
            })
        })
        .expect("an interior triangle at 8x8 subdivision");
    torn.triangles.remove(interior);

    let before: Vec<_> = torn.vertices.iter().map(|v| v.position).collect();
    let mut meshes = vec![torn, healthy];
    sew_boundaries_to_edge_polylines(&brep, &mut meshes, &[0, 1]);

    // The torn face (index 0) must be untouched — every vertex where
    // it was, including the tear's rim.
    for (v, orig) in meshes[0].vertices.iter().zip(before.iter()) {
        assert_eq!(
            v.position.x.to_bits(),
            orig.x.to_bits(),
            "torn face vertex moved"
        );
        assert_eq!(v.position.y.to_bits(), orig.y.to_bits());
        assert_eq!(v.position.z.to_bits(), orig.z.to_bits());
    }
    // The healthy face (index 1) remains a sew candidate: the pass
    // ran over it without declining (its boundary verts end ON the
    // edge polylines — for a planar cube face they already are, so
    // positions stay on the face plane).
    for v in &meshes[1].vertices {
        assert!(
            v.position.x.abs() <= 1.0 + 1.0e-9
                && v.position.y.abs() <= 1.0 + 1.0e-9
                && v.position.z.abs() <= 1.0 + 1.0e-9,
            "healthy face vertex left the cube"
        );
    }
}
