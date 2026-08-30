//! PK-387 — a boss sketched ON a face of the active body and extruded in
//! Boss mode merges into that body and the recorder returns the host.
//!
//! An on-face sketch is its own tool body (a sketch executes by
//! replacing its body's BRep with the profile face; placed in the host
//! body it erased the host). The Boss/Cut auto-merge recovers the host
//! through the on-face reference, and `Body::extrude_with` returns the
//! body the feature merged into.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_part::Plane as PartPlane;
use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn signature(part: &Part) -> (usize, usize, usize, Vec<[i64; 3]>) {
    let b = part.brep();
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut verts: Vec<[i64; 3]> = mesh
        .vertices
        .iter()
        .map(|v| {
            [
                (v.position.x * 1e6).round() as i64,
                (v.position.y * 1e6).round() as i64,
                (v.position.z * 1e6).round() as i64,
            ]
        })
        .collect();
    verts.sort();
    verts.dedup();
    (b.faces.len(), b.edges.len(), mesh.triangles.len(), verts)
}

#[test]
fn on_face_boss_merges_into_the_plate_and_returns_the_host() {
    let mut m = Model::new();
    let _plate = m
        .profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.4.into()))
        .expect("extrude");
    let boss = m
        .profile(Plane::on_face(faces().facing(Dir::Z).at_height(0.4)))
        .circle(0.35)
        .sketch("boss_sketch")
        .expect("boss sketch")
        .extrude_with(&mut m, End::Blind(0.6.into()), OperationMode::Boss, "boss")
        .expect("boss");
    assert!(boss.is_alive(&m), "returned body is the served host");
    assert!(
        m.alive_bodies().len() == 1,
        "the tool body was consumed: {:?}",
        m.alive_bodies()
    );
    let part = m.to_part().expect("serving boundary");
    let brep = part.brep();
    let cylinders = brep
        .faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .count();
    assert_eq!(cylinders, 1, "one periodic boss wall");
    assert_eq!(part.face_count(), 8, "plate 6 + boss wall + boss top");
    // Parity with the direct build at the frame's circle centre.
    let (cx, cy) = brep
        .faces
        .values()
        .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Cylinder { origin, .. }) => Some((origin.x, origin.y)),
            _ => None,
        })
        .expect("boss wall");
    let direct_plate = Profile::on(PartPlane::XY)
        .rect(2.4, 1.6)
        .extrude(0.4)
        .expect("plate");
    let direct_boss = Profile::on(PartPlane::xy_at(0.4))
        .circle(0.35)
        .at(cx, cy)
        .extrude(0.6)
        .expect("boss");
    let direct = direct_plate.union(&direct_boss).expect("union");
    assert_eq!(
        signature(&part),
        signature(&direct),
        "recorded on-face boss ≡ direct build"
    );
}
