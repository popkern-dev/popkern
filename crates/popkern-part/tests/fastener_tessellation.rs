//! PK-560: tessellated thread-flank faces lie on their own surfaces.
//!
//! The angled tapped hole's helical flank faces have trims that are
//! straight uv lines over several turns. Sampled by uv straightness,
//! their trim polylines got two points, and the 3D chord between them
//! crossed the hole. Trim-clip vertices were interpolated along that
//! chord and landed up to 2.57 off the flank, on the hole axis. The
//! watertight gate could not see it, because the whole-body pipeline
//! carries canonical edge stations; only the standalone per-face path
//! showed it. Both paths are gated here, per vertex, against the face's
//! own surface.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_part::fastener::{HoleDepth, TapOptions, ThreadSpec};
use popkern_part::prelude::*;

/// Worst distance from any of `vertices` to `surface` (closest point).
fn off_surface(surface: &SurfaceGeometry, vertices: &[popkern_mesh::MeshVertex]) -> f64 {
    vertices
        .iter()
        .map(|v| {
            surface
                .closest_point(&v.position)
                .map(|(u, w)| surface.point_at(u, w).distance_to(&v.position))
                .unwrap_or(f64::INFINITY)
        })
        .fold(0.0, f64::max)
}

fn nurbs_surface(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> Option<&SurfaceGeometry> {
    let surface = &brep.faces[face].surface.as_ref()?.geometry;
    matches!(surface, SurfaceGeometry::Nurbs(_)).then_some(surface)
}

#[test]
fn angled_tap_flank_meshes_lie_on_their_surfaces() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let plate = Profile::on(Plane::XY)
        .rect(20.0, 20.0)
        .extrude(8.0)
        .expect("block");
    let tapped = plate
        .tapped_hole(
            &spec,
            Point3::new(1.0, -0.5, 8.0),
            Vec3::new(0.35, 0.2, -1.0),
            HoleDepth::Blind(4.5),
            &TapOptions {
                clearance: 0.15,
                entry_chamfer: false,
            },
        )
        .expect("angled tapped hole");
    let brep = tapped.brep();
    let options = popkern_mesh::TessellationOptions::default();
    let gate = 1e-6;
    let mut flanks = 0;
    for face in brep.faces.keys() {
        let Some(surface) = nurbs_surface(brep, face) else {
            continue;
        };
        flanks += 1;
        let mesh = popkern_mesh::tessellate_face(brep, face, &options).expect("face mesh");
        let worst = off_surface(surface, &mesh.vertices);
        assert!(
            worst <= gate,
            "{face:?}: standalone mesh {worst:e} off its surface"
        );
    }
    assert!(
        flanks >= 4,
        "expected the thread's NURBS flanks, found {flanks}"
    );
    let body = tapped.mesh(&options).expect("body mesh");
    for face_mesh in &body.face_meshes {
        if let Some(surface) = nurbs_surface(brep, face_mesh.face) {
            let worst = off_surface(surface, &face_mesh.vertices);
            assert!(
                worst <= gate,
                "{:?}: body mesh {worst:e} off its surface",
                face_mesh.face
            );
        }
    }
}
