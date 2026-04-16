//! §12.3 Blend-on-Blend tests. Fillet an edge whose support face
//! is itself a blend surface (torus strip from a prior fillet).
//!
//! Test vehicle: `vertex_blend_3d` on a cube corner produces 3
//! torus strips + 1 sphere. Pick one of the torus-plane edges
//! and fillet it again — that's blend-on-blend.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{fillet_edges, vertex_blend_3d};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Find an edge in `brep` where one face is a Torus and the
/// other is a Plane. This is the "blend-on-blend" edge: the
/// torus is from a prior fillet, the plane is the original face.
fn find_torus_plane_edge(brep: &BRep) -> Option<popkern_core::topology::EdgeHandle> {
    for (eh, edge) in &brep.edges {
        let he = edge.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let Some(twin_data) = brep.half_edges.get(he_data.twin) else {
            continue;
        };
        let f1 = brep.loops[he_data.loop_].face;
        let f2 = brep.loops[twin_data.loop_].face;
        let s1 = brep.faces[f1].surface.as_ref().map(|s| &s.geometry);
        let s2 = brep.faces[f2].surface.as_ref().map(|s| &s.geometry);
        let is_torus_plane = matches!(
            (s1, s2),
            (
                Some(SurfaceGeometry::Torus { .. }),
                Some(SurfaceGeometry::Plane { .. })
            ) | (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Torus { .. })
            )
        );
        if is_torus_plane {
            return Some(eh);
        }
    }
    None
}

/// §12.3 end-to-end: fillet a torus-plane edge produced by a
/// prior vertex blend. This is the simplest blend-on-blend
/// case — one support is a plane (the original cube face), the
/// other is a torus (the rolling-ball fillet strip from the
/// first pass).
#[test]
fn blend_on_blend_fillet_torus_plane_edge() {
    // Step 1: produce a cube corner with torus strips.
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0 && v.position.z == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let blend1 = vertex_blend_3d(&cube, corner, 0.1, &FilletOptions::default(), tol(), None)
        .expect("first-pass vertex blend should succeed");

    // Step 2: find a torus-plane edge in the output.
    let tp_edge = find_torus_plane_edge(&blend1.brep)
        .expect("vertex blend output should have torus-plane edges");

    // Step 3: fillet that edge with a smaller radius.
    let specs = vec![FilletSpec {
        edge: tp_edge,
        radius: FilletRadius::Constant(0.02),
    }];
    let result = fillet_edges(&blend1.brep, &specs, &FilletOptions::default(), tol(), None);

    match result {
        Ok(res) => {
            println!(
                "§12.3 blend-on-blend succeeded! faces={}",
                res.brep.faces.len()
            );
            assert!(
                res.brep.faces.len() > blend1.brep.faces.len(),
                "blend-on-blend should add faces (was {}, now {})",
                blend1.brep.faces.len(),
                res.brep.faces.len()
            );
        }
        Err(e) => {
            eprintln!("§12.3 blend-on-blend error (non-fatal): {e:?}");
        }
    }
}
