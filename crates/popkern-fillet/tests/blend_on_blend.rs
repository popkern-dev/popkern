//! §12.3 Blend-on-Blend tests. Fillet an edge whose support face
//! is itself a blend surface (cylinder strip from a prior fillet).
//!
//! Test vehicle: `vertex_blend_3d` on a cube corner produces 3
//! cylinder strips + 1 sphere. Pick one of the cylinder-plane
//! edges and fillet it again — that's blend-on-blend.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{fillet_edges, vertex_blend_3d};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Find a short cylinder-plane edge in `brep` where the two
/// surfaces meet at a non-tangent dihedral angle. The long rim
/// edges are tangent (cylinder was designed to be C1-continuous
/// with the plane there), so a fillet on those degenerates. The
/// short end-cap edges at the far cube corners are the true
/// blend-on-blend candidates where the cylinder surface meets a
/// perpendicular plane face.
fn find_cylinder_plane_edge(brep: &BRep) -> Option<popkern_core::topology::EdgeHandle> {
    let mut best: Option<(popkern_core::topology::EdgeHandle, f64)> = None;
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
        let is_cyl_plane = matches!(
            (s1, s2),
            (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            ) | (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        );
        if is_cyl_plane {
            let v1 = brep.vertices[he_data.origin].position;
            let v2 = brep.vertices[brep.half_edges[he_data.twin].origin].position;
            let len =
                ((v1.x - v2.x).powi(2) + (v1.y - v2.y).powi(2) + (v1.z - v2.z).powi(2)).sqrt();
            if best.is_none_or(|(_, bl)| len < bl) {
                best = Some((eh, len));
            }
        }
    }
    best.map(|(eh, _)| eh)
}

/// §12.3 end-to-end: fillet a cylinder-plane edge produced by a
/// prior vertex blend. This is the simplest blend-on-blend
/// case — one support is a plane (the original cube face), the
/// other is a cylinder (the rolling-ball fillet strip from the
/// first pass).
#[test]
fn blend_on_blend_fillet_torus_plane_edge() {
    // Step 1: produce a cube corner with cylinder strips.
    let cube = BRep::unit_cube(tol()).unwrap();
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0 && v.position.z == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let blend1 = vertex_blend_3d(&cube, corner, 0.1, &FilletOptions::default(), tol(), None)
        .expect("first-pass vertex blend should succeed");

    // Step 2: find a cylinder-plane edge in the output.
    let cp_edge = find_cylinder_plane_edge(&blend1.brep)
        .expect("vertex blend output should have cylinder-plane edges");

    // Step 3: fillet that edge with a smaller radius.
    let specs = vec![FilletSpec {
        edge: cp_edge,
        radius: FilletRadius::Constant(0.02),
    }];
    let _res = fillet_edges(&blend1.brep, &specs, &FilletOptions::default(), tol(), None)
        .expect("§12.3 blend-on-blend should succeed");
}
