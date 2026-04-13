//! Round 5r+ post-Phase-C cleanup: merge triangular plane
//! fragments into adjacent sphere vertex-blend faces.
//!
//! Multi-chord Phase B leaves ≤4-vert triangular plane fragments
//! at cube corners — one per cube face per cube corner — that
//! should be part of the sphere vertex blend. Without this
//! cleanup a 12-edge cube fillet reports 30 planes + 12 tori +
//! 8 spheres (= 50 faces) when the topologically correct answer
//! is 6 + 12 + 8 = 26 faces.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

/// Merge small triangular plane fragments into adjacent sphere
/// blend faces via `kef`.
///
/// For each plane face with ≤4 vertices, walks its outer loop
/// and kefs the first shared edge whose twin face is a sphere.
///
/// Note: `kef` keeps the half-edge's face and kills the twin's
/// face, so we pass in the *sphere-side* half-edge to preserve
/// the sphere surface. Otherwise the merged face would inherit
/// the plane's surface and lose the sphere labelling.
///
/// Iterates to a fixpoint: each successful kef may expose a new
/// fragment that can now merge into the enlarged blend.
pub(super) fn merge_corner_fragments_into_blends(brep: &mut BRep) {
    loop {
        let mut kef_target: Option<popkern_core::topology::HalfEdgeHandle> = None;
        'faces: for fh in brep.faces.keys() {
            let face = &brep.faces[fh];
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            if !matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
                continue;
            }
            let positions = brep.face_outer_loop_positions(fh);
            if positions.len() > 4 {
                continue;
            }
            let start_he = brep.loops[face.outer_loop].half_edge;
            let mut cursor = start_he;
            let max_iters = brep.half_edges.len() + 1;
            for _ in 0..max_iters {
                let he = &brep.half_edges[cursor];
                let twin_h = he.twin;
                let twin = &brep.half_edges[twin_h];
                let twin_face = brep.loops[twin.loop_].face;
                if twin_face != fh
                    && let Some(twin_sof) = brep.faces[twin_face].surface.as_ref()
                    && matches!(twin_sof.geometry, SurfaceGeometry::Sphere { .. })
                {
                    kef_target = Some(twin_h);
                    break 'faces;
                }
                cursor = he.next;
                if cursor == start_he {
                    break;
                }
            }
        }
        let Some(he_h) = kef_target else {
            break;
        };
        if brep.kef(he_h).is_err() {
            break;
        }
    }
}
