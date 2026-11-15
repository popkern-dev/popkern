//! An opening rim vertex reseeded inside the opening (PK-1216).
//!
//! A rim vertex is solved on the kept faces' offsets and the removed
//! faces' own surfaces, seeded from where it stood. Where those surfaces
//! meet the seed square to both (a capsule opened at its dome: the
//! offset wall, radius 4, and the dome's sphere, radius 5 about the
//! wall's top, are both radial at the old rim), projection cycles
//! between the two (onto the wall at (4, 0, 10), back onto the sphere
//! at (5, 0, 10)) and Newton sees dependent normals: the vertex kept its
//! old place and the wall opened by √2. The rim lies inside the removed
//! face, where the cavity runs into it (the circle of radius 4 at
//! z = 13): seeds stepped from the vertex toward the removed faces'
//! other vertices reach it where the surfaces cross transversally.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Fractions of the way from the vertex toward each removed face's
/// other vertices at which to reseed.
const STEPS: [f64; 4] = [0.25, 0.5, 0.75, 0.1];

/// The rim vertex `vertex` of `outer` solved on `surfaces` from seeds
/// inside `removed` faces; `None` when none converges.
pub(crate) fn reseed_into_openings(
    outer: &BRep,
    vertex: VertexHandle,
    removed: &[FaceHandle],
    surfaces: &[&SurfaceGeometry],
    tolerance: &Tolerance,
) -> Option<Point3> {
    let at = outer.vertices.get(vertex)?.position;
    let mut targets: Vec<Point3> = Vec::new();
    for &face in removed {
        for p in outer.face_outer_loop_positions(face) {
            if p.distance_to(&at) > tolerance.linear()
                && !targets
                    .iter()
                    .any(|q| q.distance_to(&p) <= tolerance.linear())
            {
                targets.push(p);
            }
        }
    }
    STEPS.iter().find_map(|&f| {
        targets.iter().find_map(|&t| {
            let seed = at + (t - at) * f;
            crate::curved_recompute::solve_vertex_on_surfaces(surfaces, seed, tolerance)
        })
    })
}
