//! Pinched vertices of a served boolean result (PK-992).
//!
//! The half-edges leaving a manifold vertex form ONE cycle under
//! `h -> twin(h).next`: the faces around it close up into a single
//! fan, the vertex's link is a disk. A vertex where they form two or
//! more cycles is shared by separate fans: two cones of material
//! touching at a point that the result carries as ONE vertex, which no
//! coincident-copy test sees. Measured: main − branch of an
//! equal-radius tee whose branch passes through the main cylinder
//! leaves the main body's material in two lobes, |z − 60| > |x|,
//! meeting only at the singular points (0, ±40, 60), and the cellular
//! route served that body as two-manifold edges around a pinched
//! vertex. Under decision D1 it declines with the vertex as a point
//! contact.

use popkern_core::topology::{BRep, ContactEntity, VertexNeighbourhood};

use crate::error::{ContactKind, NonManifoldContact};

/// Vertices of `brep` whose leaving half-edges form more than one
/// fan, in vertex order, each as a point contact. The fans are core's
/// (`BRep::vertex_neighbourhoods`, the one authority on manifoldness,
/// PK-1129): a walk from the middle of an open fan (a sheet's free rim)
/// counted it twice here before.
pub(crate) fn pinched_vertices(brep: &BRep) -> Vec<(NonManifoldContact, Vec<ContactEntity>)> {
    brep.vertex_neighbourhoods()
        .into_iter()
        .filter(|(_, n)| matches!(n, VertexNeighbourhood::Fans { .. }))
        .filter_map(|(v, _)| {
            let vertex = brep.vertices.get(v)?;
            Some((
                NonManifoldContact {
                    kind: ContactKind::Vertex,
                    location: vertex.position,
                },
                vec![ContactEntity::Vertex(v)],
            ))
        })
        .collect()
}
