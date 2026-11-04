//! The pieces of a body (design: operations on the cell complex, 5.4):
//! its solids, each kept or dropped whole with the faces, edges and
//! vertices that bound it, their names and lineage unchanged.
//!
//! A split serves its pieces as one body of several solids
//! ([`crate::split`], [`crate::concatenate_solids`]), each with its own
//! topology: pieces meeting at a contact hold one copy of the contact
//! edge or vertex each, the contact recorded on the body. Taking
//! pieces out therefore removes whole components; nothing is cut. A
//! body whose solids share an entity (a vertex or edge used by two
//! solids) is not pieces of this kind and declines.

use std::collections::HashSet;

use popkern_core::topology::{BRep, ContactEntity, SolidHandle};

use crate::error::BoolError;

/// The body's pieces, its solids, in slot order (the order a split
/// served them in).
pub fn piece_handles(body: &BRep) -> Vec<SolidHandle> {
    body.solids.keys().collect()
}

/// `body` with only the solids in `keep`: every entity bounding a
/// dropped solid is removed, every kept entity keeps its handle, name
/// and geometry, and a recorded contact survives when every entity it
/// names is kept.
///
/// ```
/// use popkern_bool::{concatenate_solids, pieces};
/// use popkern_core::{geometry::Transform3, tolerance::Tolerance, topology::BRep};
///
/// let tol = Tolerance::default();
/// let a = BRep::unit_cube(tol).unwrap();
/// let b = BRep::unit_cube(tol)
///     .unwrap()
///     .transform(&Transform3::translation(3.0, 0.0, 0.0))
///     .unwrap();
/// let both = concatenate_solids(a, b).unwrap();
/// let first = pieces::piece_handles(&both)[0];
/// let one = pieces::retain_pieces(&both, &[first]).unwrap();
/// assert_eq!((one.solids.len(), one.faces.len()), (1, 6));
/// ```
///
/// # Errors
///
/// [`BoolError::DegenerateInput`] when a kept and a dropped solid share
/// a half-edge, edge or vertex (the body is not separable pieces).
pub fn retain_pieces(body: &BRep, keep: &[SolidHandle]) -> Result<BRep, BoolError> {
    let solids: HashSet<_> = keep
        .iter()
        .copied()
        .filter(|s| body.solids.contains_key(*s))
        .collect();
    let shells: HashSet<_> = body
        .shells
        .iter()
        .filter(|(_, s)| solids.contains(&s.solid))
        .map(|(h, _)| h)
        .collect();
    let faces: HashSet<_> = body
        .faces
        .iter()
        .filter(|(_, f)| shells.contains(&f.shell))
        .map(|(h, _)| h)
        .collect();
    let loops: HashSet<_> = body
        .loops
        .iter()
        .filter(|(_, l)| faces.contains(&l.face))
        .map(|(h, _)| h)
        .collect();
    let half_edges: HashSet<_> = body
        .half_edges
        .iter()
        .filter(|(_, he)| loops.contains(&he.loop_))
        .map(|(h, _)| h)
        .collect();
    let edges: HashSet<_> = half_edges
        .iter()
        .filter_map(|h| body.half_edges.get(*h).map(|he| he.edge))
        .collect();
    let vertices: HashSet<_> = half_edges
        .iter()
        .filter_map(|h| body.half_edges.get(*h).map(|he| he.origin))
        .collect();
    let shared = |what: &str| {
        BoolError::DegenerateInput(format!(
            "the body's solids share {what}: its pieces cannot be taken apart"
        ))
    };
    for h in &half_edges {
        let he = &body.half_edges[*h];
        let links = [he.twin, he.next, he.prev];
        if links
            .iter()
            .any(|l| body.half_edges.contains_key(*l) && !half_edges.contains(l))
        {
            return Err(shared("a half-edge"));
        }
    }
    for (h, he) in &body.half_edges {
        if half_edges.contains(&h) {
            continue;
        }
        if edges.contains(&he.edge) {
            return Err(shared("an edge"));
        }
        if vertices.contains(&he.origin) {
            return Err(shared("a vertex"));
        }
    }

    let mut out = body.clone();
    out.solids.retain(|h, _| solids.contains(&h));
    out.shells.retain(|h, _| shells.contains(&h));
    out.faces.retain(|h, _| faces.contains(&h));
    out.loops.retain(|h, _| loops.contains(&h));
    out.half_edges.retain(|h, _| half_edges.contains(&h));
    out.edges.retain(|h, _| edges.contains(&h));
    out.vertices.retain(|h, _| vertices.contains(&h));
    out.boundary_edges.retain(|e| edges.contains(e));
    out.nonmanifold_radial
        .retain(|a, b| half_edges.contains(a) && half_edges.contains(b));
    out.surface_boundary_hints
        .retain(|h, _| half_edges.contains(h));
    out.surface_boundary_declines
        .retain(|h, _| half_edges.contains(h));
    out.surface_boundary_audit_events
        .retain(|e| half_edges.contains(&e.half_edge));
    out.face_use_pcurves
        .retain(|k, _| faces.contains(&k.face) && half_edges.contains(&k.half_edge));
    out.topology_contacts.retain(|c| {
        c.entities.iter().all(|e| match e {
            ContactEntity::Vertex(v) => vertices.contains(v),
            ContactEntity::Edge(e) => edges.contains(e),
            ContactEntity::Face(f) => faces.contains(f),
        })
    });
    Ok(out)
}

/// Each piece of `body` as its own body (see [`retain_pieces`]), in
/// [`piece_handles`] order.
///
/// # Errors
///
/// As [`retain_pieces`].
pub fn piece_bodies(body: &BRep) -> Result<Vec<BRep>, BoolError> {
    piece_handles(body)
        .into_iter()
        .map(|s| retain_pieces(body, &[s]))
        .collect()
}

#[cfg(test)]
#[path = "pieces_tests.rs"]
mod tests;
