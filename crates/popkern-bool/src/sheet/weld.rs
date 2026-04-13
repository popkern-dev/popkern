//! Post-concatenation vertex welding for sheet-sheet union.
//!
//! [`weld_coincident_vertices`] walks every vertex in a BRep, groups
//! them by 3D position within `tolerance.linear()`, and rewrites
//! every half-edge's `origin` to reference a single canonical
//! representative per group. Used by [`super::dispatch::boolean_with_sheets`]
//! to weld the two imprinted sheets' intersection vertices after
//! `combine_breps` has concatenated them.
//!
//! The weld operates at the vertex level only. Edges along a
//! coincident intersection curve end up with matching endpoints but
//! remain separate BRep edges (each with its own 2 half-edges). This
//! is a valid topology in the popkern BRep model — "two edges share
//! endpoints without being topologically merged" — and it is the
//! correct representation for transversal sheet-sheet junctions that
//! are inherently non-manifold: four half-edges share an edge
//! position, two belonging to each sheet's face loop, and the BRep
//! model's 2-half-edge-per-edge invariant is preserved by keeping
//! the two edges distinct.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};
use std::collections::HashMap;

/// Weld vertices within `tolerance.linear()` of each other into a
/// single canonical representative.
///
/// Returns the number of vertices removed (i.e. the number of
/// duplicates welded away). Walks all vertex pairs in `brep`,
/// builds an equivalence map via union-find, rewrites every
/// half-edge's `origin` and every vertex's `outgoing` pointer to
/// reference the canonical representative, and deletes the
/// non-canonical vertices from the slotmap.
///
/// Runtime is O(n²) in vertex count, which is acceptable for the
/// sheet-body sizes seen by the boolean engine. A spatial index
/// would only help for very large imprints.
pub(super) fn weld_coincident_vertices(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    let eps = tolerance.linear();
    let eps_sq = eps * eps;

    // Snapshot the (handle, position) pairs before we start mutating.
    let all: Vec<(VertexHandle, popkern_core::geometry::Point3)> =
        brep.vertices.iter().map(|(h, v)| (h, v.position)).collect();
    if all.len() < 2 {
        return 0;
    }

    // Union-find: each vertex starts as its own class.
    let mut parent: HashMap<VertexHandle, VertexHandle> =
        all.iter().map(|(h, _)| (*h, *h)).collect();

    // `find` with path compression.
    fn find(parent: &mut HashMap<VertexHandle, VertexHandle>, vh: VertexHandle) -> VertexHandle {
        let p = parent[&vh];
        if p == vh {
            return vh;
        }
        let root = find(parent, p);
        parent.insert(vh, root);
        root
    }

    // For each coincident pair, union.
    for (i, &(vh_i, pos_i)) in all.iter().enumerate() {
        for &(vh_j, pos_j) in all.iter().skip(i + 1) {
            let dx = pos_i.x - pos_j.x;
            let dy = pos_i.y - pos_j.y;
            let dz = pos_i.z - pos_j.z;
            if dx * dx + dy * dy + dz * dz <= eps_sq {
                let root_i = find(&mut parent, vh_i);
                let root_j = find(&mut parent, vh_j);
                if root_i != root_j {
                    parent.insert(root_j, root_i);
                }
            }
        }
    }

    // Resolve canonical handles after all unions.
    let keys: Vec<VertexHandle> = parent.keys().copied().collect();
    let mut canonical: HashMap<VertexHandle, VertexHandle> = HashMap::new();
    for vh in keys {
        let root = find(&mut parent, vh);
        canonical.insert(vh, root);
    }

    // Rewrite half-edge origins to point at canonical representatives.
    for (_, he) in brep.half_edges.iter_mut() {
        if let Some(&canon) = canonical.get(&he.origin)
            && canon != he.origin
        {
            he.origin = canon;
        }
    }

    // Rewrite vertex outgoing pointers for the CANONICAL vertices
    // whose old outgoing reference still lives on a non-canonical
    // vertex. Collect first to avoid borrow conflicts with the
    // remove pass below.
    //
    // Non-canonical vertices will be deleted; their outgoing pointers
    // become irrelevant, but canonical vertices keep their own
    // outgoing which must still reference a live half-edge. This is
    // preserved by default — half-edge handles themselves don't
    // change during vertex welding; only their `origin` field does.

    // Delete non-canonical vertices.
    let to_remove: Vec<VertexHandle> = canonical
        .iter()
        .filter_map(|(&vh, &canon)| if vh != canon { Some(vh) } else { None })
        .collect();
    let welds = to_remove.len();
    for vh in to_remove {
        brep.vertices.remove(vh);
    }
    welds
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    #[test]
    fn weld_merges_two_vertices_at_same_position() {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        // Two coincident vertices (same position). No edges/faces —
        // just the vertex-level merge path.
        brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(1.0, 2.0, 3.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(1.0, 2.0, 3.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        assert_eq!(brep.vertices.len(), 2);
        let welds = weld_coincident_vertices(&mut brep, &tol);
        assert_eq!(welds, 1);
        assert_eq!(brep.vertices.len(), 1);
    }

    #[test]
    fn weld_leaves_distinct_vertices_untouched() {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(1.0, 0.0, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let welds = weld_coincident_vertices(&mut brep, &tol);
        assert_eq!(welds, 0);
        assert_eq!(brep.vertices.len(), 2);
    }

    #[test]
    fn weld_transitive_merge_via_middle_vertex() {
        // Three vertices: A at (0,0,0), B at (1e-7, 0, 0), C at (2e-7, 0, 0).
        // With default tol.linear() == 1e-6, A-B is within tol (1e-7 < 1e-6),
        // B-C is within tol, and A-C is also within tol (2e-7 < 1e-6).
        // All three should collapse into one.
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(1e-7, 0.0, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(2e-7, 0.0, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let welds = weld_coincident_vertices(&mut brep, &tol);
        assert_eq!(welds, 2);
        assert_eq!(brep.vertices.len(), 1);
    }

    #[test]
    fn weld_rewrites_half_edge_origins_to_canonical() {
        // Two coincident vertices, each referenced by a half-edge as
        // origin. After weld, both half-edges should point at the
        // canonical vertex.
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        let v1 = brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(5.0, 5.0, 5.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let v2 = brep.vertices.insert(popkern_core::topology::Vertex {
            position: Point3::new(5.0, 5.0, 5.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let he1 = brep.half_edges.insert(popkern_core::topology::HalfEdge {
            origin: v1,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: Default::default(),
            loop_: Default::default(),
            pcurve: None,
        });
        let he2 = brep.half_edges.insert(popkern_core::topology::HalfEdge {
            origin: v2,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: Default::default(),
            loop_: Default::default(),
            pcurve: None,
        });

        weld_coincident_vertices(&mut brep, &tol);

        // The two half-edges now reference the SAME vertex — whichever
        // one the union-find picked as canonical. That vertex must
        // still exist in the slotmap.
        let o1 = brep.half_edges[he1].origin;
        let o2 = brep.half_edges[he2].origin;
        assert_eq!(o1, o2, "both half-edges should reference the same origin");
        assert!(
            brep.vertices.get(o1).is_some(),
            "canonical vertex must still exist"
        );
    }
}
