//! Post-concatenation vertex and edge welding for sheet-sheet union.
//!
//! [`weld_coincident_vertices`] walks every vertex in a BRep, groups
//! them by 3D position within `tolerance.linear()`, and rewrites
//! every half-edge's `origin` to reference a single canonical
//! representative per group.
//!
//! [`weld_coincident_edges`] then walks edges whose endpoints have
//! collapsed to the same vertex pair, picks one as canonical,
//! redirects the other half-edges' `edge` back-reference to it,
//! splices the half-edges into a single non-manifold radial cycle
//! via [`popkern_core::topology::BRep::nonmanifold_radial`], and
//! deletes the duplicate edges. The result is a topologically
//! correct non-manifold BRep: one edge per geometric intersection
//! curve, with N half-edges threaded through the radial cycle so
//! downstream code that walks via [`BRep::radial_next`] visits every
//! face-use in order.
//!
//! Used by [`super::dispatch::boolean_with_sheets`] to weld two
//! imprinted sheets after `combine_breps` concatenates them.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};
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
    // The implementation lives in popkern_core::topology::brep_weld
    // as `BRep::weld_coincident_vertices`. This thin wrapper kept
    // the previous internal call sites compiling unchanged; new
    // call sites can use the inherent method directly.
    brep.weld_coincident_vertices(tolerance)
}

/// Weld edges whose endpoint vertex pairs are identical after
/// vertex welding into a single canonical edge, forming a
/// non-manifold radial cycle when multiple edges collapse together.
///
/// Algorithm:
///
/// 1. Build a map from `(v_min, v_max)` endpoint keys to the list
///    of edges connecting that pair. The key is independent of
///    direction so `E(v1 → v2)` and `E(v2 → v1)` group together.
///
/// 2. For each group with more than one edge, pick the first as
///    canonical. Collect all half-edges of the group (both via
///    `edge.half_edge.twin` AND the other edges' half-edges via
///    their `edge` back-reference), then splice them into a single
///    cyclic chain in [`BRep::nonmanifold_radial`]. The chain order
///    is the half-edge iteration order within each edge — for a
///    4-half-edge group, the cycle is
///    `he_A1 → he_A2 → he_B1 → he_B2 → he_A1`.
///
/// 3. Redirect every half-edge in the group to reference the
///    canonical edge via its `edge` field.
///
/// 4. Delete the non-canonical edges from the slotmap.
///
/// Only groups where every edge has matching curve geometry (or
/// all edges have `None` curves, i.e. implicit line segments
/// between their endpoints) are welded. Mismatched curves are left
/// untouched — welding them would silently lose geometry.
///
/// Returns the number of edges removed. Requires a preceding
/// [`weld_coincident_vertices`] call — edges are grouped by
/// canonical vertex handles, so pairs that look geometrically
/// coincident but reference different VertexHandles stay separate.
pub(super) fn weld_coincident_edges(brep: &mut BRep) -> usize {
    // Group edges by canonical (v_min, v_max) endpoint pair. The
    // VertexHandle ordering is via slotmap key comparison, which is
    // stable across snapshots of the same BRep.
    let mut groups: HashMap<(VertexHandle, VertexHandle), Vec<EdgeHandle>> = HashMap::new();
    for (eh, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let v1 = he.origin;
        let v2 = twin.origin;
        let key = if v1 <= v2 { (v1, v2) } else { (v2, v1) };
        groups.entry(key).or_default().push(eh);
    }

    let mut welds = 0;

    // Canonical group order: `groups` is a HashMap, and the losers'
    // `edges.remove()` calls below execute per group — hash order
    // would permute the edge slotmap's free list in the SERVED result
    // BRep, so later insertions would allocate run-varying handles
    // (PK-300 hash-iteration audit). Keys are canonical vertex pairs.
    let mut groups: Vec<_> = groups.into_iter().collect();
    groups.sort_unstable_by_key(|(key, _)| *key);
    for (_, edges) in groups {
        if edges.len() < 2 {
            continue;
        }

        // Verify every edge in the group has matching curve geometry.
        // Either all are None, or all share the same curve (by strict
        // equality — any drift is treated as "don't weld"). A more
        // tolerant comparison would need a curve-distance metric and
        // is deferred until a use case appears.
        let first_curve = brep.edges[edges[0]].curve.clone();
        let curves_match = edges.iter().all(|&eh| {
            let c = &brep.edges[eh].curve;
            match (&first_curve, c) {
                (None, None) => true,
                (Some(a), Some(b)) => a.t_start == b.t_start && a.t_end == b.t_end,
                _ => false,
            }
        });
        if !curves_match {
            continue;
        }

        let canonical = edges[0];

        // Collect all half-edges belonging to the group's edges via
        // the back-reference. We iterate `brep.half_edges` once and
        // pick up every half-edge whose `edge` field is in our set.
        let edge_set: std::collections::HashSet<EdgeHandle> = edges.iter().copied().collect();
        let mut group_half_edges: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter_map(|(heh, he)| {
                if edge_set.contains(&he.edge) {
                    Some(heh)
                } else {
                    None
                }
            })
            .collect();
        // Sort for deterministic cycle order — the slotmap iter
        // order is not guaranteed stable, but handle comparison is.
        group_half_edges.sort_by_key(|h| format!("{h:?}"));

        // Splice into a radial cycle in `nonmanifold_radial`. For
        // N half-edges h1..hN, set:
        //   radial_next[h1] = h2
        //   radial_next[h2] = h3
        //   ...
        //   radial_next[hN] = h1
        // This forms a closed cyclic chain. For N == 2 the cycle
        // collapses to the standard twin pair, but we still populate
        // the map for consistency — callers can always use
        // BRep::radial_next whether or not the edge is manifold.
        let n = group_half_edges.len();
        if n >= 2 {
            for i in 0..n {
                let curr = group_half_edges[i];
                let next = group_half_edges[(i + 1) % n];
                brep.nonmanifold_radial.insert(curr, next);
            }
        }

        // Redirect every half-edge in the group to the canonical
        // edge via its `edge` field.
        for heh in &group_half_edges {
            if let Some(he) = brep.half_edges.get_mut(*heh) {
                he.edge = canonical;
            }
        }

        // Delete the non-canonical edges.
        for &eh in &edges[1..] {
            brep.edges.remove(eh);
            welds += 1;
        }
    }

    welds
}

#[cfg(test)]
#[path = "weld_tests.rs"]
mod tests;
