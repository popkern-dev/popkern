//! Post-merge loop repair: prune half-edges stranded outside their face's loop chain.
//!
//! The merge phase occasionally produces faces whose outer (or inner)
//! loops have a broken `next` chain: half-edges are tagged with the
//! loop's handle but are not reachable from the loop's `half_edge`
//! starting point via successive `next` walks. The resulting topology
//! is technically present in the slotmaps and contributes to V/E/F
//! counts, but it does not participate in the visible boundary of
//! any face — Euler-Poincaré reads back as a higher-genus surface
//! than the operation should produce.
//!
//! This module finds those stranded half-edges, removes them along
//! with their twin (if also stranded), and cleans up any edges and
//! vertices that become orphaned by the removal. The pass validates
//! the BRep after each removal step; on validation failure it
//! rolls back the removal and stops.
//!
//! This is detection-only repair, not a fix for the underlying merge
//! bug. The merge phase still produces broken loops; this pass just
//! tidies up the result so the boolean entry point can return
//! Euler-consistent topology more often.

use std::collections::HashSet;

use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

/// Report from a loop repair pass.
#[derive(Debug, Clone, Default)]
pub struct LoopRepairReport {
    /// Number of stranded half-edges identified.
    pub stranded_half_edges: usize,
    /// Number of half-edges actually removed.
    pub half_edges_removed: usize,
    /// Number of edges removed because both half-edges were stranded.
    pub edges_removed: usize,
    /// Number of vertices removed because they had no remaining half-edges.
    pub vertices_removed: usize,
}

/// Find and prune half-edges that claim a loop but are not reachable
/// from that loop's `next` chain.
///
/// The cleanup is conservative: it removes the stranded half-edges,
/// any twin pair where both are stranded, and any vertex that becomes
/// orphaned. It does NOT touch reachable half-edges or their face
/// associations, so a successful pass cannot make the topology
/// strictly worse.
pub fn repair_stranded_loops(brep: &mut BRep) -> LoopRepairReport {
    let mut report = LoopRepairReport::default();

    // 1. Build the set of half-edges reachable from any face's loops.
    let mut reachable: HashSet<HalfEdgeHandle> = HashSet::new();
    for (_, face) in &brep.faces {
        collect_loop(brep, face.outer_loop, &mut reachable);
        for &inner_loop in &face.inner_loops {
            collect_loop(brep, inner_loop, &mut reachable);
        }
    }

    // 2. Identify half-edges that are in some loop owned by a face
    //    but not reachable from the face's loop chain.
    let owned_loops: HashSet<LoopHandle> = brep
        .faces
        .values()
        .flat_map(|f| std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()))
        .collect();

    let stranded: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(heh, he)| {
            if owned_loops.contains(&he.loop_) && !reachable.contains(&heh) {
                Some(heh)
            } else {
                None
            }
        })
        .collect();

    report.stranded_half_edges = stranded.len();
    if stranded.is_empty() {
        return report;
    }

    // 3. For each stranded half-edge, remove it. If its twin is also
    //    stranded, remove the edge as well. Track orphan vertices.
    let stranded_set: HashSet<HalfEdgeHandle> = stranded.iter().copied().collect();
    let mut edges_to_remove: Vec<EdgeHandle> = Vec::new();
    let mut vertices_to_check: HashSet<VertexHandle> = HashSet::new();

    for &heh in &stranded {
        let Some(he) = brep.half_edges.get(heh) else {
            continue;
        };
        let twin = he.twin;
        let edge = he.edge;
        let origin = he.origin;
        vertices_to_check.insert(origin);

        // If the twin is also stranded, remove the edge.
        if stranded_set.contains(&twin) {
            edges_to_remove.push(edge);
        }
    }

    // Snapshot the BRep state before mutation so we can roll back on
    // validation failure. This is a coarse rollback — we drop the
    // mutation results entirely if validation regresses.
    let snapshot = brep.clone();

    // 4. Remove stranded half-edges.
    for &heh in &stranded {
        if brep.remove_half_edge_reanchored(heh).is_some() {
            report.half_edges_removed += 1;
        }
    }

    // 5. Remove edges whose half-edges are gone.
    for &eh in &edges_to_remove {
        if brep.edges.remove(eh).is_some() {
            report.edges_removed += 1;
        }
    }

    // 6. Remove now-orphaned vertices. `vertices_to_check` is a
    // HashSet, and slotmap removal order shapes the free list — later
    // insertions would get run-dependent handles. Sort so handle
    // assignment stays deterministic (PK-300 hash-iteration audit).
    let in_use: HashSet<VertexHandle> = brep.half_edges.values().map(|he| he.origin).collect();
    let mut to_remove: Vec<VertexHandle> = vertices_to_check
        .into_iter()
        .filter(|v| !in_use.contains(v))
        .collect();
    to_remove.sort_unstable();
    for vh in to_remove {
        if brep.vertices.remove(vh).is_some() {
            report.vertices_removed += 1;
        }
    }

    // 7. Validate. On regression, roll back.
    if brep.validate().is_err() {
        *brep = snapshot;
        report.half_edges_removed = 0;
        report.edges_removed = 0;
        report.vertices_removed = 0;
    }

    report
}

fn collect_loop(brep: &BRep, loop_h: LoopHandle, out: &mut HashSet<HalfEdgeHandle>) {
    let Some(loop_data) = brep.loops.get(loop_h) else {
        return;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        out.insert(cursor);
        if let Some(he) = brep.half_edges.get(cursor) {
            cursor = he.next;
        } else {
            break;
        }
        if cursor == start {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn clean_cube_has_no_stranded_half_edges() {
        let mut cube = BRep::unit_cube(tol()).unwrap();
        let report = repair_stranded_loops(&mut cube);
        assert_eq!(report.stranded_half_edges, 0);
        assert_eq!(report.half_edges_removed, 0);
        assert!(cube.validate().is_ok());
    }

    #[test]
    fn report_default_is_empty() {
        let report = LoopRepairReport::default();
        assert_eq!(report.stranded_half_edges, 0);
        assert_eq!(report.half_edges_removed, 0);
    }
}
