//! Unit gates for the support-loop interval rewrite (PK-167 step 1).
//!
//! Fixture: unit cube with one side face retired as construction
//! residue — its four edges keep single boundary sides on the
//! neighbouring faces, which is exactly the production shape (rails
//! and cap edges expose twinless boundary sides after the old
//! support material is removed). The rewrite splices the top face's
//! twinless side out and twins the three remaining boundary sides
//! into the gap: the top loop detours around the opening, face
//! identity and winding preserved, and the shell closes with every
//! edge two-sided again.

use super::{SupportLoopIntervalRewritePlan, SupportLoopIntervalRewriteResult};
use crate::error::CoreError;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, TopologyRole, VertexHandle};

fn cube_with_retired_side() -> (BRep, FaceHandle) {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let side = brep
        .faces
        .keys()
        .nth(2)
        .expect("cube has at least three faces");
    brep.faces[side].role = TopologyRole::ConstructionResidual;
    brep.retire_construction_residual_faces(&[side])
        .expect("retire side face");
    (brep, side)
}

/// The single boundary side of `edge` (its surviving half-edge).
fn boundary_side(brep: &BRep, edge: crate::topology::EdgeHandle) -> HalfEdgeHandle {
    let h = brep.edges[edge].half_edge;
    assert!(
        !brep.half_edges.contains_key(brep.half_edges[h].twin),
        "expected a twinless boundary side"
    );
    h
}

fn destination(brep: &BRep, h: HalfEdgeHandle) -> VertexHandle {
    let next = brep.half_edges[h].next;
    brep.half_edges[next].origin
}

/// Build the rewrite plan: removed = one twinless side of some loop
/// around the opening; replacement = the other three boundary
/// sides, chained by vertex identity per the operator contract.
fn opening_rewrite_plan(brep: &BRep) -> SupportLoopIntervalRewritePlan {
    let boundary: Vec<_> = brep.boundary_edges.iter().copied().collect();
    assert_eq!(boundary.len(), 4, "retired side leaves 4 boundary edges");
    // Choose the removed side: any boundary side; its loop hosts the
    // rewrite. The other three sides fill the gap.
    let removed_side = boundary_side(brep, boundary[0]);
    let start = brep.half_edges[removed_side].origin;
    let end = destination(brep, removed_side);
    let mut rest: Vec<HalfEdgeHandle> = boundary[1..]
        .iter()
        .map(|&e| boundary_side(brep, e))
        .collect();
    // Chain: destination(replacement[0]) == start, then
    // origin(r[i]) == destination(r[i+1]), ending origin == end.
    let mut chain: Vec<HalfEdgeHandle> = Vec::new();
    let mut want_dest = start;
    while !rest.is_empty() {
        let pos = rest
            .iter()
            .position(|&h| destination(brep, h) == want_dest)
            .expect("boundary sides chain around the opening");
        let h = rest.remove(pos);
        want_dest = brep.half_edges[h].origin;
        chain.push(h);
    }
    assert_eq!(want_dest, end, "chain ends at the interval end vertex");
    SupportLoopIntervalRewritePlan {
        loop_: brep.half_edges[removed_side].loop_,
        removed: vec![removed_side],
        replacement: chain,
        pair_removed_edges: false,
    }
}

#[test]
fn rewrite_closes_the_opening_with_dv_df_zero() {
    let (mut brep, _) = cube_with_retired_side();
    let (v0, e0, f0, l0) = (
        brep.vertices.len(),
        brep.edges.len(),
        brep.faces.len(),
        brep.loops.len(),
    );
    let plan = opening_rewrite_plan(&brep);
    let host_loop = plan.loop_;
    let result: SupportLoopIntervalRewriteResult = brep
        .support_loop_interval_rewrite(&plan)
        .expect("rewrite commits");
    assert_eq!(result.created.len(), 3, "three sides twin into the gap");
    assert_eq!(
        result.retired_edges.len(),
        1,
        "the removed twinless side retires its edge"
    );
    assert!(result.exposed_edges.is_empty());
    assert_eq!(brep.vertices.len(), v0, "dV = 0");
    assert_eq!(brep.faces.len(), f0, "dF = 0");
    assert_eq!(brep.loops.len(), l0, "dL = 0");
    assert_eq!(brep.edges.len(), e0 - 1, "dE = -1 (retired edge)");
    assert!(
        brep.boundary_edges.is_empty(),
        "every edge is two-sided again"
    );
    // The host loop closes through the created chain with winding
    // preserved (walk sees each created side exactly once).
    let start = brep.loops[host_loop].half_edge;
    let mut cursor = start;
    let mut created_seen = 0;
    for _ in 0..brep.half_edges.len() + 1 {
        if result.created.contains(&cursor) {
            created_seen += 1;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    assert_eq!(created_seen, 3, "host loop contains the full new chain");
    // Euler-Poincare on the closed result.
    let chi = brep.vertices.len() as i64 - brep.edges.len() as i64 + brep.faces.len() as i64;
    assert_eq!(chi, 2, "closed shell Euler characteristic");
}

#[test]
fn rewrite_declines_replacement_with_live_twin() {
    let (mut brep, _) = cube_with_retired_side();
    let mut plan = opening_rewrite_plan(&brep);
    // Swap one replacement for a fully-twinned interior side.
    let twinned = brep
        .half_edges
        .iter()
        .find(|(h, data)| brep.half_edges.contains_key(data.twin) && !plan.removed.contains(h))
        .map(|(h, _)| h)
        .expect("cube has twinned sides");
    plan.replacement[0] = twinned;
    let err = brep
        .support_loop_interval_rewrite(&plan)
        .expect_err("live-twin replacement must decline");
    assert!(matches!(err, CoreError::PreconditionFailed(_)));
}

#[test]
fn rewrite_declines_broken_chain() {
    let (mut brep, _) = cube_with_retired_side();
    let mut plan = opening_rewrite_plan(&brep);
    plan.replacement.reverse();
    let err = brep
        .support_loop_interval_rewrite(&plan)
        .expect_err("mis-ordered chain must decline");
    assert!(matches!(err, CoreError::PreconditionFailed(_)));
}

#[test]
fn rewrite_declines_whole_loop_removal() {
    let (mut brep, _) = cube_with_retired_side();
    let base = opening_rewrite_plan(&brep);
    // Extend `removed` to the entire host loop.
    let start = base.removed[0];
    let mut removed = vec![start];
    let mut cursor = brep.half_edges[start].next;
    while cursor != start {
        removed.push(cursor);
        cursor = brep.half_edges[cursor].next;
    }
    let plan = SupportLoopIntervalRewritePlan {
        loop_: base.loop_,
        removed,
        replacement: base.replacement,
        pair_removed_edges: false,
    };
    let err = brep
        .support_loop_interval_rewrite(&plan)
        .expect_err("whole-loop removal must decline");
    assert!(matches!(err, CoreError::PreconditionFailed(_)));
}
