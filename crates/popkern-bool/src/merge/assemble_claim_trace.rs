//! Claim-ledger diagnostics for graph-merge loop resolution
//! (PK-494) — split from `assemble.rs` under the 500-line budget.
//! Gated behind the `uv-scramble` trace topic plus
//! `POPKERN_CLAIM_TRACE_Z`: prints every resolution of a Shared
//! graph half-edge near the query height with both slot states, the
//! ledger that exposed the thread coil's end-region claim contest.

use std::collections::HashMap;

use popkern_core::topology::{HalfEdgeHandle, LoopHandle};

use super::super::graph::{FaceRef, IntersectionGraph};
use popkern_core::constructors::builder::BRepBuilder;

#[allow(clippy::too_many_arguments)]
pub(super) fn trace_claim(
    graph: &IntersectionGraph,
    face: &FaceRef,
    ghe: usize,
    graph_twin: usize,
    flip: bool,
    hmap: &HashMap<usize, HalfEdgeHandle>,
    bld: &BRepBuilder,
    default_loop: LoopHandle,
) {
    if !popkern_core::trace::trace_enabled("uv-scramble") {
        return;
    }
    let Some(zq) = std::env::var("POPKERN_CLAIM_TRACE_Z")
        .ok()
        .and_then(|v| v.parse::<f64>().ok())
    else {
        return;
    };
    let o = graph.nodes[graph.half_edges[ghe].origin].position;
    let e = graph.nodes[graph.half_edges[graph_twin].origin].position;
    if (o.z - zq).abs() >= 0.5 && (e.z - zq).abs() >= 0.5 {
        return;
    }
    let state = |g: usize| match hmap.get(&g) {
        None => "unmapped".to_string(),
        Some(&bh) => {
            if bld.brep().half_edges[bh].loop_ == default_loop {
                "free".to_string()
            } else {
                "taken".to_string()
            }
        }
    };
    eprintln!(
        "  [claim] ghe {ghe} src {:?} flip {flip} for {face:?} span ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) slots self:{} twin:{}",
        graph.half_edges[ghe].source,
        o.x,
        o.y,
        o.z,
        e.x,
        e.y,
        e.z,
        state(ghe),
        state(graph_twin)
    );
}

/// The alt-install ledger line: which fragment takes a twin-direction
/// half-edge, who owns the primary, and the edge span.
pub(super) fn trace_alt_install(
    graph: &IntersectionGraph,
    face: &FaceRef,
    ghe: usize,
    graph_twin: usize,
    hmap: &HashMap<usize, HalfEdgeHandle>,
    bld: &BRepBuilder,
) {
    if !popkern_core::trace::trace_enabled("uv-scramble") {
        return;
    }
    let owner = hmap
        .get(&ghe)
        .map(|&bh| {
            let lp = bld.brep().half_edges[bh].loop_;
            format!("{:?}", bld.brep().loops.get(lp).map(|l| l.face))
        })
        .unwrap_or_else(|| "unmapped".into());
    let (o, e) = (
        graph.nodes[graph.half_edges[ghe].origin].position,
        graph.nodes[graph.half_edges[graph_twin].origin].position,
    );
    eprintln!(
        "  [uv-scramble] alt-install ghe {ghe} (twin {graph_twin}) src {:?} for {face:?} primary-owner {owner} span ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
        graph.half_edges[ghe].source, o.x, o.y, o.z, e.x, e.y, e.z
    );
}

/// The dup-mint ledger line (`boolean` topic): both directions of a
/// graph edge are wired, so a fresh brep pair is minted — with the
/// owners of the consumed slots named.
pub(super) fn trace_dup_mint(
    graph: &IntersectionGraph,
    face: &FaceRef,
    ghe: usize,
    graph_twin: usize,
    hmap: &HashMap<usize, HalfEdgeHandle>,
    bld: &BRepBuilder,
) {
    if !popkern_core::trace::trace_enabled("boolean") {
        return;
    }
    let owner = |g: usize| {
        hmap.get(&g)
            .map(|&bh| {
                let lp = bld.brep().half_edges[bh].loop_;
                bld.brep().loops.get(lp).map(|l| l.face)
            })
            .map(|f| format!("{f:?}"))
            .unwrap_or_else(|| "unmapped".into())
    };
    eprintln!(
        "  [assemble-dup] ghe {ghe} (twin {graph_twin}) for {face:?} src {:?} taken-by {} / {}",
        graph.half_edges[ghe].source,
        owner(ghe),
        owner(graph_twin)
    );
}
