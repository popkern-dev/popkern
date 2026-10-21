//! Removal of degenerate SELF-TWINNED stub half-edges (PK-448).
//!
//! `cleanup_unwired` deliberately self-twins the surviving half when
//! a fragment loop consumed only one side of an edge — for seam-style
//! survivors carrying real geometry (a boss/plate intersection
//! circle) that is honest containment. But a survivor with NO curve
//! and coincident endpoints is pure boolean residue riding in a live
//! loop: the mesher skips it (zero length), `validate` tolerates it,
//! and the STEP export faithfully writes it as a degenerate closed
//! FREE-BOUNDARY wire the OCC oracle rejects (measured: the
//! deep-cut 6.5-turn coil grind served 5 faces with one curveless
//! self-twin stub at a single vertex — `n_shells: 2`,
//! `free_bound_wires.closed: 1`, `passed: false`). Splice such stubs
//! out of their loops; keep every self-twin with real extent.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

/// Outcome of [`remove_self_twin_stubs`].
#[derive(Debug, Default)]
pub struct SelfTwinStubReport {
    /// Degenerate stubs spliced out and deleted.
    pub removed: usize,
    /// Self-twinned halves KEPT because they carry real extent
    /// (seam-style survivors).
    pub kept: usize,
    /// Degenerate stubs that could NOT be repaired (the splice broke
    /// the Euler ledger and rolled back — the body carries residue
    /// beyond the stub itself). A served body must not carry these:
    /// the caller declines the boolean instead.
    pub unrepaired: usize,
}

/// Splice degenerate self-twinned stub half-edges out of their loops
/// and delete them (with their edge record). Transactional: the body
/// is restored whole if the result fails validation.
pub fn remove_self_twin_stubs(brep: &mut BRep, tolerance: &Tolerance) -> SelfTwinStubReport {
    remove_self_twin_stubs_in(brep, tolerance, &crate::serve_region::ServeRegion::Whole)
}

/// [`remove_self_twin_stubs`] for the stubs in faces `region` admits
/// (a splice edits one loop and creates nothing).
pub(crate) fn remove_self_twin_stubs_in(
    brep: &mut BRep,
    tolerance: &Tolerance,
    region: &crate::serve_region::ServeRegion,
) -> SelfTwinStubReport {
    let mut report = SelfTwinStubReport::default();
    let candidates: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(h, he)| {
            he.twin == *h
                && brep
                    .loops
                    .get(he.loop_)
                    .is_some_and(|l| region.admits_face(l.face))
        })
        .map(|(h, _)| h)
        .collect();
    if candidates.is_empty() {
        return report;
    }
    let snapshot = brep.clone();
    let gate = tolerance.linear() * 10.0;
    for h in candidates {
        let Some(he) = brep.half_edges.get(h) else {
            continue;
        };
        let (next, prev, origin, loop_h, edge_h) = (he.next, he.prev, he.origin, he.loop_, he.edge);
        if next == h || prev == h {
            // Sole member of its loop — a different class (a whole
            // degenerate loop); leave it to the loop repairs.
            report.kept += 1;
            continue;
        }
        // Real-extent test: a genuine one-sided survivor carries a
        // curve with measurable extent, or spans distinct vertices.
        let end_origin = brep.half_edges[next].origin;
        let p_start = brep.vertices[origin].position;
        let p_end = brep.vertices[end_origin].position;
        let curve_extent = brep
            .edges
            .get(edge_h)
            .and_then(|e| e.curve.as_ref())
            .map(|c| {
                use popkern_core::geometry::curve::CurveEval;
                let a = c.geometry.point_at(c.t_start);
                (1..8)
                    .map(|i| {
                        let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                        c.geometry.point_at(t).distance_to(&a)
                    })
                    .fold(0.0, f64::max)
            })
            .unwrap_or(0.0);
        if p_start.distance_to(&p_end) > gate || curve_extent > gate {
            report.kept += 1;
            continue;
        }
        // Splice out of the loop.
        brep.half_edges[prev].next = next;
        brep.half_edges[next].prev = prev;
        if brep.loops.get(loop_h).map(|l| l.half_edge) == Some(h)
            && let Some(l) = brep.loops.get_mut(loop_h)
        {
            l.half_edge = next;
        }
        brep.remove_half_edge_reanchored(h);
        // The edge record is now unreferenced (the stub was its only
        // half) — drop it if nothing else claims it.
        let still_used = brep
            .half_edges
            .iter()
            .any(|(_, other)| other.edge == edge_h);
        if !still_used {
            brep.edges.remove(edge_h);
        }
        report.removed += 1;
    }
    if report.removed > 0
        && let Err(e) = brep.validate()
    {
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!("  [self-twin-stub] rollback: {e}");
        }
        *brep = snapshot;
        report.unrepaired = report.removed;
        report.removed = 0;
    }
    report
}
