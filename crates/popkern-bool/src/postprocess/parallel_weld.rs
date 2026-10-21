//! Weld the two masked halves of a collapsed near-tangential gore
//! into one proper twin pair (PK-309).
//!
//! At a near-tangential junction interface the intersection-curve
//! connector and an operand boundary sliver bound a gore of material
//! thinner than the junction slop (~2e-3). The gore's region is
//! discarded, both its bounding halves are strand-masked
//! (self-twinned, curve stripped), and the served body carries two
//! single-sided rims 2e-3 apart — OCC free-bound micro wires and open
//! mesh segments (measured: rotated torus × torus, exactly two masked
//! strands bounding one collapsed gore).
//!
//! The weld: two SELF-TWINNED halves whose endpoints coincide
//! REVERSED within the junction class are the two sides of the
//! collapsed interface — wire them as twins and retire the duplicate
//! edge. Transactional per weld: structural validation must not
//! regress.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

/// Report from the parallel weld pass.
#[derive(Debug, Clone, Copy, Default)]
pub struct ParallelWeldReport {
    /// Pairs welded into proper twins.
    pub welded: usize,
    /// Candidate pairs declined (validation regressed) and rolled back.
    pub declined: usize,
}

/// Weld mutually-reversed self-twinned half-edge pairs.
pub fn weld_parallel_self_twins(brep: &mut BRep, tolerance: &Tolerance) -> ParallelWeldReport {
    weld_parallel_self_twins_in(brep, tolerance, &crate::serve_region::ServeRegion::Whole)
}

/// [`weld_parallel_self_twins`] for pairs whose edit stays in `region`:
/// both halves' faces and every face at the two vertices the weld
/// re-origins (a weld creates nothing, so the region stays exact).
pub(crate) fn weld_parallel_self_twins_in(
    brep: &mut BRep,
    tolerance: &Tolerance,
    region: &crate::serve_region::ServeRegion,
) -> ParallelWeldReport {
    let mut report = ParallelWeldReport::default();
    // Junction-slop class: the gore's two rims are up to ~2e-3 apart
    // mid-span; their shared endpoints coincide to the same class.
    let gate = tolerance.linear() * 1.0e4;
    loop {
        let selfies: Vec<(
            HalfEdgeHandle,
            popkern_core::geometry::Point3,
            popkern_core::geometry::Point3,
        )> = brep
            .half_edges
            .iter()
            .filter(|(heh, he)| he.twin == *heh)
            .filter_map(|(heh, he)| {
                let a = brep.vertices.get(he.origin)?.position;
                let b = brep
                    .half_edges
                    .get(he.next)
                    .and_then(|n| brep.vertices.get(n.origin))
                    .map(|v| v.position)?;
                Some((heh, a, b))
            })
            .collect();
        let mut pair: Option<(HalfEdgeHandle, HalfEdgeHandle)> = None;
        'outer: for (i, x) in selfies.iter().enumerate() {
            for y in selfies.iter().skip(i + 1) {
                if x.1.distance_to(&y.2) < gate
                    && x.2.distance_to(&y.1) < gate
                    && weld_in_region(brep, x.0, y.0, region)
                {
                    pair = Some((x.0, y.0));
                    break 'outer;
                }
            }
        }
        let Some((ha, hb)) = pair else { break };
        let snapshot = brep.clone();
        let edge_a = brep.half_edges[ha].edge;
        let edge_b = brep.half_edges[hb].edge;
        brep.half_edges[ha].twin = hb;
        brep.half_edges[hb].twin = ha;
        brep.half_edges[hb].edge = edge_a;
        // VERTEX UNIFICATION (PK-483): each masked rim kept its own
        // endpoint vertices, so the freshly wired twins straddle
        // duplicate vertices a junction-slop apart. Twins must share
        // endpoint vertices — and at a tangency pinch (the grazing
        // torus micro-lens) the unified vertex is what lets the
        // touching faces read as topologically adjacent: without it
        // validate_all flagged SelfIntersection on Boundary-only
        // witnesses at the two lens tips. Re-origin every half-edge
        // leaving the duplicate onto the partner's vertex; the
        // transactional validation guard below still rolls the whole
        // weld back if the rewire regresses structure.
        {
            let va0 = brep.half_edges[ha].origin;
            let va1 = brep.half_edges[brep.half_edges[ha].next].origin;
            let vb0 = brep.half_edges[hb].origin;
            let vb1 = brep.half_edges[brep.half_edges[hb].next].origin;
            for (kill, keep) in [(vb0, va1), (vb1, va0)] {
                if kill == keep || !brep.vertices.contains_key(kill) {
                    continue;
                }
                let leaving: Vec<HalfEdgeHandle> = brep
                    .half_edges
                    .iter()
                    .filter(|(_, d)| d.origin == kill)
                    .map(|(h, _)| h)
                    .collect();
                for h in leaving {
                    brep.half_edges[h].origin = keep;
                }
                brep.vertices.remove(kill);
            }
        }
        if edge_b != edge_a && !brep.half_edges.iter().any(|(_, he)| he.edge == edge_b) {
            brep.edges.remove(edge_b);
        }
        // The unified edge needs curve geometry (the masking stripped
        // it): restore from either original edge.
        if brep.edges.get(edge_a).is_some_and(|e| e.curve.is_none())
            && let Some(c) = snapshot
                .edges
                .get(edge_a)
                .and_then(|e| e.curve.clone())
                .or_else(|| snapshot.edges.get(edge_b).and_then(|e| e.curve.clone()))
            && let Some(e) = brep.edges.get_mut(edge_a)
        {
            e.curve = Some(c);
        }
        let before = snapshot.validate_collect_structural().len();
        let after = brep.validate_collect_structural().len();
        if after > before {
            *brep = snapshot;
            report.declined += 1;
            break;
        }
        report.welded += 1;
    }
    report
}

/// Whether welding `ha` with `hb` edits only faces `region` admits:
/// their own faces and every face at `hb`'s two ends (re-origined onto
/// `ha`'s).
fn weld_in_region(
    brep: &BRep,
    ha: HalfEdgeHandle,
    hb: HalfEdgeHandle,
    region: &crate::serve_region::ServeRegion,
) -> bool {
    let face = |h: HalfEdgeHandle| {
        brep.half_edges
            .get(h)
            .and_then(|he| brep.loops.get(he.loop_))
            .map(|l| l.face)
    };
    let ends = |h: HalfEdgeHandle| {
        let he = brep.half_edges.get(h)?;
        Some([he.origin, brep.half_edges.get(he.next)?.origin])
    };
    [ha, hb]
        .into_iter()
        .all(|h| face(h).is_some_and(|f| region.admits_face(f)))
        && ends(hb).is_some_and(|vs| vs.into_iter().all(|v| region.admits_vertex(brep, v)))
}
