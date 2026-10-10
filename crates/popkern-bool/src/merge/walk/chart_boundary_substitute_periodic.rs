//! Boundary-substitute direction resolution on PERIODIC charts
//! (PK-750), the chart-free counterpart of `chart_boundary_substitute`.
//!
//! A coincident boolean can register an overlap record ON a periodic
//! face's own boundary: a box whose caps are flush with a rod's caps
//! cuts the rod's rim circles, and the near rim arc is a Shared record
//! (rod wall × box cap) while the operand's own arc there is removed as
//! its duplicate. Both record directions register on the wall, but only
//! one bounds it; the other's material side is the cap's. Left in, the
//! surplus direction unbalanced the rim junction nodes (2 in / 3 out),
//! the seam-config search declined the chart, and the wall walked on
//! the 3D-Newell fallback.
//!
//! The non-periodic resolver decides this class on UV polygons of the
//! original loops; on a periodic chart those polygons and the record
//! runs live on different branches. This resolver decides the same two
//! facts in 3D, where no branch exists:
//!   1. A Shared pair is a BOUNDARY SUBSTITUTE iff every interior
//!      sample of its curve lies on a piece of the face's ORIGINAL
//!      boundary (outer or inner loop, sampled on the edge's own
//!      parameter window) and not on an edge still registered on the
//!      face (that partial-cover class belongs to the duplicate-cover
//!      dedup).
//!   2. The kept direction is the one whose travel agrees with the
//!      original loop's traversal there (material-left by the BRep
//!      winding invariant). A closed edge's traversal per side is the
//!      PK-303 canonical-forward convention: the canonical half-edge
//!      runs its curve in increasing t.
//!
//! **Tolerance:** 100 × linear for "on the boundary", the same local
//! band as the blend's on-surface acceptance; a genuine cut that stays
//! within it over its whole run is a sliver no walk can separate.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::chart_subdivision::{FaceChartSubdivision, FilteredRecord};
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph};
use crate::merge::graph_duplicates::closest_on_piece;

/// Interior samples per record run.
const RUN_SAMPLES: usize = 7;

/// A piece of the original boundary: its curve and whether the loop
/// traverses it in increasing t.
struct BoundaryPiece<'a> {
    curve: &'a CurveOnEdge,
    forward: bool,
    seam: bool,
}

/// The original face's boundary pieces with their traversal sense.
fn boundary_pieces<'a>(brep: &'a BRep, face: &FaceRef, linear: f64) -> Vec<BoundaryPiece<'a>> {
    let Some(f) = brep.faces.get(face.face) else {
        return Vec::new();
    };
    let mut pieces = Vec::new();
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(edge) = brep.edges.get(h.edge)
                && let Some(curve) = edge.curve.as_ref()
            {
                let origin = brep.vertices.get(h.origin).map(|v| v.position);
                let a = curve.geometry.point_at(curve.t_start);
                let b = curve.geometry.point_at(curve.t_end);
                // A closed edge (or a missing vertex) follows the
                // canonical-forward convention (the canonical half runs
                // in increasing t); an open edge's sense is its
                // endpoints' against the stored window's order (a
                // reversed arc stores t_start > t_end).
                let forward = match origin {
                    Some(o) if a.distance_to(&b) > linear => {
                        (a.distance_to(&o) <= b.distance_to(&o)) == (curve.t_end >= curve.t_start)
                    }
                    _ => edge.half_edge == he,
                };
                // A SEAM edge bounds the face in both directions (one
                // per chart border): a record on it keeps both, and
                // the seam machinery owns it.
                let seam = brep
                    .half_edges
                    .get(h.twin)
                    .and_then(|tw| brep.loops.get(tw.loop_))
                    .is_some_and(|l| l.face == face.face);
                pieces.push(BoundaryPiece {
                    curve,
                    forward,
                    seam,
                });
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    pieces
}

/// The 3D curve of graph half-edge `h` and whether `h` runs it in
/// increasing t.
fn record_curve(graph: &IntersectionGraph, h: usize) -> Option<(&CurveOnEdge, bool)> {
    let he = &graph.half_edges[h];
    let curve = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())?;
    let origin = graph.nodes[he.origin].position;
    let dest = graph.nodes[graph.half_edges[he.twin].origin].position;
    let start = curve.geometry.point_at(curve.t_start);
    let from_start = start.distance_to(&origin) <= start.distance_to(&dest);
    Some((curve, from_start == (curve.t_end >= curve.t_start)))
}

/// Resolve boundary-substitute pair directions on a periodic `sub` in
/// place (see the module docs).
pub(super) fn resolve_periodic_boundary_substitutes(
    graph: &IntersectionGraph,
    face: &FaceRef,
    sub: &mut FaceChartSubdivision,
    original_brep: &BRep,
    tolerance: &Tolerance,
) {
    let pairs: Vec<(usize, usize)> = sub
        .embeddings
        .keys()
        .copied()
        .filter(|&h| {
            let tw = graph.half_edges[h].twin;
            h < tw
                && sub.embeddings.contains_key(&tw)
                && matches!(graph.half_edges[h].source, EdgeSource::Shared { .. })
        })
        .map(|h| (h, graph.half_edges[h].twin))
        .collect();
    if pairs.is_empty() {
        return;
    }
    let pieces = boundary_pieces(original_brep, face, tolerance.linear());
    if pieces.is_empty() {
        return;
    }
    let registered: Vec<&CurveOnEdge> = sub
        .embeddings
        .keys()
        .filter(|&&h| !matches!(graph.half_edges[h].source, EdgeSource::Shared { .. }))
        .filter_map(|&h| record_curve(graph, h).map(|(c, _)| c))
        .collect();
    let band = tolerance.linear() * 100.0;
    let trace = popkern_core::trace::trace_enabled("boundary-substitute");

    let mut resolved: Vec<(usize, usize)> = Vec::new(); // (kept, dropped)
    for (h, tw) in pairs {
        let Some((curve, h_forward)) = record_curve(graph, h) else {
            continue;
        };
        let (mut on_boundary, mut covered, mut agree) = (true, true, 0i64);
        for k in 1..=RUN_SAMPLES {
            let t =
                curve.t_start + (curve.t_end - curve.t_start) * k as f64 / (RUN_SAMPLES + 1) as f64;
            let p = curve.geometry.point_at(t);
            covered &= registered.iter().any(|c| closest_on_piece(c, &p).1 <= band);
            let nearest = pieces
                .iter()
                .map(|piece| (piece, closest_on_piece(piece.curve, &p)))
                .min_by(|a, b| a.1.1.total_cmp(&b.1.1));
            let Some((piece, (tp, _))) = nearest.filter(|(p, (_, d))| *d <= band && !p.seam) else {
                on_boundary = false;
                break;
            };
            let run = curve.geometry.tangent_at(t) * if h_forward { 1.0 } else { -1.0 };
            let along =
                piece.curve.geometry.tangent_at(tp) * if piece.forward { 1.0 } else { -1.0 };
            agree += run.dot(&along).signum() as i64;
        }
        if trace {
            eprintln!(
                "      [bsub-periodic] {face:?} pair ({h},{tw}): on_boundary={on_boundary} covered={covered} agree={agree}"
            );
        }
        if !on_boundary || covered {
            continue;
        }
        match agree.cmp(&0) {
            std::cmp::Ordering::Greater => resolved.push((h, tw)),
            std::cmp::Ordering::Less => resolved.push((tw, h)),
            // Indecisive: leave both registered; the balance search
            // and the 3D-Newell walk stay its recorded fallback.
            std::cmp::Ordering::Equal => {}
        }
    }
    for (kept, dropped) in resolved {
        if sub.embeddings.remove(&dropped).is_some() {
            sub.report
                .filtered
                .push(FilteredRecord::BoundarySubstituteDirection { dropped, kept });
        }
    }
}
