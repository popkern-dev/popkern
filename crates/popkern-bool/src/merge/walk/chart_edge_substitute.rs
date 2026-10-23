//! Edge-coincident boundary substitutes (PK-920).
//!
//! A record can lie ON an edge of the face it is registered on without
//! the two operands having a coincident face pair: a disc's rim circle
//! lying on a coaxial bore wall of the same radius. The graph then
//! carries the rim piece as the Shared record, both record directions
//! register on the face like an interior cut, and every node where the
//! rim piece starts or ends reads one arrival or departure too many
//! (measured: the six-axis arm's link face A 3v1, a plane bounded by a
//! 42.5-radius rim; records n17 -> n1 -> n14 and n9 -> n21 -> n12 lay on
//! it, four nodes imbalanced, 3D-Newell fallback on every pass).
//!
//! `chart_boundary_substitute` runs only with coincident face pairs,
//! against a 4-sample polyline of the original loop that cannot hold an
//! arc (measured: the 42.5 rim became a pentagon and every rim record
//! read off it, `on_poly=false`). Here the test is exact and in 3D: a pair is
//! a substitute when every interior sample of its run lies on one
//! original edge of the face, by that edge's own curve, and the kept
//! direction is the one agreeing with the face's half-edge on that
//! edge, which keeps the material on its left by the winding invariant.
//! Runs only on faces left imbalanced, so a balanced face never pays.

use popkern_core::geometry::curve::closest_point_in_window;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, SurfaceEval};
use popkern_core::topology::BRep;

use super::chart_subdivision::{FaceChartSubdivision, FilteredRecord};
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph};

/// Resolve edge-coincident substitute pairs on `sub` in place: of each
/// Shared pair whose run lies on an original edge of `face`, drop the
/// direction that runs against the face's half-edge there. `on_edge`
/// is the 3-D distance within which a sample lies on an edge curve.
pub(super) fn resolve_edge_coincident_substitutes(
    graph: &IntersectionGraph,
    face: &FaceRef,
    sub: &mut FaceChartSubdivision,
    surface: &SurfaceGeometry,
    original_brep: &BRep,
    on_edge: f64,
) {
    let Some(original_face) = original_brep.faces.get(face.face) else {
        return;
    };
    // The face's own half-edges: each edge's curve and the traversal
    // sign of this face's side along it.
    let mut sides = Vec::new();
    for l in
        std::iter::once(original_face.outer_loop).chain(original_face.inner_loops.iter().copied())
    {
        let Some(lp) = original_brep.loops.get(l) else {
            continue;
        };
        let start = lp.half_edge;
        let mut h = start;
        for _ in 0..=original_brep.half_edges.len() {
            let Some(he) = original_brep.half_edges.get(h) else {
                break;
            };
            if let Some(edge) = original_brep.edges.get(he.edge)
                && let Some(curve) = edge.curve.as_ref()
            {
                // The canonical side runs t_start -> t_end (PK-303).
                let forward = (edge.half_edge == h) == (curve.t_end >= curve.t_start);
                sides.push((curve, if forward { 1.0 } else { -1.0 }));
            }
            h = he.next;
            if h == start {
                break;
            }
        }
    }
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
    for (h, tw) in pairs {
        let points: Vec<_> = sub.embeddings[&h]
            .uv
            .iter()
            .map(|&(u, v)| surface.point_at(u, v))
            .collect();
        if points.len() < 3 {
            continue;
        }
        let interior = &points[1..points.len() - 1];
        // The original edge every interior sample lies on, and the
        // agreement of h's travel with the face's side along it.
        let verdict = sides.iter().find_map(|&(curve, sign)| {
            let (lo, hi) = (
                curve.t_start.min(curve.t_end),
                curve.t_start.max(curve.t_end),
            );
            let mut agree = 0.0;
            for (k, p) in interior.iter().enumerate() {
                let t = closest_point_in_window(&curve.geometry, p, lo, hi);
                if curve.geometry.point_at(t).distance_to(p) > on_edge {
                    return None;
                }
                let travel = points[k + 2] - points[k];
                agree += travel.dot(&curve.geometry.derivative_at(t, 1)) * sign;
            }
            Some(agree)
        });
        let Some(agree) = verdict.filter(|a| *a != 0.0) else {
            continue;
        };
        let (kept, dropped) = if agree > 0.0 { (h, tw) } else { (tw, h) };
        if sub.embeddings.remove(&dropped).is_some() {
            if popkern_core::trace::trace_enabled("boundary-substitute") {
                eprintln!(
                    "      [esub] {face:?} pair ({h},{tw}): on an original edge, kept {kept}"
                );
            }
            sub.report
                .filtered
                .push(FilteredRecord::BoundarySubstituteDirection { dropped, kept });
        }
    }
}

/// Whether some node of `sub` has unequal arrivals and departures.
pub(super) fn imbalanced(
    sub: &FaceChartSubdivision,
    ends: &dyn Fn(usize) -> (usize, usize),
) -> bool {
    let mut degree: std::collections::HashMap<usize, i64> = std::collections::HashMap::new();
    for &h in sub.embeddings.keys() {
        let (a, b) = ends(h);
        *degree.entry(a).or_default() += 1;
        *degree.entry(b).or_default() -= 1;
    }
    degree.values().any(|&d| d != 0)
}
