//! Removal of operand edges duplicated by a shared intersection edge
//! (split from `graph_split.rs` under the file-size ratchet).
//!
//! After `split_edges_at_shared_nodes`, an operand edge segment may
//! connect the same pair of nodes as a shared intersection curve on
//! the same face — a coplanar overlap's record running along the
//! operand's own boundary. The walk cannot distinguish the two, so the
//! operand duplicate is removed.
//!
//! A duplicate is the SAME CURVE, not merely the same endpoints. Two
//! distinct curves can join one node pair: on a flush cap of two
//! overlapping cylinders, the cap's outer rim arc and the other
//! cylinder's arc inside the cap both run between the two wall
//! junctions. Removing the rim arc as a "duplicate" left the cap with
//! only its lens to walk — the crescent never became a fragment and
//! the two-cylinder union served a zero-volume body (measured).
//! Endpoint coincidence stays exact; curve coincidence is tested at
//! the operand segment's midpoint against the shared curve, whose
//! parameterization may differ (a fitted record against an exact
//! rational arc), so the test is a point-to-curve distance.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};

use super::graph::{EdgeSource, IntersectionGraph};

/// Endpoint coincidence: both curves end at the same graph nodes.
const ENDPOINT_TOL: f64 = 1e-10;

/// Curve coincidence: the operand segment's midpoint lies within this
/// fraction of the segment's extent from the shared curve. Generous
/// against record fit slop; distinct curves through one node pair
/// separate by a sizeable fraction of the chord.
const MIDPOINT_FRACTION: f64 = 0.05;

/// Remove non-shared half-edges that duplicate a shared half-edge on
/// the same face: same end nodes AND the same curve (see the module
/// docs).
pub(crate) fn remove_collinear_duplicates(graph: &mut IntersectionGraph, linear: f64) {
    let mut to_remove: Vec<usize> = Vec::new();
    let mut seam_substitutes: Vec<(super::graph::FaceRef, usize)> = Vec::new();

    for (face, he_indices) in &graph.face_edges {
        let ends = |idx: usize| -> (Point3, Point3) {
            let he = &graph.half_edges[idx];
            (
                graph.nodes[he.origin].position,
                graph.nodes[graph.half_edges[he.twin].origin].position,
            )
        };
        let shared: Vec<usize> = he_indices
            .iter()
            .copied()
            .filter(|&idx| matches!(graph.half_edges[idx].source, EdgeSource::Shared { .. }))
            .collect();
        if shared.is_empty() {
            continue;
        }

        for &idx in he_indices {
            if matches!(graph.half_edges[idx].source, EdgeSource::Shared { .. }) {
                continue;
            }
            let (p0, p1) = ends(idx);
            // Adaptive squared-distance predicate: the decision is
            // order-independent at the tolerance boundary.
            let close = |a: &Point3, b: &Point3| -> bool {
                popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                    (a.x, a.y, a.z),
                    (b.x, b.y, b.z),
                    ENDPOINT_TOL * ENDPOINT_TOL,
                )
            };
            let dup_of = shared.iter().copied().find(|&s| {
                let (sp0, sp1) = ends(s);
                let same_ends = (close(&sp0, &p0) && close(&sp1, &p1))
                    || (close(&sp0, &p1) && close(&sp1, &p0));
                same_ends && same_curve(curve_of(graph, idx), curve_of(graph, s), (p0, p1), linear)
            });
            if let Some(s) = dup_of {
                to_remove.push(idx);
                // The operand edge was this face's SEAM (both halves
                // on the face): the record now carries the seam, and
                // the chart must give it two border copies, not treat
                // it as an interior cut (PK-750).
                let twin = graph.half_edges[idx].twin;
                if graph.half_edges[twin].face.as_ref() == Some(face) {
                    seam_substitutes.push((*face, s));
                }
            }
        }
    }

    for (face, s) in seam_substitutes {
        let twin = graph.half_edges[s].twin;
        graph.seam_substitutes.insert((face, s));
        graph.seam_substitutes.insert((face, twin));
    }
    for &idx in &to_remove {
        if let Some(f) = graph.half_edges[idx].face
            && let Some(edges) = graph.face_edges.get_mut(&f)
        {
            edges.retain(|&e| e != idx);
        }
        let origin = graph.half_edges[idx].origin;
        graph.nodes[origin].outgoing.retain(|&e| e != idx);
        graph.half_edges[idx].face = None;
    }
}

/// The curve of half-edge `idx`, read from its twin when the pair
/// stores it on one side only (shared record pairs do).
fn curve_of(graph: &IntersectionGraph, idx: usize) -> Option<&CurveOnEdge> {
    let he = &graph.half_edges[idx];
    he.curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
}

/// Whether the operand segment `operand` (ending at `p0`, `p1`) lies
/// on the shared curve `shared`. Without a curve on either side the
/// endpoint match alone decides (the historical straight-segment
/// contract).
///
/// The band is the smaller of the chord fraction and 100 × linear: fit
/// slop scales with the model tolerance, while a THIN lens is still
/// two curves (measured, PK-758: a boss rising 0.5 above a ring's notch
/// edge over a 13-long chord sat inside 5% of the chord, the notch edge
/// was removed as a "duplicate" of the boss record, and the ring wall's
/// walk lost a boundary piece).
fn same_curve(
    operand: Option<&CurveOnEdge>,
    shared: Option<&CurveOnEdge>,
    (p0, p1): (Point3, Point3),
    linear: f64,
) -> bool {
    let (Some(operand), Some(shared)) = (operand, shared) else {
        return true;
    };
    let mid = operand
        .geometry
        .point_at(0.5 * (operand.t_start + operand.t_end));
    let extent = p0.distance_to(&p1).max(mid.distance_to(&p0));
    let band = (MIDPOINT_FRACTION * extent).min(100.0 * linear);
    distance_to_piece(shared, &mid) <= band
}

/// The distance from `p` to the PIECE of `curve` over its own
/// parameter window, not to its whole carrier. A record on a rim circle
/// carries the full circle with a window `[0.75, 1.25]`. Projecting onto
/// the carrier put the rim's OTHER arc (same end nodes, the far side of
/// the circle) at distance 0, so a box whose caps are flush with a
/// rod's dropped the rod's far rim arc as a "duplicate" and lost half
/// its wall (PK-750).
fn distance_to_piece(curve: &CurveOnEdge, p: &Point3) -> f64 {
    closest_on_piece(curve, p).1
}

/// The closest point of `curve`'s PIECE (its own parameter window) to
/// `p`: `(t, distance)`. Dense samples over the window seed a Newton
/// refinement clamped to it.
///
/// A CLOSED piece (a full rim circle, window [0, 2π]) has two window
/// ends at one point: a foot just short of the closing end tied with
/// the opening end's sample, Newton from t = 0 stepped negative and
/// was clamped back, and the piece read 0.30 away from a point lying
/// on it (PK-777: the periodic boundary-substitute resolver then kept
/// both directions of a rim record and the wall's chart declined).
/// Both ends seed the polish there.
pub(crate) fn closest_on_piece(curve: &CurveOnEdge, p: &Point3) -> (f64, f64) {
    const SAMPLES: usize = 32;
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    let at = |t: f64| curve.geometry.point_at(t).distance_to(p);
    let mut t = lo;
    let mut best = at(lo);
    for k in 1..=SAMPLES {
        let s = lo + (hi - lo) * k as f64 / SAMPLES as f64;
        let d = at(s);
        if d < best {
            (t, best) = (s, d);
        }
    }
    let polish = |mut t: f64, mut best: f64| -> (f64, f64) {
        for _ in 0..8 {
            let diff = curve.geometry.point_at(t) - *p;
            let d1 = curve.geometry.derivative_at(t, 1);
            let fp = d1.dot(&d1) + diff.dot(&curve.geometry.derivative_at(t, 2));
            if fp.abs() < 1e-30 {
                break;
            }
            let next = (t - diff.dot(&d1) / fp).clamp(lo, hi);
            let d = at(next);
            if d >= best {
                break;
            }
            (t, best) = (next, d);
        }
        (t, best)
    };
    let found = polish(t, best);
    let closed = hi > lo
        && curve
            .geometry
            .point_at(lo)
            .distance_to(&curve.geometry.point_at(hi))
            <= 1e-12 * (1.0 + p.distance_to(&curve.geometry.point_at(lo)));
    if closed && (t == lo || t == hi) {
        let other = if t == lo { hi } else { lo };
        let alt = polish(other, at(other));
        if alt.1 < found.1 {
            return alt;
        }
    }
    found
}

#[cfg(test)]
#[path = "graph_duplicates_tests.rs"]
mod tests;
