//! Redundant collinear valence-2 vertex dissolution (boolean output
//! hygiene).
//!
//! A coincident weld can leave a valence-2 vertex joining two
//! COLLINEAR straight edges — the point where an operand boundary
//! ended on a flush interface (an L-union's wall corner keeps a
//! vertex at the old plate height, splitting one physical 35 mm
//! corner line into 6 + 29). The body is watertight, but segment
//! semantics leak downstream: chord-length edge selectors see HALF
//! of a physical corner and dress-up features act on partial edges
//! (measured: a chord∈[5,7] fillet selector caught the 6 mm piece,
//! filleted half the corner, and the collinear continuation ended
//! 3(√2−1) off the trimmed support — check 9 rejected the body).
//!
//! Scope: STRAIGHT edges only (Line curves or curveless chords).
//! Curved-edge splits are deliberately out of scope — periodic rim
//! circles keep their seam vertices as period anchors, and merging
//! arcs would re-author curve geometry this pass has no business
//! deciding. Each dissolution is one transactional `kve` (the core
//! operator validates structure and declines rather than corrupt),
//! followed by re-authoring the surviving edge's line to span the
//! merged extent.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// True when the edge is geometrically straight within `tol`: a Line
/// curve, no curve at all (curveless edges are chords by
/// construction), or a sampled curve whose interior stays on the
/// endpoint chord — intersection records author STRAIGHT flush-weld
/// traces as NURBS, and a variant test alone would skip exactly the
/// class this pass exists for.
fn is_straight(brep: &BRep, eh: EdgeHandle, tol: f64) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    match brep.edges.get(eh).and_then(|e| e.curve.as_ref()) {
        None => true,
        Some(c) => {
            if matches!(c.geometry, CurveGeometry::Line { .. }) {
                return true;
            }
            let p0 = c.geometry.point_at(c.t_start);
            let p1 = c.geometry.point_at(c.t_end);
            let chord = p1 - p0;
            let l2 = chord.dot(&chord);
            if l2 < tol * tol {
                return false; // zero-extent: not a mergeable straight
            }
            (1..8).all(|i| {
                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                let p = c.geometry.point_at(t);
                let s = ((p - p0).dot(&chord) / l2).clamp(0.0, 1.0);
                (p0 + chord * s).distance_to(&p) <= tol
            })
        }
    }
}

/// The unordered face pair an edge bounds, or `None` when either
/// side is unreadable.
fn face_pair(
    brep: &BRep,
    eh: EdgeHandle,
) -> Option<(
    popkern_core::topology::FaceHandle,
    popkern_core::topology::FaceHandle,
)> {
    let he = brep.edges.get(eh)?.half_edge;
    let h = brep.half_edges.get(he)?;
    let t = brep.half_edges.get(h.twin)?;
    let f1 = brep.loops.get(h.loop_)?.face;
    let f2 = brep.loops.get(t.loop_)?.face;
    Some((f1.min(f2), f1.max(f2)))
}

/// Dissolve every valence-2 vertex whose two straight edges are
/// collinear, merging each pair into one edge spanning the full
/// extent. Returns the number of vertices dissolved. Pcurves on
/// surviving edges are re-derived by the caller's recompute pass.
pub(crate) fn dissolve_collinear_valence2_vertices(
    brep: &mut BRep,
    tolerance: &Tolerance,
) -> usize {
    let tol = tolerance.linear();
    let mut dissolved = 0usize;
    // Fixpoint: dissolving one vertex can expose the next along the
    // same line (a boundary split into three pieces). Bounded by the
    // vertex count — each iteration removes one vertex or stops.
    for _ in 0..brep.vertices.len().max(1) {
        let Some((m, keep_he, survivor_edge, a_pos, b_pos)) = find_candidate(brep, tol) else {
            break;
        };
        if brep.kve(keep_he).is_err() {
            // Transactional decline (would break structure): record
            // nothing here — the vertex simply stays, exactly the
            // pre-pass behavior.
            break;
        }
        let dir = b_pos - a_pos;
        let len = dir.length();
        if len > tol
            && let Some(edge) = brep.edges.get_mut(survivor_edge)
        {
            edge.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: a_pos,
                    direction: dir * (1.0 / len),
                },
                t_start: 0.0,
                t_end: len,
            });
        }
        // The surviving edge's stored pcurves still parameterize the
        // old sub-span; clear them so the recompute pass re-derives
        // both from the fresh line.
        let hes: Vec<_> = brep
            .half_edges
            .iter()
            .filter(|(_, he)| he.edge == survivor_edge)
            .map(|(h, _)| h)
            .collect();
        for h in hes {
            brep.half_edges[h].pcurve = None;
        }
        let _ = m;
        dissolved += 1;
    }
    dissolved
}

/// Find one dissolvable vertex: valence 2, two distinct straight
/// edges bounding the same face pair, collinear through the vertex.
///
/// Returns `(vertex, half_edge_to_contract, surviving_edge, far_a,
/// far_b)` where the half-edge runs FROM the far endpoint of the
/// edge being killed INTO the vertex (kve keeps the origin).
#[allow(clippy::type_complexity)]
fn find_candidate(
    brep: &BRep,
    tol: f64,
) -> Option<(
    VertexHandle,
    popkern_core::topology::HalfEdgeHandle,
    EdgeHandle,
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    // Handle order, not hash order: this scan returns the FIRST passing
    // vertex, and the fixpoint caller breaks on the first kve
    // decline — hash order would make the dissolved set and the
    // surviving-edge authoring order run-to-run nondeterministic
    // (PK-300 hash-iteration audit). A dense map iterates in the
    // handles' order (the BTreeMap's, slot index then version) without
    // a tree or a vector per vertex (PK-815: every boolean's whole
    // result); only each vertex's first two outgoing half-edges and
    // its count are kept, as only valence 2 is read.
    let mut outgoing: slotmap::SecondaryMap<
        VertexHandle,
        (usize, [popkern_core::topology::HalfEdgeHandle; 2]),
    > = slotmap::SecondaryMap::with_capacity(brep.vertices.len());
    for (h, he) in brep.half_edges.iter() {
        match outgoing.get_mut(he.origin) {
            Some((n, first)) => {
                if *n < 2 {
                    first[*n] = h;
                }
                *n += 1;
            }
            None => {
                outgoing.insert(he.origin, (1, [h, h]));
            }
        }
    }
    for (m, &(n, hes)) in &outgoing {
        let m = &m;
        if n != 2 {
            continue;
        }
        let (h1, h2) = (hes[0], hes[1]);
        let (e1, e2) = (brep.half_edges[h1].edge, brep.half_edges[h2].edge);
        if e1 == e2 {
            continue; // closed rim's seam vertex — load-bearing
        }
        if !is_straight(brep, e1, tol) || !is_straight(brep, e2, tol) {
            continue;
        }
        match (face_pair(brep, e1), face_pair(brep, e2)) {
            (Some(p1), Some(p2)) if p1 == p2 => {}
            _ => continue,
        }
        // Far endpoints: h1/h2 depart M, so their twins' origins are
        // the far ends.
        let a = brep.half_edges[brep.half_edges[h1].twin].origin;
        let b = brep.half_edges[brep.half_edges[h2].twin].origin;
        if a == b || a == *m || b == *m {
            continue;
        }
        let (pa, pb, pm) = (
            brep.vertices[a].position,
            brep.vertices[b].position,
            brep.vertices[*m].position,
        );
        // Collinear-through test: M lies on segment AB within the
        // linear tolerance (positions are the authority, matching
        // the weld that authored the split).
        let ab = pb - pa;
        let l2 = ab.dot(&ab);
        if l2 < tol * tol {
            continue;
        }
        let t = (pm - pa).dot(&ab) / l2;
        if !(0.0..=1.0).contains(&t) {
            continue;
        }
        let foot = pa + ab * t;
        if foot.distance_to(&pm) > tol {
            continue;
        }
        // Contract e1 via the half-edge running A → M (kve keeps the
        // origin A and merges M into it); e2 survives spanning A → B.
        let a_to_m = brep.half_edges[h1].twin;
        return Some((*m, a_to_m, e2, pa, pb));
    }
    None
}
