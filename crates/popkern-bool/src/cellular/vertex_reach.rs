//! A decided vertex records how far its edges' ends lie from it
//! (PK-1181).
//!
//! The vertex table joins candidates within their tolerances into one
//! vertex at the group's anchor, growing the group's tolerance to cover
//! every member; a record split there by a member it passes within the
//! reach of ends at the anchor though its curve passes beside it. The
//! assembly builds the vertex at the anchor and the edges on their
//! curves, so an edge's end can lie off its vertex by up to the group's
//! tolerance: measured, two balls whose rims cut a box's top exactly a
//! tolerance apart, the second ball's rim split at the first ball's
//! seam vertex 1.0000006e-6 off it, and the result's loop read open
//! against the bare tolerance. Each vertex at an anchor takes, as its
//! tolerance, the farthest of its edges' ends (no farther than the
//! group's tolerance: the resolution's decision, not a repair), as an
//! imported curve's end off its vertex records the file's precision
//! (PK-639) and a glued contact's end its distance (PK-1170); every
//! edge ending there takes at least that tolerance (an edge's tolerance
//! contains its vertices').

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// Give each vertex of `body` at an anchor of `reaches` the tolerance
/// its edges' ends need, within the anchor's grown tolerance.
pub(super) fn record_vertex_reaches(body: &mut BRep, reaches: &[(Point3, f64)]) {
    let linear = body.tolerance.linear();
    if reaches.iter().all(|&(_, t)| t <= linear) {
        return;
    }
    // The farthest end of an edge from each of its vertices.
    let mut needed: std::collections::HashMap<VertexHandle, f64> = Default::default();
    let mut ends_at: std::collections::HashMap<VertexHandle, Vec<EdgeHandle>> = Default::default();
    for (e, edge) in body.edges.iter() {
        let (Some(c), Some((va, vb))) = (edge.curve.as_ref(), body.edge_vertices(e)) else {
            continue;
        };
        let ends = [c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end)];
        for v in [va, vb] {
            let Some(p) = body.vertices.get(v).map(|x| x.position) else {
                continue;
            };
            let d = ends
                .iter()
                .map(|q| q.distance_to(&p))
                .fold(f64::INFINITY, f64::min);
            let slot = needed.entry(v).or_insert(0.0);
            *slot = slot.max(d);
            ends_at.entry(v).or_default().push(e);
        }
    }
    let at_anchor = |p: &Point3| {
        reaches
            .iter()
            .filter(|(a, _)| a.distance_to(p) <= 1e-3 * linear)
            .map(|&(_, t)| t)
            .fold(None, |m: Option<f64>, t| Some(m.map_or(t, |m| m.max(t))))
    };
    for (v, d) in needed {
        if d <= linear {
            continue;
        }
        let Some(cap) = body.vertices.get(v).and_then(|x| at_anchor(&x.position)) else {
            continue;
        };
        if d > cap {
            continue;
        }
        let at = d * (1.0 + 1e-9);
        if let Some(vx) = body.vertices.get_mut(v) {
            vx.tolerance = Some(vx.tolerance.map_or(at, |t| t.max(at)));
        }
        for &e in ends_at.get(&v).into_iter().flatten() {
            let raised = body.edge_tolerance(e).max(at);
            if let Some(edge) = body.edges.get_mut(e) {
                edge.tolerance = Some(raised);
            }
        }
    }
}

#[cfg(test)]
#[path = "vertex_reach_tests.rs"]
mod tests;
