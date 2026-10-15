//! Inset scaffolds at oblique mitre apex corners (PK-841).
//!
//! On a support face whose corner angle differs from its neighbours'
//! (a d10 kite's 128° equator corner, a square pyramid's base), the
//! per-edge strips end at their own setback stations, so the served
//! sliver between anchors `a_i` and `a_j` carries extra rail vertices:
//! `V — a_i — u — C — w — a_j`, where `C` is the rail crossing (on both
//! bevel planes) and `u`, `w` are the points where the anchor edges
//! land on the rails (each on one bevel plane). When the rails stop
//! short instead (`V — a_i — u — w — a_j`, no vertex on both bevels),
//! `C` is the point where the support plane and the two bevel planes
//! meet. Either way every scaffold vertex lies on its rail's line, so
//! contracting it into `C` keeps each rail straight.
//!
//! **Repair, by Euler operators:** `kve` contracts each scaffold
//! vertex into the survivor (`C`, or the first scaffold vertex moved to
//! the three-plane point), and every line edge at the survivor is
//! re-authored between its end vertices. The cluster then has the
//! shape the mitre apex construction expects (one inset corner reached
//! from two anchors). Pcurves of the touched edges go stale and are
//! refilled by the chamfer driver's post-apply ladder.
//!
//! **Preconditions, verified before any surgery:** every scaffold
//! vertex lies on the sliver's support plane and on one or two bevel
//! planes; at most one lies on two. Anything else is returned as a
//! reason and the caller declines the corner.

use super::mitre_apex_corner::{Plane, destination, outgoing, plane_distance};
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

/// The point where three planes meet, `None` when they share a line or
/// two are parallel.
fn three_plane_point(a: &Plane, b: &Plane, c: &Plane) -> Option<Point3> {
    let (n1, n2, n3) = (a.1, b.1, c.1);
    let det = n1.dot(&n2.cross(&n3));
    if det.abs() < 1e-9 {
        return None;
    }
    let origin = Point3::new(0.0, 0.0, 0.0);
    let (d1, d2, d3) = (
        n1.dot(&(a.0 - origin)),
        n2.dot(&(b.0 - origin)),
        n3.dot(&(c.0 - origin)),
    );
    let v: Vec3 = n2.cross(&n3) * d1 + n3.cross(&n1) * d2 + n1.cross(&n2) * d3;
    Some(Point3::new(v.x / det, v.y / det, v.z / det))
}

fn position(brep: &BRep, v: VertexHandle) -> Result<Point3, String> {
    brep.vertices
        .get(v)
        .map(|d| d.position)
        .ok_or_else(|| format!("stale vertex {v:?}"))
}

/// Indices of the planes in `planes` holding `p`.
fn planes_holding(planes: &[Plane], p: &Point3, tol: f64) -> Vec<usize> {
    (0..planes.len())
        .filter(|&i| plane_distance(&planes[i], p) <= tol)
        .collect()
}

/// The loop ring starting at `he`'s origin.
fn ring_from(brep: &BRep, he: HalfEdgeHandle) -> Result<Vec<VertexHandle>, String> {
    let mut ring = Vec::new();
    let mut cursor = he;
    for _ in 0..=brep.half_edges.len() {
        let data = brep.half_edges.get(cursor).ok_or("stale ring half-edge")?;
        ring.push(data.origin);
        cursor = data.next;
        if cursor == he {
            return Ok(ring);
        }
    }
    Err("unterminated sliver loop".into())
}

/// Re-author every line edge at `v` between its end vertices.
fn reauthor_lines_at(out: &mut BRep, v: VertexHandle, tol: f64) -> Result<(), String> {
    for he in outgoing(out, v) {
        let edge = out.half_edges.get(he).ok_or("stale rail half-edge")?.edge;
        let edge_data = out.edges.get(edge).ok_or("stale rail edge")?;
        let is_line = edge_data
            .curve
            .as_ref()
            .is_some_and(|c| matches!(c.geometry, CurveGeometry::Line { .. }));
        if !is_line {
            continue;
        }
        // Author along the edge's own half-edge (its curve sense).
        let own = edge_data.half_edge;
        let from = out
            .half_edges
            .get(own)
            .ok_or("stale rail half-edge")?
            .origin;
        let to = destination(out, own).ok_or("rail without destination")?;
        let curve = CurveOnEdge::line_segment(position(out, from)?, position(out, to)?, tol)
            .ok_or_else(|| format!("degenerate rail at inset corner {v:?}"))?;
        out.edges.get_mut(edge).ok_or("stale rail edge")?.curve = Some(curve);
    }
    Ok(())
}

/// Contract every sliver's inset scaffold around `vertex` onto its
/// exact inset corner. Returns how many scaffold vertices were
/// contracted (0 on a corner whose rails already meet at the anchors'
/// shared inset).
pub(crate) fn contract_inset_scaffolds(
    out: &mut BRep,
    vertex: VertexHandle,
    bevels: &[Plane],
    supports: &[Plane],
    tol: f64,
) -> Result<usize, String> {
    let mut contracted = 0;
    for start in outgoing(out, vertex) {
        // V, a_i, scaffold…, a_j.
        let ring = ring_from(out, start)?;
        if ring.len() <= 4 {
            continue;
        }
        let scaffold = &ring[2..ring.len() - 1];
        let mut on_bevels = Vec::with_capacity(scaffold.len());
        for &s in scaffold {
            on_bevels.push(planes_holding(bevels, &position(out, s)?, tol));
        }
        let first = position(out, scaffold[0])?;
        let support = planes_holding(supports, &first, tol)
            .into_iter()
            .find(|&k| {
                scaffold.iter().all(|&s| {
                    position(out, s).is_ok_and(|p| plane_distance(&supports[k], &p) <= tol)
                })
            })
            .ok_or_else(|| format!("sliver scaffold at {vertex:?} is off every support plane"))?;
        let crossings: Vec<usize> = (0..scaffold.len())
            .filter(|&i| on_bevels[i].len() >= 2)
            .collect();
        if on_bevels.iter().any(|b| b.is_empty()) || crossings.len() > 1 {
            return Err(format!(
                "sliver scaffold at {vertex:?} is not a rail chain (bevel hits {on_bevels:?})"
            ));
        }
        let (survivor, moved_to) = match crossings.first() {
            Some(&i) => (scaffold[i], None),
            None => {
                let (bi, bj) = (on_bevels[0][0], on_bevels[scaffold.len() - 1][0]);
                if bi == bj {
                    return Err(format!(
                        "sliver scaffold at {vertex:?} lies on one bevel, not two"
                    ));
                }
                let c = three_plane_point(&supports[support], &bevels[bi], &bevels[bj])
                    .ok_or_else(|| format!("support and bevels share a line at {vertex:?}"))?;
                (scaffold[0], Some(c))
            }
        };
        // Contract the chain into the survivor one neighbour at a time.
        let mut rest: Vec<VertexHandle> = scaffold
            .iter()
            .copied()
            .filter(|&s| s != survivor)
            .collect();
        while !rest.is_empty() {
            let step = outgoing(out, survivor).into_iter().find_map(|he| {
                let d = destination(out, he)?;
                rest.iter().position(|&r| r == d).map(|k| (he, k))
            });
            let Some((he, k)) = step else {
                return Err(format!(
                    "scaffold vertex {:?} at {vertex:?} is not chained to its inset",
                    rest[0]
                ));
            };
            out.kve(he)
                .map_err(|e| format!("kve scaffold at {vertex:?}: {e}"))?;
            rest.swap_remove(k);
            contracted += 1;
        }
        if let Some(c) = moved_to {
            out.vertices
                .get_mut(survivor)
                .ok_or("stale inset survivor")?
                .position = c;
        }
        reauthor_lines_at(out, survivor, tol)?;
    }
    Ok(contracted)
}
