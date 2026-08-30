//! Trihedral mitre corner for planar chamfers (PK-378).
//!
//! When all three edges at a valence-3 vertex are chamfered with
//! planar bevels, the exact result is the half-space intersection: the
//! three bevel planes concur at the mitre apex `A` and each bevel ends
//! in a point there (a cube with every edge chamfered has 18 faces —
//! no corner triangle). The per-edge strip pipeline cannot express
//! that: measured on the unit cube, the served cluster keeps the
//! original vertex `V` with three stub edges to the anchors `a_i`
//! (the per-edge setback points on the chamfered edges themselves),
//! three sliver pieces of the support faces, and each bevel's loop
//! running `C — a_i — C'` through an anchor that lies `d/√2` off its
//! plane (`C` are the inset corners of the support faces, which ARE on
//! two bevels each).
//!
//! **Construction, by Euler operators on the served cluster:**
//! `kef` × 2 and `kev` collapse the three slivers and `V` into one
//! hexagon `a_x C_z a_y C_x a_z C_y`; `mev` inserts `A` from one inset
//! corner and `mef` × 2 split the hexagon into three quads, each
//! holding one anchor between the two inset corners its bevel owns;
//! per anchor, `kef` merges the quad into that bevel and `kev` deletes
//! the anchor. Net `−3V −6E −3F`: `V`, `a_x`, `a_y`, `a_z` → `A` and
//! nine cluster edges → three mitre edges `A–C_i`, authored as exact
//! line segments. Every bevel keeps its face (kef survivor = bevel
//! side), so names and surfaces are untouched.
//!
//! **Transactional:** the cluster is verified before surgery (valence
//! pattern, inset corners on two bevels each, planes concurrent) and
//! the body is validated after every corner; any mismatch declines
//! with a structured `InvalidTopology` naming the vertex, and the
//! caller keeps its input.

use super::EdgeGeometry;
use crate::error::BlendFailure;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};
use std::collections::BTreeMap;

/// Outcome of the stage.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct TrihedralMitreReport {
    /// Corners rebuilt with a mitre apex.
    pub corners_settled: usize,
}

/// A plane as `(point, unit normal)`.
type Plane = (Point3, Vec3);

fn plane_of(surface: &SurfaceGeometry) -> Option<Plane> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*origin, normal.normalized()?)),
        _ => None,
    }
}

fn plane_distance(plane: &Plane, p: &Point3) -> f64 {
    (*p - plane.0).dot(&plane.1).abs()
}

/// Intersection of three planes; `None` when they do not meet at one
/// point (parallel pair or a common line).
fn three_plane_apex(a: &Plane, b: &Plane, c: &Plane) -> Option<Point3> {
    let (n1, n2, n3) = (a.1, b.1, c.1);
    let det = n1.dot(&n2.cross(&n3));
    if det.abs() < 1e-9 {
        return None;
    }
    let (d1, d2, d3) = (
        n1.dot(&(a.0 - Point3::new(0.0, 0.0, 0.0))),
        n2.dot(&(b.0 - Point3::new(0.0, 0.0, 0.0))),
        n3.dot(&(c.0 - Point3::new(0.0, 0.0, 0.0))),
    );
    let v = n2.cross(&n3) * d1 + n3.cross(&n1) * d2 + n1.cross(&n2) * d3;
    Some(Point3::new(v.x / det, v.y / det, v.z / det))
}

/// Outgoing half-edges of `vertex`, in slotmap order (derived state is
/// never trusted for a structural check).
fn outgoing(brep: &BRep, vertex: VertexHandle) -> Vec<HalfEdgeHandle> {
    brep.half_edges
        .iter()
        .filter(|(_, he)| he.origin == vertex)
        .map(|(h, _)| h)
        .collect()
}

fn destination(brep: &BRep, he: HalfEdgeHandle) -> Option<VertexHandle> {
    let twin = brep.half_edges.get(he)?.twin;
    Some(brep.half_edges.get(twin)?.origin)
}

fn half_edge_in_loop_from(
    brep: &BRep,
    loop_: LoopHandle,
    origin: VertexHandle,
) -> Option<HalfEdgeHandle> {
    let start = brep.loops.get(loop_)?.half_edge;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let he = brep.half_edges.get(cursor)?;
        if he.origin == origin {
            return Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            return None;
        }
    }
    None
}

fn failure(vertex: VertexHandle, position: Option<Point3>, reason: String) -> BlendFailure {
    BlendFailure::InvalidTopology {
        reason: format!("trihedral mitre corner at {vertex:?}: {reason}"),
        location: position,
    }
}

/// Corners of `input` where every incident edge is chamfered with a
/// planar bevel between planar supports: `(vertex, [geometry index; 3])`.
fn candidate_corners(
    input: &BRep,
    edge_geometries: &[EdgeGeometry],
) -> Vec<(VertexHandle, [usize; 3])> {
    let by_edge: BTreeMap<EdgeHandle, usize> = edge_geometries
        .iter()
        .enumerate()
        .map(|(i, g)| (g.analysis.edge, i))
        .collect();
    let mut out = Vec::new();
    for (vh, _) in input.vertices.iter() {
        let mut edges: Vec<EdgeHandle> = outgoing(input, vh)
            .into_iter()
            .filter_map(|he| input.half_edges.get(he).map(|h| h.edge))
            .collect();
        edges.sort();
        edges.dedup();
        if edges.len() != 3 {
            continue;
        }
        let Some(idx): Option<Vec<usize>> = edges.iter().map(|e| by_edge.get(e).copied()).collect()
        else {
            continue;
        };
        let all_planar = idx.iter().all(|&i| {
            let g = &edge_geometries[i];
            plane_of(&g.surface).is_some()
                && plane_of(&g.analysis.surface_1).is_some()
                && plane_of(&g.analysis.surface_2).is_some()
        });
        if all_planar {
            out.push((vh, [idx[0], idx[1], idx[2]]));
        }
    }
    out
}

/// The served cluster at a candidate corner, verified.
struct Cluster {
    /// Stub half-edges `V → a_i`.
    stubs: [HalfEdgeHandle; 3],
    /// Anchors `a_i` (destinations of the stubs).
    anchors: [VertexHandle; 3],
    /// Inset corners: the three vertices reached from two anchors each.
    insets: [VertexHandle; 3],
}

fn read_cluster(
    out: &BRep,
    vertex: VertexHandle,
    bevels: &[Plane; 3],
    tol: f64,
) -> Result<Cluster, String> {
    let stubs = outgoing(out, vertex);
    if stubs.len() != 3 {
        return Err(format!(
            "expected 3 stub edges at the corner, found {}",
            stubs.len()
        ));
    }
    let mut anchors = Vec::new();
    let mut inset_hits: BTreeMap<VertexHandle, usize> = BTreeMap::new();
    for &stub in &stubs {
        let a = destination(out, stub).ok_or("stub without destination")?;
        let a_out = outgoing(out, a);
        if a_out.len() != 3 {
            return Err(format!(
                "anchor {a:?} has valence {}, expected 3",
                a_out.len()
            ));
        }
        for he in a_out {
            let d = destination(out, he).ok_or("anchor edge without destination")?;
            if d != vertex {
                *inset_hits.entry(d).or_insert(0) += 1;
            }
        }
        anchors.push(a);
    }
    if inset_hits.len() != 3 || inset_hits.values().any(|&n| n != 2) {
        return Err(format!(
            "expected 3 inset corners reached twice each, found {inset_hits:?}"
        ));
    }
    let insets: Vec<VertexHandle> = inset_hits.keys().copied().collect();
    for &c in &insets {
        let p = out.vertices.get(c).ok_or("stale inset corner")?.position;
        let on = bevels
            .iter()
            .filter(|pl| plane_distance(pl, &p) <= tol)
            .count();
        if on < 2 {
            return Err(format!(
                "inset corner {c:?} lies on {on} bevel planes, expected 2"
            ));
        }
    }
    Ok(Cluster {
        stubs: [stubs[0], stubs[1], stubs[2]],
        anchors: [anchors[0], anchors[1], anchors[2]],
        insets: [insets[0], insets[1], insets[2]],
    })
}

/// Rebuild one corner. Works on `out` directly; the caller stages.
fn rebuild_corner(
    out: &mut BRep,
    vertex: VertexHandle,
    cluster: &Cluster,
    apex: Point3,
    tol: f64,
) -> Result<(), String> {
    // 1. Slivers + V → hexagon. The kef survivor is the half-edge's
    //    face; both slivers die anyway, so either side works.
    out.kef(cluster.stubs[0])
        .map_err(|e| format!("kef stub 0: {e}"))?;
    out.kef(cluster.stubs[1])
        .map_err(|e| format!("kef stub 1: {e}"))?;
    let into_v = out
        .half_edges
        .get(cluster.stubs[2])
        .ok_or("stale stub 2")?
        .twin;
    let kev = out
        .kev(into_v)
        .map_err(|e| format!("kev corner vertex: {e}"))?;
    let hexagon = kev.loop_;
    if hexagon_len(out, hexagon) != Some(6) {
        return Err(format!(
            "collapsed sliver loop has {:?} half-edges, expected 6",
            hexagon_len(out, hexagon)
        ));
    }
    // 2. Apex spur from C_0, then two mef splits to C_1 and C_2.
    let [c0, c1, c2] = cluster.insets;
    let he_c0 = half_edge_in_loop_from(out, hexagon, c0).ok_or("C0 not on hexagon")?;
    let mev = out.mev(he_c0, apex).map_err(|e| format!("mev apex: {e}"))?;
    let apex_v = mev.vertex;
    let mut mitre_edges = vec![mev.edge];
    for c in [c1, c2] {
        // Find a loop holding both an A-origin half-edge and `c`.
        let mut done = false;
        for he_a in outgoing(out, apex_v) {
            let loop_ = out
                .half_edges
                .get(he_a)
                .ok_or("stale apex half-edge")?
                .loop_;
            if let Some(he_c) = half_edge_in_loop_from(out, loop_, c) {
                let mef = out
                    .mef(he_a, he_c)
                    .map_err(|e| format!("mef apex-{c:?}: {e}"))?;
                mitre_edges.push(mef.edge);
                done = true;
                break;
            }
        }
        if !done {
            return Err(format!("no loop joins the apex and inset corner {c:?}"));
        }
    }
    // 3. Merge each quad into its bevel, delete the anchor.
    for &a in &cluster.anchors {
        let a_out = outgoing(out, a);
        if a_out.len() != 2 {
            return Err(format!(
                "anchor {a:?} valence {} before merge, expected 2",
                a_out.len()
            ));
        }
        // The quad is the face touching the apex; kef from the OTHER
        // side so the bevel survives.
        let he_a = a_out[0];
        let he_data = out.half_edges.get(he_a).ok_or("stale anchor half-edge")?;
        let twin = he_data.twin;
        let quad_side = if half_edge_in_loop_from(out, he_data.loop_, apex_v).is_some() {
            he_a
        } else {
            twin
        };
        let bevel_side = out.half_edges.get(quad_side).ok_or("stale quad side")?.twin;
        out.kef(bevel_side)
            .map_err(|e| format!("kef quad into bevel at {a:?}: {e}"))?;
        // `a` is now a valence-1 spur; kill it via the half-edge that
        // arrives at it.
        let remaining = outgoing(out, a);
        if remaining.len() != 1 {
            return Err(format!(
                "anchor {a:?} valence {} after merge, expected 1",
                remaining.len()
            ));
        }
        let into_a = out.half_edges.get(remaining[0]).ok_or("stale spur")?.twin;
        out.kev(into_a)
            .map_err(|e| format!("kev anchor {a:?}: {e}"))?;
    }
    // 4. Exact line curves on the mitre edges.
    for edge in mitre_edges {
        let he = out.edges.get(edge).ok_or("stale mitre edge")?.half_edge;
        let v0 = out
            .half_edges
            .get(he)
            .ok_or("stale mitre half-edge")?
            .origin;
        let v1 = destination(out, he).ok_or("mitre edge without destination")?;
        let p0 = out.vertices.get(v0).ok_or("stale mitre vertex")?.position;
        let p1 = out.vertices.get(v1).ok_or("stale mitre vertex")?.position;
        let curve = CurveOnEdge::line_segment(p0, p1, tol)
            .ok_or_else(|| format!("degenerate mitre edge {edge:?}"))?;
        out.edges.get_mut(edge).ok_or("stale mitre edge")?.curve = Some(curve);
    }
    let _ = vertex;
    Ok(())
}

fn hexagon_len(brep: &BRep, loop_: LoopHandle) -> Option<usize> {
    let start = brep.loops.get(loop_)?.half_edge;
    let mut n = 0;
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        n += 1;
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(n);
        }
    }
    None
}

/// Settle every trihedral all-chamfered corner of `input` in `out`
/// (the post-apply body). Transactional per corner: on any verification
/// or Euler failure the body is left as it was and the failure names
/// the corner; the caller declines the whole operation.
pub(crate) fn settle_trihedral_mitre_corners(
    input: &BRep,
    out: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<TrihedralMitreReport, (EdgeHandle, BlendFailure)> {
    let tol = tolerance.linear().max(1e-9) * 10.0;
    let mut report = TrihedralMitreReport::default();
    for (vertex, idx) in candidate_corners(input, edge_geometries) {
        let position = input.vertices.get(vertex).map(|v| v.position);
        let owner_edge = edge_geometries[idx[0]].analysis.edge;
        let decline = |reason: String| (owner_edge, failure(vertex, position, reason));
        let bevels: Vec<Plane> = idx
            .iter()
            .filter_map(|&i| plane_of(&edge_geometries[i].surface))
            .collect();
        let bevels: [Plane; 3] = match bevels.as_slice() {
            [a, b, c] => [*a, *b, *c],
            _ => return Err(decline("bevel surfaces are not all planes".into())),
        };
        let apex = three_plane_apex(&bevels[0], &bevels[1], &bevels[2])
            .ok_or_else(|| decline("bevel planes do not meet at a single apex".into()))?;
        let cluster = read_cluster(out, vertex, &bevels, tol).map_err(decline)?;
        let mut staged = out.clone();
        rebuild_corner(&mut staged, vertex, &cluster, apex, tol).map_err(decline)?;
        if let Err(e) = staged.validate() {
            return Err(decline(format!("rebuilt corner fails validation: {e}")));
        }
        *out = staged;
        report.corners_settled += 1;
    }
    Ok(report)
}
