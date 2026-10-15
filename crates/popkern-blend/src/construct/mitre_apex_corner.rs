//! Mitre apex corners for planar chamfers (PK-378, PK-841).
//!
//! When every edge at a valence-`N` vertex (`N ≥ 3`) is chamfered with
//! a planar bevel and the `N` bevel planes concur, the exact result is
//! the half-space intersection: each bevel ends in a point at the
//! mitre apex `A` (a cube with every edge chamfered has 18 faces, an
//! octahedron 26, an icosahedron 62 — no corner polygon). The per-edge
//! strip pipeline cannot express that: measured on the unit cube, the
//! served cluster keeps the original vertex `V` with `N` stub edges to
//! the anchors `a_i` (the per-edge setback points on the chamfered
//! edges themselves), `N` sliver pieces of the support faces, and each
//! bevel's loop running `C — a_i — C'` through an anchor that lies off
//! its plane (`d/√2` on the cube, `d/√3` on the octahedron; `C` are the
//! inset corners of the support faces, which ARE on two bevels each).
//!
//! **Construction, by Euler operators on the served cluster:**
//! `kef` × (N−1) and `kev` collapse the `N` slivers and `V` into one
//! 2N-gon alternating anchors and inset corners; `mev` inserts `A`
//! from one inset corner and `mef` × (N−1) split it into `N` quads,
//! each holding one anchor between the two inset corners its bevel
//! owns; per anchor, `kef` merges the quad into that bevel and `kev`
//! deletes the anchor. `V` and the `N` anchors → `A`, and `3N` cluster
//! edges → `N` mitre edges `A–C_i`, authored as exact line segments.
//! Every bevel keeps its face (kef survivor = bevel side), so names and
//! surfaces are untouched.
//!
//! **Concurrency:** three independent planes always meet at a point;
//! for `N > 3` the apex is the least-squares point of the bevel planes
//! and every plane must pass within the linear band of it (regular
//! solids with equal distances concur by symmetry). Non-concurrent
//! bevels need a corner polygon face, which this stage does not build
//! — it declines with the worst plane residual (PK-848).
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
pub struct MitreApexReport {
    /// Corners rebuilt with a mitre apex.
    pub corners_settled: usize,
}

/// A plane as `(point, unit normal)`.
pub(super) type Plane = (Point3, Vec3);

fn plane_of(surface: &SurfaceGeometry) -> Option<Plane> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*origin, normal.normalized()?)),
        _ => None,
    }
}

pub(super) fn plane_distance(plane: &Plane, p: &Point3) -> f64 {
    (*p - plane.0).dot(&plane.1).abs()
}

/// The point where all `planes` meet, as the least-squares solution of
/// `n_i · x = n_i · o_i` (Cramer on the 3×3 normal equations; exact for
/// three planes). `Err` names why there is no single apex: a common
/// line or parallel pair (singular system), or the worst plane's
/// distance from the best point when it exceeds `tol`.
fn concurrent_apex(planes: &[Plane], tol: f64) -> Result<Point3, String> {
    let mut m = [[0.0f64; 3]; 3];
    let mut b = [0.0f64; 3];
    for (o, n) in planes {
        let n = [n.x, n.y, n.z];
        let d = n[0] * o.x + n[1] * o.y + n[2] * o.z;
        for r in 0..3 {
            for c in 0..3 {
                m[r][c] += n[r] * n[c];
            }
            b[r] += n[r] * d;
        }
    }
    let det3 = |a: &[[f64; 3]; 3]| {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let det = det3(&m);
    // Unit normals bound every eigenvalue by N, so the scale-free
    // singularity test is relative to N³.
    let n = planes.len() as f64;
    if det.abs() < 1e-9 * n * n * n {
        return Err("bevel planes do not meet at a single apex (common line or parallel)".into());
    }
    let mut x = [0.0f64; 3];
    for (col, xc) in x.iter_mut().enumerate() {
        let mut mc = m;
        for (row, bv) in b.iter().enumerate() {
            mc[row][col] = *bv;
        }
        *xc = det3(&mc) / det;
    }
    let apex = Point3::new(x[0], x[1], x[2]);
    let worst = planes
        .iter()
        .map(|pl| plane_distance(pl, &apex))
        .fold(0.0f64, f64::max);
    if worst > tol {
        return Err(format!(
            "{} bevel planes are not concurrent: worst plane lies {worst:.3e} off the \
             best apex (band {tol:.1e}); a corner polygon face is needed",
            planes.len()
        ));
    }
    Ok(apex)
}

/// Outgoing half-edges of `vertex`, in slotmap order (derived state is
/// never trusted for a structural check).
pub(super) fn outgoing(brep: &BRep, vertex: VertexHandle) -> Vec<HalfEdgeHandle> {
    brep.half_edges
        .iter()
        .filter(|(_, he)| he.origin == vertex)
        .map(|(h, _)| h)
        .collect()
}

pub(super) fn destination(brep: &BRep, he: HalfEdgeHandle) -> Option<VertexHandle> {
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
        reason: format!("mitre apex corner at {vertex:?}: {reason}"),
        location: position,
    }
}

/// Corners of `input` (valence ≥ 3) where every incident edge is
/// chamfered with a planar bevel between planar supports:
/// `(vertex, geometry indices)`.
fn candidate_corners(
    input: &BRep,
    edge_geometries: &[EdgeGeometry],
) -> Vec<(VertexHandle, Vec<usize>)> {
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
        if edges.len() < 3 {
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
            out.push((vh, idx));
        }
    }
    out
}

/// The served cluster at a candidate corner, verified.
struct Cluster {
    /// Stub half-edges `V → a_i`, one per chamfered edge.
    stubs: Vec<HalfEdgeHandle>,
    /// Anchors `a_i` (destinations of the stubs).
    anchors: Vec<VertexHandle>,
    /// Inset corners: the `N` vertices reached from two anchors each.
    insets: Vec<VertexHandle>,
}

fn read_cluster(
    out: &BRep,
    vertex: VertexHandle,
    bevels: &[Plane],
    tol: f64,
) -> Result<Cluster, String> {
    let n = bevels.len();
    let stubs = outgoing(out, vertex);
    if stubs.len() != n {
        return Err(format!(
            "expected {n} stub edges at the corner, found {}",
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
    if inset_hits.len() != n || inset_hits.values().any(|&hits| hits != 2) {
        return Err(format!(
            "expected {n} inset corners reached twice each, found {inset_hits:?}"
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
        stubs,
        anchors,
        insets,
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
    // 1. Slivers + V → 2N-gon. The N slivers form a cycle around V
    //    whose dual edges are the stubs; any N−1 of them form a path,
    //    so each kef joins two distinct faces. The kef survivor is the
    //    half-edge's face; every sliver dies anyway.
    let n = cluster.stubs.len();
    let (last, merged) = cluster.stubs.split_last().ok_or("no stubs")?;
    for (i, &stub) in merged.iter().enumerate() {
        out.kef(stub).map_err(|e| format!("kef stub {i}: {e}"))?;
    }
    let into_v = out.half_edges.get(*last).ok_or("stale last stub")?.twin;
    let kev = out
        .kev(into_v)
        .map_err(|e| format!("kev corner vertex: {e}"))?;
    let polygon = kev.loop_;
    if loop_len(out, polygon) != Some(2 * n) {
        return Err(format!(
            "collapsed sliver loop has {:?} half-edges, expected {}",
            loop_len(out, polygon),
            2 * n
        ));
    }
    // 2. Apex spur from C_0, then N−1 mef splits to the other insets.
    let (&c0, rest) = cluster.insets.split_first().ok_or("no inset corners")?;
    let he_c0 = half_edge_in_loop_from(out, polygon, c0).ok_or("C0 not on the polygon")?;
    let mev = out.mev(he_c0, apex).map_err(|e| format!("mev apex: {e}"))?;
    let apex_v = mev.vertex;
    let mut mitre_edges = vec![mev.edge];
    for &c in rest {
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

fn loop_len(brep: &BRep, loop_: LoopHandle) -> Option<usize> {
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

/// Settle every all-chamfered planar corner of `input` in `out`
/// (the post-apply body). Transactional per corner: on any verification
/// or Euler failure the body is left as it was and the failure names
/// the corner; the caller declines the whole operation.
pub(crate) fn settle_mitre_apex_corners(
    input: &BRep,
    out: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<MitreApexReport, (EdgeHandle, BlendFailure)> {
    let tol = tolerance.linear().max(1e-9) * 10.0;
    let mut report = MitreApexReport::default();
    for (vertex, idx) in candidate_corners(input, edge_geometries) {
        let position = input.vertices.get(vertex).map(|v| v.position);
        let owner_edge = edge_geometries[idx[0]].analysis.edge;
        let decline = |reason: String| (owner_edge, failure(vertex, position, reason));
        let bevels: Vec<Plane> = idx
            .iter()
            .filter_map(|&i| plane_of(&edge_geometries[i].surface))
            .collect();
        if bevels.len() != idx.len() {
            return Err(decline("bevel surfaces are not all planes".into()));
        }
        let apex = concurrent_apex(&bevels, tol).map_err(decline)?;
        let supports: Vec<Plane> = idx
            .iter()
            .flat_map(|&i| {
                let a = &edge_geometries[i].analysis;
                [plane_of(&a.surface_1), plane_of(&a.surface_2)]
            })
            .flatten()
            .collect();
        let mut staged = out.clone();
        super::mitre_apex_scaffold::contract_inset_scaffolds(
            &mut staged,
            vertex,
            &bevels,
            &supports,
            tol,
        )
        .map_err(decline)?;
        let cluster = read_cluster(&staged, vertex, &bevels, tol).map_err(decline)?;
        rebuild_corner(&mut staged, vertex, &cluster, apex, tol).map_err(decline)?;
        if let Err(e) = staged.validate() {
            return Err(decline(format!("rebuilt corner fails validation: {e}")));
        }
        *out = staged;
        report.corners_settled += 1;
    }
    Ok(report)
}
