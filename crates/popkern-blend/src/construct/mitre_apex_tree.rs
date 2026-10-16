//! Mitre trees for non-concurrent planar chamfer corners (PK-848).
//!
//! When the `N` bevel planes at an all-chamfered corner do not meet in
//! one point, the exact half-space intersection near the corner is
//! still made of bevel faces only, but their corner pieces meet along a
//! TREE of mitre edges: its interior nodes are the feasible triple
//! points of the bevel planes, its edges the ridges between nodes that
//! share two bevels, and its leaves the `N` inset corners (each on two
//! bevels). A stretched octahedron's four bevels at a short vertex meet
//! in two triple points joined by a ridge; a single concurrent apex is
//! the one-node star `mitre_apex_corner` builds.
//!
//! **Construction**, on the collapsed `2N`-gon of the served cluster
//! (anchors alternating with inset corners): every tree edge lies on two
//! bevels `i`, `j`, so it separates the cells of bevels `i` and `j`, and
//! is inserted (`mev` to a new node, else `mef` between built ones) in
//! the loop that still holds both anchors `a_i`, `a_j`. Edges go in
//! breadth-first order from the inset leaves. Each finished cell holds
//! one anchor and is merged into its bevel by the caller's existing
//! anchor pass.

use super::mitre_apex_corner::{Plane, outgoing, plane_distance};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};
use std::collections::BTreeSet;

/// A tree node: its position and the bevels it lies on.
#[derive(Debug, Clone)]
struct Node {
    position: Point3,
    bevels: BTreeSet<usize>,
    /// The cluster vertex of an inset leaf.
    vertex: Option<VertexHandle>,
}

fn solve3(planes: [&Plane; 3]) -> Option<Point3> {
    let rows: Vec<([f64; 3], f64)> = planes
        .iter()
        .map(|(o, n)| ([n.x, n.y, n.z], n.x * o.x + n.y * o.y + n.z * o.z))
        .collect();
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let m = [rows[0].0, rows[1].0, rows[2].0];
    let d = det(m);
    // Unit normals: a determinant this small is a common line.
    if d.abs() < 1e-9 {
        return None;
    }
    let mut x = [0.0; 3];
    for (col, xc) in x.iter_mut().enumerate() {
        let mut mc = m;
        for row in 0..3 {
            mc[row][col] = rows[row].1;
        }
        *xc = det(mc) / d;
    }
    Some(Point3::new(x[0], x[1], x[2]))
}

/// The mitre tree of a corner: `(nodes, edges as node index pairs,
/// the bevel pair each edge lies on)`, or why there is none.
#[allow(clippy::type_complexity)]
fn mitre_tree(
    corner: Point3,
    bevels: &[Plane],
    insets: &[(VertexHandle, Point3)],
    tol: f64,
) -> Result<(Vec<Node>, Vec<(usize, usize, (usize, usize))>), String> {
    // Material lies on the far side of every bevel from the cut corner.
    let side: Vec<f64> = bevels
        .iter()
        .map(|(o, n)| (corner - *o).dot(n).signum())
        .collect();
    let feasible = |p: &Point3| {
        bevels
            .iter()
            .zip(&side)
            .all(|((o, n), s)| (*p - *o).dot(n) * s <= tol)
    };
    let mut nodes: Vec<Node> = Vec::new();
    for &(vertex, position) in insets {
        let on: BTreeSet<usize> = (0..bevels.len())
            .filter(|&k| plane_distance(&bevels[k], &position) <= tol)
            .collect();
        if on.len() != 2 {
            return Err(format!(
                "inset corner {vertex:?} lies on {} bevels, expected 2",
                on.len()
            ));
        }
        nodes.push(Node {
            position,
            bevels: on,
            vertex: Some(vertex),
        });
    }
    let reach = insets
        .iter()
        .map(|(_, p)| p.distance_to(&corner))
        .fold(0.0f64, f64::max);
    let n = bevels.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let Some(t) = solve3([&bevels[i], &bevels[j], &bevels[k]]) else {
                    continue;
                };
                // A feasible triple point of THIS corner lies within the
                // corner's own reach (the insets bound its bevel pieces).
                if !feasible(&t) || t.distance_to(&corner) > 2.0 * reach {
                    continue;
                }
                if let Some(node) = nodes
                    .iter_mut()
                    .skip(insets.len())
                    .find(|node| node.position.distance_to(&t) <= tol)
                {
                    node.bevels.extend([i, j, k]);
                } else {
                    nodes.push(Node {
                        position: t,
                        bevels: [i, j, k].into_iter().collect(),
                        vertex: None,
                    });
                }
            }
        }
    }
    // Ridges: along each bevel pair's line, consecutive nodes on it.
    let mut edges = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let mut on: Vec<usize> = (0..nodes.len())
                .filter(|&k| nodes[k].bevels.contains(&i) && nodes[k].bevels.contains(&j))
                .collect();
            if on.len() < 2 {
                continue;
            }
            let dir: Vec3 = bevels[i].1.cross(&bevels[j].1);
            let origin = nodes[on[0]].position;
            on.sort_by(|&a, &b| {
                let ta = (nodes[a].position - origin).dot(&dir);
                let tb = (nodes[b].position - origin).dot(&dir);
                ta.total_cmp(&tb)
            });
            for w in on.windows(2) {
                edges.push((w[0], w[1], (i, j)));
            }
        }
    }
    // A tree whose leaves are exactly the insets.
    if edges.len() + 1 != nodes.len() {
        return Err(format!(
            "mitre graph has {} nodes and {} edges, not a tree",
            nodes.len(),
            edges.len()
        ));
    }
    for (k, node) in nodes.iter().enumerate() {
        let degree = edges.iter().filter(|e| e.0 == k || e.1 == k).count();
        if node.vertex.is_some() != (degree == 1) {
            return Err(format!(
                "mitre tree node at {:?} has degree {degree}",
                node.position
            ));
        }
    }
    Ok((nodes, edges))
}

fn loop_vertices(brep: &BRep, loop_: LoopHandle) -> Option<Vec<(VertexHandle, HalfEdgeHandle)>> {
    let start = brep.loops.get(loop_)?.half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let he = brep.half_edges.get(cursor)?;
        out.push((he.origin, cursor));
        cursor = he.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}

/// Insert the mitre tree of a non-concurrent corner into the collapsed
/// `2N`-gon `polygon` (anchors alternating with insets). Returns the
/// tree's edges (for exact line curves) and its interior vertices.
pub(super) fn build_mitre_tree(
    out: &mut BRep,
    polygon: LoopHandle,
    corner: Point3,
    bevels: &[Plane],
    anchors: &[VertexHandle],
    insets: &[VertexHandle],
    tol: f64,
) -> Result<(Vec<EdgeHandle>, Vec<VertexHandle>), String> {
    let inset_positions: Vec<(VertexHandle, Point3)> = insets
        .iter()
        .map(|&c| {
            out.vertices
                .get(c)
                .map(|v| (c, v.position))
                .ok_or_else(|| format!("stale inset {c:?}"))
        })
        .collect::<Result<_, _>>()?;
    let (nodes, edges) = mitre_tree(corner, bevels, &inset_positions, tol)?;
    // Each anchor's bevel: the one its two polygon neighbours share.
    let ring = loop_vertices(out, polygon).ok_or("stale polygon")?;
    let at = |v: VertexHandle| ring.iter().position(|(r, _)| *r == v);
    let bevel_of_inset = |v: VertexHandle| {
        nodes
            .iter()
            .find(|node| node.vertex == Some(v))
            .map(|node| node.bevels.clone())
    };
    let mut anchor_of_bevel = vec![None; bevels.len()];
    for &a in anchors {
        let k = at(a).ok_or_else(|| format!("anchor {a:?} not on the polygon"))?;
        let prev = ring[(k + ring.len() - 1) % ring.len()].0;
        let next = ring[(k + 1) % ring.len()].0;
        let (Some(p), Some(q)) = (bevel_of_inset(prev), bevel_of_inset(next)) else {
            return Err(format!("anchor {a:?} not between two inset corners"));
        };
        let shared: Vec<usize> = p.intersection(&q).copied().collect();
        let [b] = shared.as_slice() else {
            return Err(format!(
                "anchor {a:?}'s insets share {} bevels",
                shared.len()
            ));
        };
        anchor_of_bevel[*b] = Some(a);
    }
    let mut built: Vec<Option<VertexHandle>> = nodes.iter().map(|n| n.vertex).collect();
    let mut pending: Vec<(usize, usize, (usize, usize))> = edges;
    let mut tree_edges = Vec::new();
    while !pending.is_empty() {
        let Some(pos) = pending
            .iter()
            .position(|e| built[e.0].is_some() || built[e.1].is_some())
        else {
            return Err("mitre tree is not connected to the polygon".into());
        };
        let (p, q, (i, j)) = pending.remove(pos);
        let (u, w) = if built[p].is_some() { (p, q) } else { (q, p) };
        let (Some(ai), Some(aj)) = (anchor_of_bevel[i], anchor_of_bevel[j]) else {
            return Err(format!("bevels {i}/{j} own no anchor"));
        };
        let uv = built[u].ok_or("unbuilt origin")?;
        // The loop holding u and both anchors is the face this edge splits.
        let loop_ = outgoing(out, uv)
            .into_iter()
            .filter_map(|he| out.half_edges.get(he).map(|h| h.loop_))
            .find(|&l| {
                loop_vertices(out, l).is_some_and(|vs| {
                    vs.iter().any(|(v, _)| *v == ai) && vs.iter().any(|(v, _)| *v == aj)
                })
            })
            .ok_or_else(|| format!("no loop holds node {u} and anchors {ai:?}, {aj:?}"))?;
        let vs = loop_vertices(out, loop_).ok_or("stale loop")?;
        let he_u = vs
            .iter()
            .find(|(v, _)| *v == uv)
            .map(|(_, h)| *h)
            .ok_or("origin off its loop")?;
        if let Some(wv) = built[w] {
            let he_w = vs
                .iter()
                .find(|(v, _)| *v == wv)
                .map(|(_, h)| *h)
                .ok_or_else(|| format!("node {w} off the loop it closes"))?;
            let mef = out
                .mef(he_u, he_w)
                .map_err(|e| format!("mef mitre ridge: {e}"))?;
            tree_edges.push(mef.edge);
        } else {
            let mev = out
                .mev(he_u, nodes[w].position)
                .map_err(|e| format!("mev mitre node: {e}"))?;
            built[w] = Some(mev.vertex);
            tree_edges.push(mev.edge);
        }
    }
    let interior: Vec<VertexHandle> = nodes
        .iter()
        .zip(&built)
        .filter(|(n, _)| n.vertex.is_none())
        .filter_map(|(_, b)| *b)
        .collect();
    Ok((tree_edges, interior))
}
