//! Ridge trees of planar half-space corners (PK-848, PK-997).
//!
//! Near a corner where `N` planes bound the material and do not meet in
//! one point, the exact half-space intersection is made of pieces of
//! those planes meeting along a TREE: its interior nodes are the
//! feasible triple points of the planes, its edges the ridges between
//! nodes that share two planes, and its leaves are given points that
//! each lie on a known pair of the planes (where the corner's pieces
//! continue away from it). The mitre apex stage uses it with bevel
//! planes and inset corners as leaves; the truncated apex stage with a
//! vertex fan's face planes and the far ends of its edges.

use super::mitre_apex_corner::{Plane, plane_distance};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};
use std::collections::BTreeSet;

/// A tree node: its position and the planes it lies on.
#[derive(Debug, Clone)]
pub(super) struct RidgeNode {
    pub(super) position: Point3,
    pub(super) planes: BTreeSet<usize>,
    /// The index of the leaf this node is, `None` for a triple point.
    pub(super) leaf: Option<usize>,
}

/// A tree edge: two node indices and the plane pair it lies on.
pub(super) type Ridge = (usize, usize, (usize, usize));

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

/// The ridge tree of a corner at `corner`: `planes` with the sign of
/// each plane's REMOVED side (`(p − o) · n` has that sign off the
/// material), and `leaves` as `(position, plane pair)`. Returns the
/// nodes (leaves first, in the given order) and the ridges, or why the
/// planes do not form a tree with exactly those leaves.
pub(super) fn ridge_tree(
    corner: Point3,
    planes: &[Plane],
    removed_side: &[f64],
    leaves: &[(Point3, (usize, usize))],
    tol: f64,
) -> Result<(Vec<RidgeNode>, Vec<Ridge>), String> {
    let feasible = |p: &Point3| {
        planes
            .iter()
            .zip(removed_side)
            .all(|((o, n), s)| (*p - *o).dot(n) * s <= tol)
    };
    let mut nodes: Vec<RidgeNode> = leaves
        .iter()
        .enumerate()
        .map(|(k, &(position, (i, j)))| RidgeNode {
            position,
            planes: [i, j].into_iter().collect(),
            leaf: Some(k),
        })
        .collect();
    let reach = leaves
        .iter()
        .map(|(p, _)| p.distance_to(&corner))
        .fold(0.0f64, f64::max);
    let n = planes.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let Some(t) = solve3([&planes[i], &planes[j], &planes[k]]) else {
                    continue;
                };
                // A feasible triple point of THIS corner lies within the
                // corner's own reach (the leaves bound its pieces).
                if !feasible(&t) || t.distance_to(&corner) > 2.0 * reach {
                    continue;
                }
                if nodes[..leaves.len()]
                    .iter()
                    .any(|leaf| leaf.position.distance_to(&t) <= tol)
                {
                    return Err(format!(
                        "a triple point coincides with a leaf at {t:?} (a zero-length piece)"
                    ));
                }
                if let Some(node) = nodes
                    .iter_mut()
                    .skip(leaves.len())
                    .find(|node| node.position.distance_to(&t) <= tol)
                {
                    node.planes.extend([i, j, k]);
                } else {
                    nodes.push(RidgeNode {
                        position: t,
                        planes: [i, j, k].into_iter().collect(),
                        leaf: None,
                    });
                }
            }
        }
    }
    // Ridges: along each plane pair's line, consecutive nodes on it.
    let mut edges = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let mut on: Vec<usize> = (0..nodes.len())
                .filter(|&k| nodes[k].planes.contains(&i) && nodes[k].planes.contains(&j))
                .collect();
            if on.len() < 2 {
                continue;
            }
            let dir: Vec3 = planes[i].1.cross(&planes[j].1);
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
    // A tree whose leaves are exactly the given ones.
    if edges.len() + 1 != nodes.len() {
        return Err(format!(
            "ridge graph has {} nodes and {} edges, not a tree",
            nodes.len(),
            edges.len()
        ));
    }
    for (k, node) in nodes.iter().enumerate() {
        let degree = edges.iter().filter(|e| e.0 == k || e.1 == k).count();
        if node.leaf.is_some() != (degree == 1) {
            return Err(format!(
                "ridge tree node at {:?} has degree {degree}",
                node.position
            ));
        }
    }
    Ok((nodes, edges))
}

/// Indices of the planes in `planes` holding `p`.
pub(super) fn planes_through(planes: &[Plane], p: &Point3, tol: f64) -> Vec<usize> {
    (0..planes.len())
        .filter(|&i| plane_distance(&planes[i], p) <= tol)
        .collect()
}

/// Split the fan of `vertex` (outgoing `half_edges` in successor order,
/// edge `k` between planes `k` and `k + 1`) into a ridge tree whose
/// first `half_edges.len()` nodes are the fan's edges as leaves: the
/// first interior node stays `vertex`, every interior ridge on planes
/// `(i, j)` is one `split_vertex` moving the contiguous sector of fan
/// edges between them to the new node at `positions[node]`. Returns
/// each node's vertex (`None` for leaves); the inner `Err` = a ridge
/// that does not separate a contiguous sector (nothing split yet only
/// when it is the first), the outer `Err` = an operator failure.
pub(super) fn split_fan(
    out: &mut BRep,
    half_edges: &[HalfEdgeHandle],
    vertex: VertexHandle,
    positions: &[Point3],
    ridges: &[Ridge],
) -> Result<Result<Vec<Option<VertexHandle>>, String>, String> {
    let m = half_edges.len();
    let Some(root) = (m < positions.len()).then_some(m) else {
        return Ok(Err("the tree has no interior node".into()));
    };
    let mut built: Vec<Option<VertexHandle>> = vec![None; positions.len()];
    built[root] = Some(vertex);
    let mut pending: Vec<Ridge> = ridges
        .iter()
        .copied()
        .filter(|&(a, b, _)| a >= m && b >= m)
        .collect();
    let sector = |from: usize, to: usize| -> BTreeSet<usize> {
        let mut s = BTreeSet::new();
        let mut k = from;
        while k != to {
            s.insert(k);
            k = (k + 1) % m;
        }
        s
    };
    while !pending.is_empty() {
        let Some(pos) = pending
            .iter()
            .position(|e| built[e.0].is_some() || built[e.1].is_some())
        else {
            return Ok(Err("the ridge tree is not connected".into()));
        };
        let (p, q, (i, j)) = pending.remove(pos);
        let (u, w) = if built[p].is_some() { (p, q) } else { (q, p) };
        let side = subtree_leaves(ridges, w, u, m);
        // The child's sector: fan edges i..j or j..i, whichever it owns.
        let (first, last) = if sector(i, j) == side {
            (i, (j + m - 1) % m)
        } else if sector(j, i) == side {
            (j, (i + m - 1) % m)
        } else {
            return Ok(Err(format!(
                "ridge on planes {i}/{j} does not separate a contiguous sector"
            )));
        };
        let split = out
            .split_vertex(half_edges[first], half_edges[last], positions[w])
            .map_err(|e| format!("split_vertex ridge {i}/{j}: {e}"))?;
        built[w] = Some(split.vertex);
    }
    Ok(Ok(built))
}

/// Leaf indices in the subtree of `node` entered from `parent`.
fn subtree_leaves(
    ridges: &[Ridge],
    node: usize,
    parent: usize,
    leaf_count: usize,
) -> BTreeSet<usize> {
    let mut leaves = BTreeSet::new();
    let mut stack = vec![(node, parent)];
    while let Some((at, from)) = stack.pop() {
        if at < leaf_count {
            leaves.insert(at);
            continue;
        }
        for &(a, b, _) in ridges {
            let next = if a == at {
                b
            } else if b == at {
                a
            } else {
                continue;
            };
            if next != from {
                stack.push((next, at));
            }
        }
    }
    leaves
}
