//! The topology events of bevel growth (PK-1001).
//!
//! At the instant an edge's length reaches zero, the vanishing piece is
//! merged into one vertex: a face of that edge with no area left is
//! consumed (`kef` into the neighbour across its longest side, then
//! `kve` of its vanished sides; a triangle leaves one vertex, a
//! sliver's two ends leave two), otherwise the edge must vanish alone
//! and is contracted (`kve`). Every merged vertex on more than three
//! planes is then split into the ridge tree its planes make just after
//! the instant ([`expand`]). An edge flip is the four-plane case: two
//! nodes joined by a ridge between the two planes the old edge did not
//! lie on.

use super::super::mitre_apex_corner::{Plane, destination, outgoing};
use super::super::mitre_apex_scaffold::three_plane_point;
use super::super::plane_ridge_tree::{ridge_tree, split_fan};
use super::{BevelGrowthReport, Motion, edge_ends, signed_length};
use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};
use std::collections::BTreeSet;

/// Resolve the event at `edge` (vertices already placed at `tau`).
pub(super) fn resolve(
    out: &mut BRep,
    motion: &Motion,
    edge: EdgeHandle,
    tau: f64,
    sigma: f64,
    gate: f64,
    report: &mut BevelGrowthReport,
) -> Result<(), String> {
    let (p, q, _, _) = edge_ends(out, edge).ok_or("stale edge")?;
    let merged = if let Some(face) = [p, q].into_iter().find(|&f| degenerate(out, f, gate)) {
        let merged = consume(out, face, gate)?;
        report.consumed_faces.push(face);
        merged
    } else {
        let (vertices, edges) = vanishing_cluster(out, edge, gate)?;
        if (vertices.len(), edges.len()) != (2, 1) {
            return Err(format!(
                "{} edges over {} vertices vanish together with no face consumed",
                edges.len(),
                vertices.len()
            ));
        }
        let he = out.edges.get(edge).ok_or("stale edge")?.half_edge;
        report.edge_flips += 1;
        vec![
            out.kve(he)
                .map_err(|e| format!("contract kve: {e}"))?
                .vertex,
        ]
    };
    for v in merged {
        if out.vertices.contains_key(v) {
            expand(out, motion, v, tau, sigma, gate)?;
        }
    }
    Ok(())
}

/// Whether `face` has no area left: every corner within `gate` of the
/// line through its two farthest corners (a sliver), or those two
/// within `gate` of each other (a point).
fn degenerate(out: &BRep, face: FaceHandle, gate: f64) -> bool {
    let ring = out.face_outer_loop_positions(face);
    let mut far = (0, 0, 0.0f64);
    for i in 0..ring.len() {
        for j in i + 1..ring.len() {
            let d = ring[i].distance_to(&ring[j]);
            if d > far.2 {
                far = (i, j, d);
            }
        }
    }
    let Some(axis) = (ring.len() >= 2)
        .then(|| (ring[far.1] - ring[far.0]).normalized())
        .flatten()
    else {
        return true;
    };
    if far.2 <= gate {
        return true;
    }
    ring.iter().all(|p| {
        let w = *p - ring[far.0];
        (w - axis * w.dot(&axis)).length() <= gate
    })
}

/// The vertices and edges of zero length connected to `edge`.
fn vanishing_cluster(
    out: &BRep,
    edge: EdgeHandle,
    gate: f64,
) -> Result<(BTreeSet<VertexHandle>, BTreeSet<EdgeHandle>), String> {
    let (_, _, o, d) = edge_ends(out, edge).ok_or("stale edge")?;
    let mut vertices: BTreeSet<VertexHandle> = [o, d].into_iter().collect();
    let mut edges: BTreeSet<EdgeHandle> = [edge].into_iter().collect();
    let mut frontier = vec![o, d];
    while let Some(at) = frontier.pop() {
        let here = out.vertices.get(at).ok_or("stale vertex")?.position;
        for he in outgoing(out, at) {
            let Some(next) = destination(out, he) else {
                continue;
            };
            let there = out.vertices.get(next).ok_or("stale vertex")?.position;
            if here.distance_to(&there) > gate {
                continue;
            }
            edges.insert(out.half_edges.get(he).ok_or("stale half-edge")?.edge);
            if vertices.insert(next) {
                frontier.push(next);
            }
        }
    }
    Ok((vertices, edges))
}

/// The face left of `he`.
fn face_of(out: &BRep, he: HalfEdgeHandle) -> Result<FaceHandle, String> {
    let h = out.half_edges.get(he).ok_or("stale half-edge")?;
    Ok(out.loops.get(h.loop_).ok_or("stale loop")?.face)
}

/// The fan of `vertex` in successor order (`twin.next`).
fn fan(out: &BRep, vertex: VertexHandle) -> Result<Vec<HalfEdgeHandle>, String> {
    let first = *outgoing(out, vertex).first().ok_or("an isolated vertex")?;
    let mut fan = vec![first];
    loop {
        let last = *fan.last().ok_or("empty fan")?;
        let twin = out.half_edges.get(last).ok_or("stale half-edge")?.twin;
        let next = out.half_edges.get(twin).ok_or("stale twin")?.next;
        if next == first {
            return Ok(fan);
        }
        if fan.len() > outgoing(out, vertex).len() {
            return Err("the fan does not close".into());
        }
        fan.push(next);
    }
}

/// Consume `face`, which has shrunk to a point or a sliver: merge it
/// into the neighbour across its longest side (`kef`), then contract
/// its vanished sides (`kve`). A sliver's remaining long side becomes
/// the ridge between the neighbours on either side of it. Returns the
/// merged vertices.
fn consume(out: &mut BRep, face: FaceHandle, gate: f64) -> Result<Vec<VertexHandle>, String> {
    let f = out.faces.get(face).ok_or("stale face")?;
    if !f.inner_loops.is_empty() {
        return Err("a consumed face has a hole".into());
    }
    let start = out.loops.get(f.outer_loop).ok_or("stale loop")?.half_edge;
    let ring = loop_half_edges(out, start)?;
    let mut longest = (ring[0], f64::NEG_INFINITY);
    let mut vanished = Vec::new();
    for &h in &ring {
        let o = out.half_edges.get(h).ok_or("stale half-edge")?.origin;
        let d = destination(out, h).ok_or("an edge without a destination")?;
        let l = out
            .vertices
            .get(o)
            .ok_or("stale vertex")?
            .position
            .distance_to(&out.vertices.get(d).ok_or("stale vertex")?.position);
        if l > longest.1 {
            longest = (h, l);
        }
        if l <= gate {
            vanished.push(h);
        }
    }
    vanished.retain(|&h| h != longest.0);
    let neighbour_side = out.half_edges.get(longest.0).ok_or("stale half-edge")?.twin;
    out.kef(neighbour_side)
        .map_err(|e| format!("consume kef: {e}"))?;
    let mut merged = Vec::new();
    for h in vanished {
        merged.push(out.kve(h).map_err(|e| format!("consume kve: {e}"))?.vertex);
    }
    merged.dedup();
    Ok(merged)
}

/// The half-edges of the loop of `he`, from `he`.
fn loop_half_edges(out: &BRep, he: HalfEdgeHandle) -> Result<Vec<HalfEdgeHandle>, String> {
    let mut ring = vec![he];
    loop {
        let last = *ring.last().ok_or("empty loop")?;
        let next = out.half_edges.get(last).ok_or("stale half-edge")?.next;
        if next == he {
            return Ok(ring);
        }
        if ring.len() > out.half_edges.len() {
            return Err("a loop does not close".into());
        }
        ring.push(next);
    }
}

/// Split `vertex` (at `tau`, on more than three planes) into the ridge
/// tree its planes make just after `tau`. That tree is read from the
/// planes' VELOCITIES: a point leaving the vertex at velocity `w` stays
/// on plane `k` when `n_k · w = ċ_k`, so the result's tangent cone at
/// the vertex is the half-space intersection of the velocity planes,
/// with the fan's edges as leaves far out along their directions. Exact
/// at the instant, with no small step to choose. Planes that keep
/// concurring give one node and no split.
pub(super) fn expand(
    out: &mut BRep,
    motion: &Motion,
    vertex: VertexHandle,
    tau: f64,
    sigma: f64,
    gate: f64,
) -> Result<(), String> {
    let fan = fan(out, vertex)?;
    let m = fan.len();
    let faces: Vec<FaceHandle> = fan
        .iter()
        .map(|&h| face_of(out, h))
        .collect::<Result<_, _>>()?;
    if faces.iter().collect::<BTreeSet<_>>().len() != m {
        return Err("a face appears twice around a merged vertex".into());
    }
    if m <= 3 {
        return Ok(());
    }
    let corner = out.vertices.get(vertex).ok_or("stale vertex")?.position;
    let origin = Point3::new(0.0, 0.0, 0.0);
    let mut planes: Vec<Plane> = Vec::with_capacity(m);
    let mut velocity: Vec<Plane> = Vec::with_capacity(m);
    for &f in &faces {
        let (o, n) = motion
            .plane(f, tau)
            .ok_or("a face at a merged vertex is not planar")?;
        planes.push((o, n));
        velocity.push((origin + n * motion.rate(f).unwrap_or(0.0), n));
    }
    let mut fars = Vec::with_capacity(m);
    let mut dirs = Vec::with_capacity(m);
    for &h in &fan {
        let far = destination(out, h)
            .and_then(|d| out.vertices.get(d))
            .ok_or("stale far end")?
            .position;
        dirs.push(
            (far - corner)
                .normalized()
                .ok_or("a vanished edge is left at a merged vertex")?,
        );
        fars.push(far);
    }
    // The removed side of each plane, from the far ends not on it (all
    // on the material side of a convex fan).
    let mut removed = Vec::with_capacity(m);
    for (k, (o, n)) in planes.iter().enumerate() {
        let mut sign = 0.0f64;
        for (j, far) in fars.iter().enumerate() {
            if j == k || (j + 1) % m == k {
                continue;
            }
            let s = (*far - *o).dot(n);
            if s.abs() <= gate || (sign != 0.0 && s.signum() != sign) {
                return Err("the fan at a merged vertex is not convex".into());
            }
            sign = s.signum();
        }
        removed.push(-sign);
    }
    // Leaves: far out along each edge's line in velocity space.
    let mut base = Vec::with_capacity(m);
    for k in 0..m {
        let cut = (origin, dirs[k]);
        base.push(
            three_plane_point(&velocity[k], &velocity[(k + 1) % m], &cut)
                .ok_or("a fan edge's velocity line is undefined")?,
        );
    }
    let scale = base
        .iter()
        .map(|b| (*b - origin).length())
        .chain(faces.iter().map(|&f| motion.rate(f).unwrap_or(0.0).abs()))
        .fold(0.0f64, f64::max);
    if scale <= 0.0 {
        return Err("a merged vertex does not move".into());
    }
    let leaves: Vec<(Point3, (usize, usize))> = (0..m)
        .map(|k| (base[k] + dirs[k] * (1e3 * scale), (k, (k + 1) % m)))
        .collect();
    let (nodes, ridges) = ridge_tree(origin, &velocity, &removed, &leaves, 1e-9 * scale)?;
    let positions = vec![corner; nodes.len()];
    let built = split_fan(out, &fan, vertex, &positions, &ridges)??;
    // Every new ridge must grow after the instant.
    for &(a, b, _) in ridges.iter().filter(|&&(a, b, _)| a >= m && b >= m) {
        let (Some(u), Some(w)) = (built[a], built[b]) else {
            return Err("a tree node was never built".into());
        };
        let edge = outgoing(out, u)
            .into_iter()
            .find(|&he| destination(out, he) == Some(w))
            .and_then(|he| out.half_edges.get(he).map(|h| h.edge))
            .ok_or("a ridge edge is missing after the split")?;
        if signed_length(out, motion, edge, 1.0, sigma)? <= gate {
            return Err(format!("a ridge made at τ = {tau:.6} does not grow"));
        }
    }
    Ok(())
}
