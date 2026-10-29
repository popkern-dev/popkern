//! Closed rim chamfer chains (PK-135).
//!
//! When every edge of one planar face's outer loop is chamfered — the
//! "cap" of a prism-like body — the per-edge strip pipeline cannot
//! serve the request: each edge's surgery would retrim the shared cap
//! independently and the strips never meet at the corners (measured:
//! four overlapping caps and 93 open mesh edges on a square rim).
//! But this configuration is fully determined by construction: the
//! cap outer loop shrinks to the mitred inset polygon, each rim edge
//! carries one planar bevel quad, and adjacent bevels share a mitre
//! edge. Convex and reflex corners are served uniformly — the mitre
//! is the intersection of adjacent bevel planes, which exists for any
//! corner angle bounded away from 180°.
//!
//! **Invariants**: topology changes go through `mev`/`mef` only (each
//! op verifies its own Euler postconditions); the cap keeps its face
//! handle, surface, and inner loops; wall faces keep their handles and
//! loops (their rim vertices move); rim edges keep their handles and
//! become the wall-side bevel boundaries. **Tolerance**: all
//! coincidence checks use the caller's explicit linear tolerance.
//! **Fallback**: `detect_closed_rim_chains` declines (returns no
//! plan) for anything outside the class — non-planar cap or walls,
//! curved edges, partial loops, vertices with valence ≠ 3, or
//! near-180° corners — and the caller falls through to the per-edge
//! pipeline. A chain whose bevels consume the cap
//! ([`chain_consumes_cap`]) declines whole and grows from a smaller
//! ring (PK-1001); one whose corner setbacks differ (walls meeting the
//! cap at different angles) also falls through to the per-edge path.
//! **Rollback**: the caller applies this on a clone and discards it on
//! `Err`.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::EdgeGeometry;

#[path = "rim_chain_cap.rs"]
mod cap;
pub(crate) use cap::chain_consumes_cap;

/// One edge of a detected chain, in cap-loop traversal order.
struct RimChainEntry {
    /// Index into the caller's `edge_geometries`.
    geometry_index: usize,
    /// The cap-loop half-edge running along this rim edge.
    cap_half_edge: HalfEdgeHandle,
    /// Rim vertex at this half-edge's origin (corner `i`).
    corner_vertex: VertexHandle,
}

/// A complete closed chamfer chain around one cap face.
pub(crate) struct RimChainPlan {
    cap_face: FaceHandle,
    entries: Vec<RimChainEntry>,
}

impl RimChainPlan {
    /// Edge-geometry indices consumed by this chain.
    pub(crate) fn geometry_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.entries.iter().map(|e| e.geometry_index)
    }

    pub(crate) fn cap_face(&self) -> FaceHandle {
        self.cap_face
    }
}

fn face_plane(brep: &BRep, face: FaceHandle) -> Option<()> {
    matches!(
        brep.faces.get(face)?.surface.as_ref()?.geometry,
        SurfaceGeometry::Plane { .. }
    )
    .then_some(())
}

fn edge_is_linear(brep: &BRep, edge: EdgeHandle) -> bool {
    match brep.edges.get(edge).and_then(|e| e.curve.as_ref()) {
        Some(c) => matches!(c.geometry, CurveGeometry::Line { .. }),
        // No stored curve: positions define a segment; acceptable.
        None => true,
    }
}

/// A straight trim curve as an infinite (point, direction) line.
/// Accepts any NURBS whose control points are collinear within `tol`
/// (chamfer trims are often degree-1 polylines with interior points).
fn trim_line(curve: &popkern_core::geometry::NurbsCurve, tol: f64) -> Option<(Point3, Vec3)> {
    let first = *curve.control_points.first()?;
    let last = *curve.control_points.last()?;
    let d = (last - first).normalized()?;
    for p in &curve.control_points {
        let v = *p - first;
        let off = (v - d * v.dot(&d)).length();
        if off > tol {
            return None;
        }
    }
    Some((first, d))
}

/// Closest-approach intersection of two lines; `None` if they are
/// near-parallel or miss by more than `tol`.
fn line_intersection(a: (Point3, Vec3), b: (Point3, Vec3), tol: f64) -> Option<Point3> {
    let (p1, d1) = a;
    let (p2, d2) = b;
    let cross = d1.cross(&d2);
    let denom = cross.length_squared();
    if denom < 1.0e-14 {
        return None;
    }
    let r = p2 - p1;
    let t1 = r.cross(&d2).dot(&cross) / denom;
    let t2 = r.cross(&d1).dot(&cross) / denom;
    let q1 = p1 + d1 * t1;
    let q2 = p2 + d2 * t2;
    (q1.distance_to(&q2) <= tol).then(|| {
        Point3::new(
            (q1.x + q2.x) / 2.0,
            (q1.y + q2.y) / 2.0,
            (q1.z + q2.z) / 2.0,
        )
    })
}

fn outer_loop_half_edges(brep: &BRep, face: FaceHandle) -> Option<Vec<HalfEdgeHandle>> {
    let loop_ = brep.loops.get(brep.faces.get(face)?.outer_loop)?;
    let start = loop_.half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}

/// Detect complete closed chamfer chains. Every returned plan's cap
/// outer loop consists solely of chamfered straight edges between
/// planar faces, with valence-3 linear-edge corners.
pub(crate) fn detect_closed_rim_chains(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
) -> Vec<RimChainPlan> {
    use std::collections::BTreeMap;
    let mut by_edge: BTreeMap<EdgeHandle, usize> = BTreeMap::new();
    for (i, eg) in edge_geometries.iter().enumerate() {
        by_edge.insert(eg.analysis.edge, i);
    }

    // Candidate caps: faces adjacent to chamfered edges.
    let mut candidates: Vec<FaceHandle> = Vec::new();
    for eg in edge_geometries {
        for f in [eg.analysis.face_1, eg.analysis.face_2] {
            if !candidates.contains(&f) {
                candidates.push(f);
            }
        }
    }

    let mut used: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut plans = Vec::new();
    'cap: for cap in candidates {
        if face_plane(brep, cap).is_none() {
            continue;
        }
        let Some(loop_hes) = outer_loop_half_edges(brep, cap) else {
            continue;
        };
        if loop_hes.len() < 3 {
            continue;
        }
        let mut entries = Vec::with_capacity(loop_hes.len());
        for he in &loop_hes {
            let Some(he_data) = brep.half_edges.get(*he) else {
                continue 'cap;
            };
            let Some(&gi) = by_edge.get(&he_data.edge) else {
                continue 'cap; // loop edge not chamfered → not a full chain
            };
            if used.contains(&gi) {
                continue 'cap;
            }
            let eg = &edge_geometries[gi];
            // The chain face must be one side of this chamfer.
            if eg.analysis.face_1 != cap && eg.analysis.face_2 != cap {
                continue 'cap;
            }
            let wall = if eg.analysis.face_1 == cap {
                eg.analysis.face_2
            } else {
                eg.analysis.face_1
            };
            if face_plane(brep, wall).is_none()
                || !matches!(eg.surface, SurfaceGeometry::Plane { .. })
                || !edge_is_linear(brep, he_data.edge)
            {
                continue 'cap;
            }
            entries.push(RimChainEntry {
                geometry_index: gi,
                cap_half_edge: *he,
                corner_vertex: he_data.origin,
            });
        }
        // Corner requirements: valence 3 (two rim edges + one linear
        // corner edge shared by the two walls).
        for e in &entries {
            let Some(edges) = vertex_edges(brep, e.corner_vertex) else {
                continue 'cap;
            };
            if edges.len() != 3 || !edges.iter().all(|&ed| edge_is_linear(brep, ed)) {
                continue 'cap;
            }
        }
        for e in &entries {
            used.insert(e.geometry_index);
        }
        plans.push(RimChainPlan {
            cap_face: cap,
            entries,
        });
    }
    plans
}

fn vertex_edges(brep: &BRep, vertex: VertexHandle) -> Option<Vec<EdgeHandle>> {
    let start = brep.vertices.get(vertex)?.outgoing;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let he = brep.half_edges.get(cursor)?;
        if !out.contains(&he.edge) {
            out.push(he.edge);
        }
        cursor = brep.half_edges.get(he.twin)?.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}

/// The corner edge at `vertex`: the incident edge that is neither of
/// the two rim edges.
fn corner_edge_line(
    brep: &BRep,
    vertex: VertexHandle,
    rim_a: EdgeHandle,
    rim_b: EdgeHandle,
) -> Result<(EdgeHandle, (Point3, Vec3)), String> {
    let edges = vertex_edges(brep, vertex)
        .ok_or_else(|| format!("vertex ring walk failed at {vertex:?}"))?;
    let corner = edges
        .into_iter()
        .find(|&e| e != rim_a && e != rim_b)
        .ok_or_else(|| format!("no corner edge at {vertex:?}"))?;
    let (start_v, end_v) = edge_vertices(brep, corner)
        .ok_or_else(|| format!("corner edge {corner:?} endpoints unavailable"))?;
    let p0 = brep.vertices[start_v].position;
    let p1 = brep.vertices[end_v].position;
    let d = (p1 - p0)
        .normalized()
        .ok_or_else(|| format!("degenerate corner edge at {vertex:?}"))?;
    Ok((corner, (p0, d)))
}

/// Apply one chain. On success returns the new bevel faces (one per
/// rim edge, in chain order). On error the brep may be partially
/// modified — the caller works on a clone and discards it.
pub(crate) fn apply_rim_chain_chamfer(
    brep: &mut BRep,
    plan: &RimChainPlan,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, String> {
    let n = plan.entries.len();
    let lin_tol = tolerance.linear();
    let join_tol = (lin_tol * 10.0).max(1.0e-9);

    // ---- Geometry: inset corners C_i and wall corners W_i. --------
    // Entry i's cap half-edge runs corner i -> corner i+1. The trim
    // on the cap side of edge i and of edge i-1 meet at C_i; the two
    // wall-side trims meet the corner edge at W_i.
    let mut inset_corners = Vec::with_capacity(n); // C_i
    let mut wall_corners = Vec::with_capacity(n); // W_i
    // Corner edges, read while every corner still has valence 3: after
    // the spurs, "the edge that is not a rim edge" is ambiguous (the
    // ring walk can meet the spur first, leaving the corner edge's
    // curve at the vertex's old position).
    let mut corner_edges = Vec::with_capacity(n);
    for i in 0..n {
        let prev = &plan.entries[(i + n - 1) % n];
        let cur = &plan.entries[i];
        let eg_prev = &edge_geometries[prev.geometry_index];
        let eg_cur = &edge_geometries[cur.geometry_index];
        let (cap_prev, wall_prev) = cap_wall_trims(eg_prev, plan.cap_face);
        let (cap_cur, wall_cur) = cap_wall_trims(eg_cur, plan.cap_face);
        let cap_line_prev = trim_line(cap_prev, join_tol)
            .ok_or_else(|| format!("edge {i}: cap trim not a line"))?;
        let cap_line_cur =
            trim_line(cap_cur, join_tol).ok_or_else(|| format!("edge {i}: cap trim not a line"))?;
        let c = line_intersection(cap_line_prev, cap_line_cur, join_tol)
            .ok_or_else(|| format!("corner {i}: cap trims do not meet (near-180° corner?)"))?;

        let (corner_edge, corner_line) = corner_edge_line(
            brep,
            cur.corner_vertex,
            eg_prev.analysis.edge,
            eg_cur.analysis.edge,
        )?;
        corner_edges.push(corner_edge);
        let wall_line_prev = trim_line(wall_prev, join_tol)
            .ok_or_else(|| format!("edge {i}: wall trim not a line"))?;
        let wall_line_cur = trim_line(wall_cur, join_tol)
            .ok_or_else(|| format!("edge {i}: wall trim not a line"))?;
        let w_a = line_intersection(wall_line_prev, corner_line, join_tol)
            .ok_or_else(|| format!("corner {i}: prev wall trim misses the corner edge"))?;
        let w_b = line_intersection(wall_line_cur, corner_line, join_tol)
            .ok_or_else(|| format!("corner {i}: next wall trim misses the corner edge"))?;
        if w_a.distance_to(&w_b) > join_tol {
            return Err(format!(
                "corner {i}: wall setbacks disagree by {:.3e} (non-uniform chain)",
                w_a.distance_to(&w_b)
            ));
        }
        let w = Point3::new(
            (w_a.x + w_b.x) / 2.0,
            (w_a.y + w_b.y) / 2.0,
            (w_a.z + w_b.z) / 2.0,
        );
        if w.distance_to(&c) <= join_tol {
            return Err(format!("corner {i}: degenerate mitre (W == C)"));
        }
        inset_corners.push(c);
        wall_corners.push(w);
    }

    // ---- Topology: n spur mev + n closing mef (inset pattern). ----
    let mut spurs = Vec::with_capacity(n); // (C_i vertex, spur edge, spur back-he C_i -> v_i)
    for (i, entry) in plan.entries.iter().enumerate() {
        let r = brep
            .mev(entry.cap_half_edge, inset_corners[i])
            .map_err(|e| format!("corner {i}: mev failed: {e}"))?;
        spurs.push((r.vertex, r.edge, r.half_edge_twin));
    }
    let mut bevel_faces = Vec::with_capacity(n);
    let mut inset_edges = Vec::with_capacity(n);
    for i in 0..n {
        let he1 = spurs[i].2;
        let next_c = spurs[(i + 1) % n].0;
        let he2 = outer_loop_half_edges(brep, plan.cap_face)
            .ok_or("cap loop walk failed mid-surgery")?
            .into_iter()
            .find(|he| brep.half_edges.get(*he).is_some_and(|h| h.origin == next_c))
            .ok_or_else(|| format!("corner {i}: no cap half-edge at next inset corner"))?;
        let r = brep
            .mef(he1, he2)
            .map_err(|e| format!("bevel {i}: mef failed: {e}"))?;
        bevel_faces.push(r.face);
        inset_edges.push(r.edge);
    }

    // ---- Vertex moves + fresh analytic curves. --------------------
    for (i, entry) in plan.entries.iter().enumerate() {
        brep.vertices[entry.corner_vertex].position = wall_corners[i];
    }
    for (i, entry) in plan.entries.iter().enumerate() {
        let eg = &edge_geometries[entry.geometry_index];
        let w0 = wall_corners[i];
        let w1 = wall_corners[(i + 1) % n];
        set_line_curve(brep, eg.analysis.edge, w0, w1, lin_tol)?;
        set_line_curve(brep, spurs[i].1, w0, inset_corners[i], lin_tol)?;
        set_line_curve(
            brep,
            inset_edges[i],
            inset_corners[i],
            inset_corners[(i + 1) % n],
            lin_tol,
        )?;
        let corner_edge = corner_edges[i];
        let (sv, ev) = edge_vertices(brep, corner_edge)
            .ok_or_else(|| format!("corner edge {corner_edge:?} endpoints unavailable"))?;
        let (p0, p1) = (brep.vertices[sv].position, brep.vertices[ev].position);
        set_line_curve(brep, corner_edge, p0, p1, lin_tol)?;
    }

    // ---- Bevel surfaces (orientation from the loop's Newell normal). ----
    for (i, entry) in plan.entries.iter().enumerate() {
        let eg = &edge_geometries[entry.geometry_index];
        let SurfaceGeometry::Plane { normal, .. } = eg.surface else {
            return Err(format!("bevel {i}: chamfer surface is not a plane"));
        };
        let loop_normal = newell_normal(brep, bevel_faces[i])
            .ok_or_else(|| format!("bevel {i}: degenerate loop normal"))?;
        brep.faces[bevel_faces[i]].surface = Some(SurfaceOnFace {
            geometry: eg.surface.clone(),
            orientation: loop_normal.dot(&normal) >= 0.0,
        });
    }

    // ---- Invalidate stale pcurves on every touched face. ----------
    let mut touched: Vec<FaceHandle> = vec![plan.cap_face];
    touched.extend(bevel_faces.iter().copied());
    for entry in &plan.entries {
        let eg = &edge_geometries[entry.geometry_index];
        for f in [eg.analysis.face_1, eg.analysis.face_2] {
            if !touched.contains(&f) {
                touched.push(f);
            }
        }
    }
    for face in &touched {
        if let Some(hes) = outer_loop_half_edges(brep, *face) {
            for he in hes {
                if let Some(h) = brep.half_edges.get_mut(he) {
                    h.pcurve = None;
                }
            }
        }
    }
    Ok(bevel_faces)
}

fn cap_wall_trims(
    eg: &EdgeGeometry,
    cap: FaceHandle,
) -> (
    &popkern_core::geometry::NurbsCurve,
    &popkern_core::geometry::NurbsCurve,
) {
    if eg.analysis.face_1 == cap {
        (&eg.trim_3d_1, &eg.trim_3d_2)
    } else {
        (&eg.trim_3d_2, &eg.trim_3d_1)
    }
}

fn set_line_curve(
    brep: &mut BRep,
    edge: EdgeHandle,
    start: Point3,
    end: Point3,
    tol: f64,
) -> Result<(), String> {
    let curve = CurveOnEdge::line_segment(start, end, tol)
        .ok_or_else(|| format!("degenerate segment on {edge:?}"))?;
    brep.edges
        .get_mut(edge)
        .ok_or_else(|| format!("stale edge {edge:?}"))?
        .curve = Some(curve);
    Ok(())
}

fn newell_normal(brep: &BRep, face: FaceHandle) -> Option<Vec3> {
    let pts = brep.face_outer_loop_positions(face);
    if pts.len() < 3 {
        return None;
    }
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        nx += (a.y - b.y) * (a.z + b.z);
        ny += (a.z - b.z) * (a.x + b.x);
        nz += (a.x - b.x) * (a.y + b.y);
    }
    Vec3::new(nx, ny, nz).normalized()
}

fn edge_vertices(brep: &BRep, edge: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let next = brep.half_edges.get(half_edge.next)?;
    Some((half_edge.origin, next.origin))
}
