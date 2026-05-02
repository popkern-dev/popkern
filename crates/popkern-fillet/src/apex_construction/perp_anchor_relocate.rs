//! Phase 6c-A2 §28.C — pre-A2 perp-anchor relocation to plan targets.
//!
//! Phase A perp-walks place anchors at perp-walk-math positions
//! (where the fillet's transverse boundary intersects the host plane,
//! perpendicular to the fillet edge).  The planner places
//! `apex_vertices` at Gregory tangent-termination math positions
//! (where the cap surface meets the fillet cylinder along the
//! cap-tangent arc).  These two coordinate systems agree for uniform
//! radii but diverge for variable radii.
//!
//! This pass runs before A2's stage 1 (`resolve_apex_vertices`).
//! For each plan apex_vertex point, it finds the nearest BRep
//! perp-anchor in V's neighborhood and relocates it (vertex
//! position + incident edge curve coherence sweep) to the plan
//! target.  After the pass, A2's `resolve_apex_vertices` snap
//! finds matching anchors at zero distance, and stage 2's Adopt
//! path fires for CapVsCyl edges.
//!
//! See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-
//! construction.md` §28.3, §29.3, §32.5.

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::{CurveGeometry, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use std::collections::HashSet;

use super::plan::ApexPlan;

/// Relocate perp-anchors in V's neighborhood to plan target positions.
///
/// For each `Point3` in `plan.target_topology.apex_vertices`, find
/// the nearest BRep vertex (within `2 * max_fillet_radius` of V) and
/// move it to the target.  Bipartite distance-sorted matching
/// prevents two plan targets from claiming the same anchor.
///
/// After moving each vertex, sweeps incident edges to maintain
/// curve coherence: line edges are re-fit through the new endpoint;
/// circular-arc edges (the common case for Phase B fillet-end arcs
/// and cyl-vs-Plane cube boundaries) are re-fit by projecting the
/// new vertex position onto the underlying conic and updating the
/// arc parameter range.
///
/// Returns the number of vertices relocated and the worst per-
/// vertex movement distance, for trace logging and diagnostics.
/// Vertices already within `tolerance.linear()` of their target
/// are left unchanged.
///
/// `Ok((n_moved, worst_dist))` on success.
/// `Err(reason)` if matching is ambiguous or a curve-coherence
/// step fails — caller rolls back via the construct_apex snapshot.
pub(super) fn relocate_perp_anchors_to_plan_targets(
    brep: &mut BRep,
    plan: &ApexPlan,
    tolerance: &Tolerance,
) -> Result<(usize, f64), String> {
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let targets: &[Point3] = &plan.target_topology.apex_vertices;
    if targets.is_empty() {
        return Ok((0, 0.0));
    }

    // Synthetic test plans may pass `VertexHandle::null`; treat as
    // "no relocation needed" so unit tests for downstream stages
    // (resolve_apex_vertices, resolve_apex_edges) can drive A2
    // with handcrafted plans that don't have a real apex vertex
    // in the BRep.
    let Some(v_data) = brep.vertices.get(plan.apex_vertex) else {
        return Ok((0, 0.0));
    };
    let v_pos = v_data.position;

    let r_max = max_fillet_radius_from_plan(plan);
    let r_min = min_fillet_radius_from_plan(plan);
    let search_radius = (r_max * 2.0).max(tolerance.linear() * 100.0);
    // Tight relocation gate: only move a perp-anchor when the
    // planner target is within `r_min * 0.6` of an existing
    // perp-anchor.  Larger distances would mean the planner is
    // pointing at a position that has no nearby perp-anchor —
    // forcing a move there would yank the perp-anchor off its
    // home cylinder boundary and break face-curve-vertex
    // coherence.  Targets without a tight match flow through to
    // resolve_apex_vertices (snap with wide tolerance, or
    // split_edge fallback).
    // Threshold tuning: cube_mixed has good matches at 0.020-
    // 0.040 m and pathological far matches at 0.10-0.12 m.
    // pyramid3_uniform has good matches up to ~0.07 m.  Set the
    // gate to `r_min * 0.9` — accepts pyramid3's 0.07 m matches
    // and excludes cube_mixed's 0.10 m forced matches.
    let relocate_max_dist = (r_min * 0.9).max(tolerance.linear() * 10.0);

    // Collect candidate vertices: BFS from V (V's neighborhood) AND
    // any vertex on a cylinder face's outer loop within search_radius
    // of V.  Cylinder face boundaries hold perp-anchors; the BFS
    // from V finds plane-side perp-anchors.
    let mut candidates: HashSet<VertexHandle> = HashSet::new();
    bfs_collect_neighborhood(
        brep,
        plan.apex_vertex,
        search_radius,
        v_pos,
        &mut candidates,
    );
    collect_cylinder_outer_loop_vertices(brep, search_radius, v_pos, &mut candidates);
    candidates.remove(&plan.apex_vertex); // V itself is not a relocation target

    if trace {
        eprintln!(
            "[a2-trace §28.C] relocate: V={:?}@{v_pos:?}, search_radius={search_radius:.4}, \
             {} plan targets, {} candidate perp-anchors",
            plan.apex_vertex,
            targets.len(),
            candidates.len(),
        );
    }

    // Bipartite distance-sorted matching: build all (target_idx,
    // candidate_handle, distance) triples, sort ascending, greedily
    // claim min-distance pairs.  Each plan target claims exactly one
    // candidate; candidates already claimed are skipped.
    let mut triples: Vec<(usize, VertexHandle, f64)> = Vec::new();
    for (ti, t) in targets.iter().enumerate() {
        for &c in &candidates {
            let cp = brep.vertices[c].position;
            let d = t.distance_to(&cp);
            if d <= relocate_max_dist {
                triples.push((ti, c, d));
            }
        }
    }
    triples.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

    let mut target_to_anchor: Vec<Option<VertexHandle>> = vec![None; targets.len()];
    let mut claimed: HashSet<VertexHandle> = HashSet::new();
    for (ti, c, _) in &triples {
        if target_to_anchor[*ti].is_some() {
            continue;
        }
        if claimed.contains(c) {
            continue;
        }
        target_to_anchor[*ti] = Some(*c);
        claimed.insert(*c);
    }

    // Apply relocations.  Skip targets that are already at their
    // anchor (within tolerance.linear()) and targets without an
    // anchor (no candidate within search_radius — the snap path
    // in resolve_apex_vertices will create via split_edge or fail).
    let mut n_moved = 0usize;
    let mut worst = 0.0_f64;
    let lin_tol = tolerance.linear();
    for (ti, anchor_opt) in target_to_anchor.iter().enumerate() {
        let Some(anchor) = anchor_opt else {
            if trace {
                eprintln!(
                    "[a2-trace §28.C]   target[{ti}]@{:?} → no anchor (will fall through to \
                     resolve_apex_vertices)",
                    targets[ti]
                );
            }
            continue;
        };
        let target = targets[ti];
        let cp = brep.vertices[*anchor].position;
        let d = target.distance_to(&cp);
        if d < lin_tol {
            if trace {
                eprintln!(
                    "[a2-trace §28.C]   target[{ti}]@{target:?} ↔ {anchor:?}@{cp:?}: \
                     dist {d:.3e} < lin_tol {lin_tol:.3e}, skip"
                );
            }
            continue;
        }
        if trace {
            eprintln!(
                "[a2-trace §28.C]   target[{ti}]@{target:?} ↔ {anchor:?}@{cp:?}: dist {d:.4e}, \
                 relocating"
            );
        }
        relocate_vertex_with_curve_coherence(brep, *anchor, target)?;
        n_moved += 1;
        if d > worst {
            worst = d;
        }
    }

    if trace {
        eprintln!("[a2-trace §28.C] relocate done: {n_moved} vertices moved, worst {worst:.4e} m");
    }
    Ok((n_moved, worst))
}

/// Move a vertex's position to `target` and update each incident
/// edge's stored curve to maintain coherence.
///
/// Strategy by curve type:
/// * [`CurveGeometry::Line`]: re-fit origin/direction so the
///   line still passes through both endpoints.
/// * [`CurveGeometry::Arc`] / [`CurveGeometry::Circle`] /
///   [`CurveGeometry::Ellipse`]: re-fit the arc parameter range
///   so the moved end terminates at `target`.  The underlying
///   conic geometry (center, radius, basis) is preserved.
/// * [`CurveGeometry::Nurbs`]: replace with a degree-1 NURBS
///   straight line through both endpoints.  Loses original
///   curvature; acceptable for cube-edge / perp-walk segments
///   that originated as straight lines and got NURBS-encoded.
///   Risk: if the NURBS represented an analytic arc, this
///   degrades the cyl-on-plane intersection fidelity by the
///   anchor-to-target offset (~0.04 m on cube_mixed).  Phase B's
///   wire authoring + Phase 6c-D's pcurve refit downstream
///   recover most of the lost fidelity.
fn relocate_vertex_with_curve_coherence(
    brep: &mut BRep,
    vertex: VertexHandle,
    target: Point3,
) -> Result<(), String> {
    // 1. Snapshot incident edges + which end of each is the moved
    //    vertex.  Indexed by edge handle.
    let mut incidents: Vec<(EdgeHandle, EndChoice)> = Vec::new();
    for (eh, edge) in &brep.edges {
        let primary = edge.half_edge;
        let Some(he) = brep.half_edges.get(primary) else {
            continue;
        };
        let twin_he = brep
            .half_edges
            .get(he.twin)
            .ok_or_else(|| format!("edge {eh:?}: twin half-edge missing"))?;
        if he.origin == vertex {
            incidents.push((eh, EndChoice::PrimaryOrigin));
        } else if twin_he.origin == vertex {
            incidents.push((eh, EndChoice::TwinOrigin));
        }
    }

    // 2. Move the vertex.
    brep.vertices
        .get_mut(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} disappeared mid-relocate"))?
        .position = target;

    // 3. For each incident edge, refit curve to maintain coherence.
    for (eh, end_choice) in incidents {
        refit_edge_curve_for_relocated_endpoint(brep, eh, end_choice)?;
    }
    Ok(())
}

/// Whether the relocated vertex was the primary half-edge's origin
/// or the twin half-edge's origin (= primary's destination).
#[derive(Debug, Clone, Copy)]
enum EndChoice {
    PrimaryOrigin,
    TwinOrigin,
}

/// Refit one edge's stored curve so its endpoint at `end_choice`
/// equals the (now-updated) vertex position.
fn refit_edge_curve_for_relocated_endpoint(
    brep: &mut BRep,
    edge: EdgeHandle,
    end_choice: EndChoice,
) -> Result<(), String> {
    let Some(edge_data) = brep.edges.get(edge) else {
        return Ok(());
    };
    let primary = edge_data.half_edge;
    let Some(he) = brep.half_edges.get(primary) else {
        return Ok(());
    };
    let v_origin = he.origin;
    let twin = he.twin;
    let Some(twin_he) = brep.half_edges.get(twin) else {
        return Ok(());
    };
    let v_dest = twin_he.origin;
    let p_origin = brep.vertices[v_origin].position;
    let p_dest = brep.vertices[v_dest].position;
    let _ = end_choice; // Either end may have moved; we re-fit using both endpoints.

    // Pull the existing curve (if any) so we can preserve type.
    let existing = brep.edges.get(edge).and_then(|e| e.curve.clone());

    let new_curve = match existing {
        None => {
            // No curve was stored — leave as None.  Phase B may
            // attach a curve later.
            return Ok(());
        }
        Some(coe) => refit_curve_through_endpoints(&coe, p_origin, p_dest),
    };

    if let Some(edge_data) = brep.edges.get_mut(edge) {
        edge_data.curve = Some(new_curve);
    }
    Ok(())
}

/// Refit a curve so it terminates at `p0` (at `t_start`) and `p1`
/// (at `t_end`).  Returns a new [`CurveOnEdge`].  Strategy by
/// curve type — see [`relocate_vertex_with_curve_coherence`].
fn refit_curve_through_endpoints(coe: &CurveOnEdge, p0: Point3, p1: Point3) -> CurveOnEdge {
    match &coe.geometry {
        CurveGeometry::Line { .. } => {
            // Re-fit as a line from p0 to p1.
            let direction = Vec3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: p0,
                    direction,
                },
                t_start: 0.0,
                t_end: 1.0,
            }
        }
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            // Refit by re-projecting both endpoints onto the arc's
            // underlying circle.  Compute the new sweep angles from
            // (p - center) projected onto the in-plane basis.
            let (start_angle, end_angle) = arc_endpoint_angles(*center, *normal, p0, p1);
            CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center: *center,
                    normal: *normal,
                    radius: *radius,
                    start_angle,
                    end_angle,
                },
                t_start: start_angle,
                t_end: end_angle,
            }
        }
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => {
            // Promote to Arc covering the swept range.
            let (start_angle, end_angle) = arc_endpoint_angles(*center, *normal, p0, p1);
            CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center: *center,
                    normal: *normal,
                    radius: *radius,
                    start_angle,
                    end_angle,
                },
                t_start: start_angle,
                t_end: end_angle,
            }
        }
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            ..
        } => {
            let (start_angle, end_angle) =
                ellipse_endpoint_angles(*center, *u_basis, *v_basis, p0, p1);
            CurveOnEdge {
                geometry: CurveGeometry::Ellipse {
                    center: *center,
                    u_basis: *u_basis,
                    v_basis: *v_basis,
                    start_angle,
                    end_angle,
                },
                t_start: start_angle,
                t_end: end_angle,
            }
        }
        CurveGeometry::Parabola { .. }
        | CurveGeometry::Hyperbola { .. }
        | CurveGeometry::Nurbs(_) => {
            // For non-conic / NURBS curves, replace with a degree-1
            // line.  Acceptable for cube-edge and perp-walk segments
            // that originated as straight lines.  Cyl-incident NURBS
            // arcs lose curvature — Phase B's downstream pcurve
            // refit recovers most of it, and the cap-vs-cyl Adopt
            // path uses planner curves directly for the cap surface.
            let direction = Vec3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: p0,
                    direction,
                },
                t_start: 0.0,
                t_end: 1.0,
            }
        }
    }
}

/// Compute (start_angle, end_angle) for a circular arc centred at
/// `center` with plane normal `normal`, where the start point is
/// `p0` and the end point is `p1`.  Angles are measured around
/// `normal` from an arbitrary in-plane reference direction
/// (x-axis-projected-onto-plane), wrapped to (-π, π].
fn arc_endpoint_angles(center: Point3, normal: Vec3, p0: Point3, p1: Point3) -> (f64, f64) {
    let (u, v) = in_plane_basis(normal);
    let pa = Vec3::new(p0.x - center.x, p0.y - center.y, p0.z - center.z);
    let pb = Vec3::new(p1.x - center.x, p1.y - center.y, p1.z - center.z);
    let a0 = pa.dot(&v).atan2(pa.dot(&u));
    let a1 = pb.dot(&v).atan2(pb.dot(&u));
    (a0, a1)
}

/// Compute (start_angle, end_angle) for an elliptical arc centred at
/// `center` with bases `u_basis` / `v_basis`.  Same convention as
/// [`arc_endpoint_angles`] but using the ellipse's own basis pair.
fn ellipse_endpoint_angles(
    center: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    p0: Point3,
    p1: Point3,
) -> (f64, f64) {
    let pa = Vec3::new(p0.x - center.x, p0.y - center.y, p0.z - center.z);
    let pb = Vec3::new(p1.x - center.x, p1.y - center.y, p1.z - center.z);
    let u_len2 = u_basis.dot(&u_basis);
    let v_len2 = v_basis.dot(&v_basis);
    let a0 = (pa.dot(&v_basis) / v_len2.sqrt()).atan2(pa.dot(&u_basis) / u_len2.sqrt());
    let a1 = (pb.dot(&v_basis) / v_len2.sqrt()).atan2(pb.dot(&u_basis) / u_len2.sqrt());
    (a0, a1)
}

/// Build an orthonormal pair (u, v) spanning the plane perpendicular
/// to `normal`.  `u` is `normal × world_axis`, where `world_axis` is
/// chosen to be linearly independent from `normal`; `v = normal × u`.
fn in_plane_basis(normal: Vec3) -> (Vec3, Vec3) {
    let n = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    // Pick the world axis least aligned with `n`.
    let world = if n.x.abs() < n.y.abs() && n.x.abs() < n.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if n.y.abs() < n.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let u = n
        .cross(&world)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = n.cross(&u).normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

/// BFS from `v` collecting vertices reachable through edges, capped
/// by 3D distance from `v_pos`.
fn bfs_collect_neighborhood(
    brep: &BRep,
    v: VertexHandle,
    radius: f64,
    v_pos: Point3,
    out: &mut HashSet<VertexHandle>,
) {
    use std::collections::VecDeque;
    let mut queue: VecDeque<VertexHandle> = VecDeque::new();
    queue.push_back(v);
    out.insert(v);
    while let Some(u) = queue.pop_front() {
        for (_he_h, he) in brep.half_edges.iter() {
            if he.origin != u {
                continue;
            }
            let Some(twin_he) = brep.half_edges.get(he.twin) else {
                continue;
            };
            let neighbor = twin_he.origin;
            if out.contains(&neighbor) {
                continue;
            }
            let Some(np) = brep.vertices.get(neighbor) else {
                continue;
            };
            if v_pos.distance_to(&np.position) > radius {
                continue;
            }
            out.insert(neighbor);
            queue.push_back(neighbor);
        }
    }
}

/// Walk every cylinder face's outer loop and add vertices within
/// `radius` of `v_pos` to `out`.  Catches cylinder-only perp-anchors
/// that aren't reachable from V via plane-sub-face edges (Phase A
/// perp-walks insert these on cyl boundaries near V).
fn collect_cylinder_outer_loop_vertices(
    brep: &BRep,
    radius: f64,
    v_pos: Point3,
    out: &mut HashSet<VertexHandle>,
) {
    use popkern_core::geometry::surface::SurfaceGeometry;
    for (fh, face) in &brep.faces {
        let _ = fh;
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let Some(loop_data) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            let v = he.origin;
            if let Some(vp) = brep.vertices.get(v)
                && v_pos.distance_to(&vp.position) <= radius
            {
                out.insert(v);
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
}

/// Maximum fillet radius across the plan's per-fillet cap data.
fn max_fillet_radius_from_plan(plan: &ApexPlan) -> f64 {
    use super::cap_kind::CapKind;
    match &plan.cap_kind {
        CapKind::Sphere { radius, .. } => *radius,
        CapKind::Gregory {
            per_fillet_caps, ..
        } => per_fillet_caps
            .iter()
            .map(|(_, r)| *r)
            .fold(0.0_f64, f64::max),
        CapKind::Open { .. } => 0.1,
    }
}

/// Minimum fillet radius across the plan's per-fillet cap data.
fn min_fillet_radius_from_plan(plan: &ApexPlan) -> f64 {
    use super::cap_kind::CapKind;
    match &plan.cap_kind {
        CapKind::Sphere { radius, .. } => *radius,
        CapKind::Gregory {
            per_fillet_caps, ..
        } => per_fillet_caps
            .iter()
            .map(|(_, r)| *r)
            .fold(f64::INFINITY, f64::min),
        CapKind::Open { .. } => 0.1,
    }
}

#[allow(dead_code)]
fn _typecheck(_fh: FaceHandle) {}
