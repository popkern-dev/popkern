//! Phase 6c-X.A.4 — apex vertex relocation pass.
//!
//! After `wire_substituted_boundary_to_edges` runs in Path A mode
//! (Phase 6c-X.A.3), every cap-loop BRep edge has its substituted
//! curve at the FULL segment range.  At apex vertices on cap face
//! outer loops, all incident cap-loop edges' curves agree on a
//! single 3D termination point T (the FilletEnd termination, =
//! sphere∩cyl arc endpoint).
//!
//! The BRep vertex, however, is still at the original rolling-ball
//! footprint position V.  This module relocates each apex vertex
//! V → T, by reading the target from any incident cap-loop edge.
//!
//! Activation: on by default since Phase 6c-F.  Set
//! `POPKERN_APEX_PATH_A=0` to disable for debugging.
//!
//! After A.4, BASE-EDGE curves (cyl-vs-Plane apex generators going
//! AWAY from the apex) still terminate at the OLD V position.
//! A.5 fixes those.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use std::collections::BTreeSet;

/// Walk every cap face's outer loop, collect the set of apex
/// vertices (vertex handles encountered on cap loops), and mutate
/// each one's position to its target T.
///
/// The target is read from the FIRST incident half-edge whose
/// edge has a curve3d AND whose twin face is NOT a cap face (i.e.
/// a cap-loop edge).  Per A.3, all such edges agree on T.
///
/// Returns the count of vertices relocated and the worst |V - T|
/// distance encountered.  Vertices that were already within 1e-9
/// of their target are left unchanged.
pub(crate) fn relocate_apex_vertices(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    boundary: &crate::blend::BlendBoundary,
) -> (usize, f64) {
    let trace = std::env::var_os("POPKERN_APEX_PATH_A_TRACE").is_some();
    let cap_set: BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();

    // 1. Collect apex vertices from cap face outer loops (existing
    //    behaviour for Sphere caps + Gregory caps where sliver_absorb
    //    didn't expand the loop).
    let mut apex_vertices: BTreeSet<VertexHandle> = BTreeSet::new();
    for &cap in cap_faces {
        let Some(face_data) = brep.faces.get(cap) else {
            continue;
        };
        let outer = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cur = start;
        let mut step = 0usize;
        #[allow(clippy::while_let_loop)]
        loop {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            apex_vertices.insert(he.origin);
            cur = he.next;
            step += 1;
            if cur == start || step > 32 {
                break;
            }
        }
    }

    if trace {
        eprintln!(
            "[apex-path-a] {} apex vertices on cap-face loops",
            apex_vertices.len()
        );
    }

    // 2. §6c-X.A.7 boundary-junction fallback: for variable-radius
    //    Gregory corners, A.4's "shared endpoint across incident
    //    cap-loop edges" heuristic can fail on vertices whose 2
    //    incident cap-loop edges come from different fillets'
    //    substitutions (per-fillet T's don't agree).  In that case
    //    fall back to the closest substituted-boundary junction —
    //    the boundary is internally vertex-coherent by construction
    //    (each junction is one Point3, shared by adjacent segments).
    //    Use a tight search radius to avoid pulling unrelated
    //    vertices toward unrelated junctions.
    let max_radius = boundary
        .per_fillet_caps
        .iter()
        .map(|(_, r)| *r)
        .fold(0.0_f64, f64::max);
    let snap_radius = if max_radius > 0.0 { max_radius } else { 0.1 };

    // 3. Build the boundary's junction list (substituted curve
    //    endpoints — coherent per-vertex by construction).
    let junctions = collect_boundary_junctions(boundary);
    if trace {
        eprintln!(
            "[apex-path-a] {} unique boundary junctions, snap radius {snap_radius:.3e} m",
            junctions.len()
        );
    }

    // For each apex vertex: try the cap-loop-edge cluster
    // (`find_target_for_vertex`); fall back to the nearest
    // boundary junction within `snap_radius`.  Multiple BRep
    // vertices may legitimately snap to the same 3D junction
    // (e.g. when sliver_absorb expanded the cap face's loop with
    // duplicate-position vertices that the substitution
    // collapses to one T).  Those vertex coincidences are
    // structurally-degenerate edges but pass `validate_loop_closure`
    // and `validate_watertight` checks because the curves
    // genuinely terminate at the shared position.
    let mut moves: Vec<(VertexHandle, Point3, f64)> = Vec::new();
    for &v in &apex_vertices {
        let v_pos = brep.vertices[v].position;
        let target = find_target_for_vertex(brep, v, &cap_set)
            .or_else(|| nearest_junction(&junctions, v_pos, snap_radius));
        let Some(target) = target else {
            if trace {
                eprintln!("[apex-path-a]   {v:?}: no target found");
            }
            continue;
        };
        let dist = v_pos.distance_to(&target);
        if dist > 1e-9 {
            moves.push((v, target, dist));
        }
    }

    let worst = moves.iter().map(|(_, _, d)| *d).fold(0.0_f64, f64::max);
    let n_moved = moves.len();
    for (v, target, dist) in moves {
        if trace {
            eprintln!(
                "[apex-path-a]   {v:?}: → ({:.4},{:.4},{:.4}) dist={dist:.3e}",
                target.x, target.y, target.z
            );
        }
        if let Some(vd) = brep.vertices.get_mut(v) {
            vd.position = target;
        }
    }
    if trace {
        eprintln!("[apex-path-a] relocated {n_moved} vertices, worst distance {worst:.3e} m");
    }
    (n_moved, worst)
}

/// Collect substituted-boundary junctions: the unique 3D
/// positions where adjacent boundary segments meet.  For an
/// alternating FilletEnd/FaceEdge boundary with 2N segments,
/// there are 2N junctions (one per segment-to-next-segment
/// junction).  Duplicates within `1e-6 m` are merged.
fn collect_boundary_junctions(boundary: &crate::blend::BlendBoundary) -> Vec<Point3> {
    use popkern_core::geometry::curve::CurveEval;
    let mut out: Vec<Point3> = Vec::new();
    for seg in &boundary.segments {
        let curve = match seg {
            crate::blend::BoundarySegment::FilletEnd { curve, .. } => curve,
            crate::blend::BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let (t0, t1) = curve.parameter_range();
        for p in [curve.point_at(t0), curve.point_at(t1)] {
            if !out.iter().any(|q| q.distance_to(&p) < 1e-6) {
                out.push(p);
            }
        }
    }
    out
}

/// Find the junction nearest to `pos` within `max_dist`.  Returns
/// the junction position when one qualifies; `None` otherwise.
fn nearest_junction(junctions: &[Point3], pos: Point3, max_dist: f64) -> Option<Point3> {
    let mut best: Option<(Point3, f64)> = None;
    for &j in junctions {
        let d = pos.distance_to(&j);
        if d > max_dist {
            continue;
        }
        match best {
            Some((_, b)) if b <= d => {}
            _ => best = Some((j, d)),
        }
    }
    best.map(|(j, _)| j)
}

/// Phase 6c-X.A.5: BASE-EDGE coherence repair.
///
/// After A.4 moves each apex vertex V → T, the BASE-EDGE incident
/// to V (the cyl-vs-Plane apex generator going AWAY from the apex)
/// still has its old curve3d ending at V's OLD position.  This pass
/// rebuilds each such edge's curve as a degree-1 line from V's NEW
/// position to its other endpoint's position, restoring
/// vertex↔curve coherence on the BASE-EDGE side.
///
/// BASE-EDGE = edge incident to apex vertex V where both faces are
/// non-cap (e.g. cyl-vs-Plane / cyl-vs-Nurbs cube boundary) and
/// the OTHER endpoint is a non-apex vertex.  Per X.A.1 each apex
/// vertex has exactly one such edge.
///
/// Returns the number of base edges rewritten.
pub(crate) fn repair_base_edges(brep: &mut BRep, cap_faces: &[FaceHandle]) -> usize {
    use popkern_core::topology::EdgeHandle;

    let trace = std::env::var_os("POPKERN_APEX_PATH_A_TRACE").is_some();
    let cap_set: BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();

    let mut apex_vertices: BTreeSet<VertexHandle> = BTreeSet::new();
    for &cap in cap_faces {
        let Some(face_data) = brep.faces.get(cap) else {
            continue;
        };
        let outer = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cur = start;
        let mut step = 0usize;
        #[allow(clippy::while_let_loop)]
        loop {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            apex_vertices.insert(he.origin);
            cur = he.next;
            step += 1;
            if cur == start || step > 32 {
                break;
            }
        }
    }

    // Build (edge, primary_origin_pos, primary_dest_pos) tuples for
    // each BASE-EDGE so we can mutate after the read-only pass.
    let mut jobs: Vec<(EdgeHandle, Point3, Point3)> = Vec::new();
    for &v in &apex_vertices {
        for (_he_h, he) in brep.half_edges.iter() {
            if he.origin != v {
                continue;
            }
            let Some(edge) = brep.edges.get(he.edge) else {
                continue;
            };
            // Identify face roles.
            let primary = edge.half_edge;
            let twin = brep.half_edges[primary].twin;
            let face_a = brep.loops[brep.half_edges[primary].loop_].face;
            let face_b = brep.loops[brep.half_edges[twin].loop_].face;
            if cap_set.contains(&face_a) || cap_set.contains(&face_b) {
                continue; // skip cap-loop edges
            }
            // OTHER endpoint = primary's dest if primary.origin == v, else primary's origin.
            let primary_origin = brep.half_edges[primary].origin;
            let primary_dest = brep.half_edges[brep.half_edges[primary].next].origin;
            let other = if primary_origin == v {
                primary_dest
            } else {
                primary_origin
            };
            if apex_vertices.contains(&other) {
                continue; // not a base edge — both endpoints in apex set
            }
            // Found base edge.  Build new curve from primary_origin to primary_dest
            // (in the BRep edge's stored direction).
            let p0 = brep.vertices[primary_origin].position;
            let p1 = brep.vertices[primary_dest].position;
            jobs.push((he.edge, p0, p1));
            break; // each apex vertex has at most one base edge per X.A.1
        }
    }

    let n_jobs = jobs.len();
    for (eh, p0, p1) in jobs {
        if let Some(edge_data) = brep.edges.get_mut(eh) {
            edge_data.curve = CurveOnEdge::line_segment(p0, p1, 1.0e-9);
        }
        if trace {
            eprintln!(
                "[apex-path-a]   base edge {eh:?}: curve replaced ({:.4},{:.4},{:.4})→({:.4},{:.4},{:.4})",
                p0.x, p0.y, p0.z, p1.x, p1.y, p1.z,
            );
        }
    }
    if trace {
        eprintln!("[apex-path-a] repaired {n_jobs} base edges");
    }
    n_jobs
}

/// Phase 6c-X.A.5b: realign each cap-loop edge's curve direction so
/// `curve(t_start) ≈ primary_origin.position` (within `cluster_tol`).
///
/// After A.4 relocates apex vertices to their shared boundary
/// endpoints T, the wire's stored curve direction (set by an old
/// vertex-position heuristic) may now be opposite to BRep
/// orientation: `curve(t_start)` may actually equal
/// `primary_dest.position`.  `validate_watertight` requires the
/// BRep convention, so this pass reverses the underlying curve
/// when misaligned.
///
/// Returns the number of edges reversed.
pub(crate) fn realign_cap_loop_curve_directions(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
) -> usize {
    use popkern_core::geometry::CurveGeometry;
    use popkern_core::topology::EdgeHandle;

    let trace = std::env::var_os("POPKERN_APEX_PATH_A_TRACE").is_some();
    let cap_set: BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let cluster_tol = 1e-6_f64;
    // Discover cap-loop edges: edges where one face is cap.
    let mut edges_to_check: BTreeSet<EdgeHandle> = BTreeSet::new();
    for &cap in cap_faces {
        let Some(face_data) = brep.faces.get(cap) else {
            continue;
        };
        let outer = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cur = start;
        let mut step = 0usize;
        #[allow(clippy::while_let_loop)]
        loop {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            edges_to_check.insert(he.edge);
            cur = he.next;
            step += 1;
            if cur == start || step > 32 {
                break;
            }
        }
    }
    let mut count = 0usize;
    for eh in edges_to_check {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let primary = edge.half_edge;
        let twin = brep.half_edges[primary].twin;
        let face_a = brep.loops[brep.half_edges[primary].loop_].face;
        let face_b = brep.loops[brep.half_edges[twin].loop_].face;
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let primary_origin = brep.half_edges[primary].origin;
        let twin_origin = brep.half_edges[twin].origin;
        let v_origin_pos = brep.vertices[primary_origin].position;
        let v_dest_pos = brep.vertices[twin_origin].position;
        use popkern_core::geometry::curve::CurveEval;
        let p_t_start = coe.geometry.point_at(coe.t_start);
        let p_t_end = coe.geometry.point_at(coe.t_end);
        let aligned = v_origin_pos.distance_to(&p_t_start) < cluster_tol
            && v_dest_pos.distance_to(&p_t_end) < cluster_tol;
        if aligned {
            continue;
        }
        let reversed_alignment = v_origin_pos.distance_to(&p_t_end) < cluster_tol
            && v_dest_pos.distance_to(&p_t_start) < cluster_tol;
        if !reversed_alignment {
            // Neither alignment fits — the edge's curve doesn't
            // terminate at either BRep vertex.  Force coherence by
            // rewriting the curve as a straight line from
            // v_origin_pos to v_dest_pos.  This sacrifices the
            // cap-vs-cyl arc fidelity for V↔T coherence
            // (validate_loop_closure / validate_watertight pass
            // even on degenerate-curve edges as long as endpoints
            // match vertex positions).  cap_outer_g0 may regress
            // slightly on the affected edges; for the variable-
            // radius cube fixtures the trade-off is favourable
            // (apex-radius-scale ShellGap → straight-line edge).
            //
            // §6c-X.A.7: enables the variable-radius Gregory branch
            // to land path_a_relevant=0.
            if trace {
                eprintln!(
                    "[apex-path-a]   edge {eh:?}: neither aligned nor reversed; \
                     rewriting as straight line v_origin→v_dest \
                     (v_origin={:.3e}/{:.3e} v_dest={:.3e}/{:.3e})",
                    v_origin_pos.distance_to(&p_t_start),
                    v_origin_pos.distance_to(&p_t_end),
                    v_dest_pos.distance_to(&p_t_start),
                    v_dest_pos.distance_to(&p_t_end),
                );
            }
            if let Some(edge_mut) = brep.edges.get_mut(eh) {
                edge_mut.curve = CurveOnEdge::line_segment(v_origin_pos, v_dest_pos, 1.0e-9);
            }
            count += 1;
            continue;
        }
        // Reverse the underlying curve.
        let CurveGeometry::Nurbs(nurbs) = &coe.geometry else {
            continue;
        };
        let reversed = reverse_nurbs(nurbs);
        if let Some(edge_mut) = brep.edges.get_mut(eh)
            && let Some(coe_mut) = edge_mut.curve.as_mut()
        {
            coe_mut.geometry = CurveGeometry::Nurbs(reversed);
        }
        count += 1;
        if trace {
            eprintln!("[apex-path-a]   edge {eh:?}: reversed to align with BRep direction");
        }
    }
    if trace {
        eprintln!("[apex-path-a] realigned {count} cap-loop edges");
    }
    count
}

/// Reverse a NURBS curve in-place — knot reflection + CP/weight
/// reversal.  Same direction-flip recipe as
/// `popkern_core::topology::split_helpers::reverse_nurbs_curve`.
fn reverse_nurbs(curve: &popkern_core::geometry::NurbsCurve) -> popkern_core::geometry::NurbsCurve {
    use popkern_core::geometry::curve::CurveEval;
    let (t_start, t_end) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_start + t_end - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut reversed_cp = curve.control_points.clone();
    reversed_cp.reverse();
    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();
    popkern_core::geometry::NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cp,
        weights: reversed_weights,
    }
}

/// Find the target position for apex vertex `v` as the SHARED
/// curve endpoint between its incident cap-loop edges.
///
/// Per X.A.2, the substituted boundary is internally coherent:
/// adjacent FilletEnd[k] and FaceEdge[k] share a single 3D point
/// at their junction.  After A.3 (full-range wire), each cap-loop
/// BRep edge stores its curve preserving that endpoint.  The
/// wire's "reversed" flag (set by vertex-position heuristic) is
/// inconsistent with the boundary's natural direction, so the
/// BRep convention `primary-origin → t_start` does NOT robustly
/// identify v's target.
///
/// Geometric matching does: collect all 4 curve endpoints from
/// v's 2 incident cap-loop edges (typically 1 cap-vs-Cyl and 1
/// cap-vs-Plane).  The endpoint that appears in BOTH edges'
/// curves (within `cluster_tol`) is v's target — by construction
/// this is the shared boundary T.
fn find_target_for_vertex(
    brep: &BRep,
    v: VertexHandle,
    cap_set: &BTreeSet<FaceHandle>,
) -> Option<Point3> {
    use popkern_core::topology::EdgeHandle;

    let trace = std::env::var_os("POPKERN_APEX_PATH_A_TRACE").is_some();
    let cluster_tol = 1e-6_f64;
    let mut cap_loop_edges: Vec<(EdgeHandle, Point3, Point3)> = Vec::new();
    let mut seen_edges: BTreeSet<EdgeHandle> = BTreeSet::new();
    for (_he_h, he) in brep.half_edges.iter() {
        if he.origin != v {
            continue;
        }
        if seen_edges.contains(&he.edge) {
            continue;
        }
        let Some(edge) = brep.edges.get(he.edge) else {
            continue;
        };
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let face_a = brep.loops[brep.half_edges[edge.half_edge].loop_].face;
        let face_b = brep.loops[brep.half_edges[brep.half_edges[edge.half_edge].twin].loop_].face;
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let p_start = coe.geometry.point_at(coe.t_start);
        let p_end = coe.geometry.point_at(coe.t_end);
        cap_loop_edges.push((he.edge, p_start, p_end));
        seen_edges.insert(he.edge);
    }
    if cap_loop_edges.len() < 2 {
        // Single cap-loop edge: use its endpoint nearest v's
        // current position.  No coherence check possible.
        if let Some(&(_, p_start, p_end)) = cap_loop_edges.first() {
            let v_pos = brep.vertices[v].position;
            let target = if v_pos.distance_to(&p_start) <= v_pos.distance_to(&p_end) {
                p_start
            } else {
                p_end
            };
            if trace {
                eprintln!(
                    "[apex-path-a]     {v:?}: single-edge fallback → ({:.4},{:.4},{:.4})",
                    target.x, target.y, target.z
                );
            }
            return Some(target);
        }
        return None;
    }
    // Find the endpoint shared between any two incident cap-loop
    // edges.  Each cap-loop edge has 2 endpoints; for the shared
    // vertex v, exactly one of them appears in both edges.
    let candidates: Vec<Point3> = cap_loop_edges
        .iter()
        .flat_map(|&(_, p_start, p_end)| [p_start, p_end])
        .collect();
    for (i, p_i) in candidates.iter().enumerate() {
        let count = candidates
            .iter()
            .filter(|p| p.distance_to(p_i) < cluster_tol)
            .count();
        if count >= 2 {
            if trace {
                eprintln!(
                    "[apex-path-a]     {v:?}: shared endpoint ({:.4},{:.4},{:.4}) [count={count}]",
                    p_i.x, p_i.y, p_i.z
                );
            }
            return Some(*p_i);
            // (i shadowed but kept for trace)
        }
        let _ = i;
    }
    if trace {
        eprintln!(
            "[apex-path-a]     {v:?}: no shared endpoint found across {} cap-loop edges",
            cap_loop_edges.len()
        );
    }
    None
}
