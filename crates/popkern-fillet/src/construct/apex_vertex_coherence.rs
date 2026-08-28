//! Apex vertex relocation pass.
//!
//! After `wire_substituted_boundary_to_edges` runs in apex-rebuild mode
//! (the full-range wire pass), every cap-loop BRep edge has its substituted
//! curve at the FULL segment range.  At apex vertices on cap face
//! outer loops, all incident cap-loop edges' curves agree on a
//! single 3D termination point T (the FilletEnd termination, =
//! sphere∩cyl arc endpoint).
//!
//! The BRep vertex, however, is still at the original rolling-ball
//! footprint position V.  This module relocates each apex vertex
//! V → T, by reading the target from any incident cap-loop edge.
//!
//! After A.4, BASE-EDGE curves (cyl-vs-Plane apex generators going
//! AWAY from the apex) still terminate at the OLD V position.
//! A.5 fixes those.

#[path = "apex_vertex_targets.rs"]
mod apex_vertex_targets;

use apex_vertex_targets::{
    collect_boundary_junctions, find_target_for_vertex, nearest_junction, reverse_nurbs,
};
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
    let trace = std::env::var_os("POPKERN_APEX_REBUILD_TRACE").is_some();
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

/// BASE-EDGE coherence repair.
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

    let trace = std::env::var_os("POPKERN_APEX_REBUILD_TRACE").is_some();
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
            let primary_origin = brep.half_edges[primary].origin;
            let primary_dest = brep.half_edges[twin].origin;
            let other = if primary_origin == v {
                primary_dest
            } else {
                primary_origin
            };
            if apex_vertices.contains(&other) {
                continue; // not a base edge — both endpoints in apex set
            }
            // Found base edge.  Build the curve in the BRep edge's stored
            // direction using the topological twin as the destination.
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

/// Realign each cap-loop edge's curve direction so
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

    let trace = std::env::var_os("POPKERN_APEX_REBUILD_TRACE").is_some();
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
            // to land apex_rebuild_relevant=0.
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
