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
//! Activation: `POPKERN_APEX_PATH_A=1`.  Default off until the
//! Phase 6c-F flag flip.
//!
//! After A.4, BASE-EDGE curves (cyl-vs-Plane apex generators going
//! AWAY from the apex) still terminate at the OLD V position.
//! A.5 fixes those.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
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
pub(crate) fn relocate_apex_vertices(brep: &mut BRep, cap_faces: &[FaceHandle]) -> (usize, f64) {
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
    if trace {
        eprintln!(
            "[apex-path-a] {} apex vertices to relocate",
            apex_vertices.len()
        );
    }

    let mut moves: Vec<(VertexHandle, Point3, f64)> = Vec::new();
    for &v in &apex_vertices {
        let v_pos = brep.vertices[v].position;
        let target = find_target_for_vertex(brep, v, &cap_set);
        let Some(target) = target else {
            if trace {
                eprintln!("[apex-path-a]   {v:?}: no cap-loop curve found");
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

/// Find the target position for an apex vertex by inspecting one
/// incident cap-loop edge's curve.  Returns the curve's endpoint
/// nearest the BRep vertex's current position (which, after A.3,
/// is the boundary's true T at this vertex).
fn find_target_for_vertex(
    brep: &BRep,
    v: VertexHandle,
    cap_set: &BTreeSet<FaceHandle>,
) -> Option<Point3> {
    let trace = std::env::var_os("POPKERN_APEX_PATH_A_TRACE").is_some();
    let v_pos = brep.vertices[v].position;
    let mut n_he_with_origin_v = 0;
    let mut n_with_curve = 0;
    let n_cap_loop = 0;
    for (he_h, he) in brep.half_edges.iter() {
        if he.origin != v {
            continue;
        }
        n_he_with_origin_v += 1;
        let Some(edge) = brep.edges.get(he.edge) else {
            continue;
        };
        let Some(coe) = edge.curve.as_ref() else {
            if trace {
                eprintln!(
                    "[apex-path-a]     {v:?} HE {he_h:?} edge={:?}: no curve",
                    he.edge
                );
            }
            continue;
        };
        n_with_curve += 1;
        // Cap-loop edge: exactly one side is a cap face.
        let face_a = brep.loops[brep.half_edges[edge.half_edge].loop_].face;
        let face_b = brep.loops[brep.half_edges[brep.half_edges[edge.half_edge].twin].loop_].face;
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue; // both cap (inner spoke) or both non-cap (base edge)
        }
        let _ = n_cap_loop;
        // Pick the endpoint nearest v_pos.  Per A.3 + boundary
        // coherence, every incident cap-loop edge agrees on this
        // endpoint to numerical precision.
        let p_start = coe.geometry.point_at(coe.t_start);
        let p_end = coe.geometry.point_at(coe.t_end);
        let target = if v_pos.distance_to(&p_start) <= v_pos.distance_to(&p_end) {
            p_start
        } else {
            p_end
        };
        return Some(target);
    }
    if trace {
        eprintln!(
            "[apex-path-a]     {v:?}: HEs with origin=V: {n_he_with_origin_v}, with curve: {n_with_curve}, cap-loop: {n_cap_loop}"
        );
    }
    None
}
