//! Geometry helpers for loop walking.
//!
//! Angular comparison, tangent computation, and surface normal lookup
//! used by the walk module to select the correct next half-edge at
//! intersection nodes.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{CurveEval, Vec3};
use popkern_core::topology::BRep;

use crate::merge::graph::{FaceRef, IntersectionGraph, Operand};

/// Pick the next edge from candidates using per-face CW angle.
///
/// Projects edge tangent vectors onto the surface tangent plane at the
/// node, then picks the candidate with the smallest clockwise angle from
/// the incoming direction.  Using the surface normal (instead of the
/// cross product of edge directions) is essential for curved surfaces
/// where edges are not coplanar.
pub(super) fn pick_by_angle(
    graph: &IntersectionGraph,
    he_idx: usize,
    candidates: &[usize],
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<usize> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let node_pos = graph.nodes[end_node].position;

    let origin_node = graph.half_edges[he_idx].origin;
    // For self-loop edges (origin == end), use the curve tangent at the end
    // instead of the zero-length position difference.
    let incoming_dir = if origin_node == end_node {
        incoming_tangent_for_self_loop(graph, he_idx)
    } else {
        graph.nodes[origin_node].position - node_pos
    };

    // Try to get the surface normal from the face's geometry.
    // For curved surfaces this is critical for correct angle computation.
    let normal =
        surface_normal_at_node(graph, end_node, face, brep_a, brep_b).unwrap_or_else(|| {
            // Fallback: compute from incoming direction and candidates.
            // For self-loop candidates, use curve tangent instead of position.
            let mut n = Vec3::new(0.0, 0.0, 1.0);
            for &c in candidates {
                let twin_c = graph.half_edges[c].twin;
                let is_self_loop_c = graph.half_edges[c].origin == graph.half_edges[twin_c].origin;
                let out = if is_self_loop_c {
                    outgoing_tangent_for_self_loop(graph, c)
                } else {
                    graph.nodes[graph.half_edges[twin_c].origin].position - node_pos
                };
                if let Some(computed) = incoming_dir.cross(&out).normalized() {
                    n = computed;
                    break;
                }
            }
            n
        });

    let in_proj = incoming_dir - normal * incoming_dir.dot(&normal);
    let ref_dir = match in_proj.normalized() {
        Some(d) => d,
        None => return Some(candidates[0]),
    };

    let mut best = None;
    let mut best_angle = f64::MAX;

    for &c in candidates {
        let twin_c = graph.half_edges[c].twin;
        let end_pos = graph.nodes[graph.half_edges[twin_c].origin].position;
        // For self-loop candidates (origin == end), use the curve tangent
        // at the start parameter as the outgoing direction.
        let out_dir = if graph.half_edges[c].origin == graph.half_edges[twin_c].origin {
            outgoing_tangent_for_self_loop(graph, c)
        } else {
            end_pos - node_pos
        };
        let out_proj = out_dir - normal * out_dir.dot(&normal);

        let cos_val = out_proj.dot(&ref_dir);
        let sin_val = ref_dir.cross(&out_proj).dot(&normal);
        let mut angle = sin_val.atan2(cos_val);

        if angle <= 1e-10 {
            angle += std::f64::consts::TAU;
        }

        if angle < best_angle {
            best_angle = angle;
            best = Some(c);
        }
    }

    best
}

/// Compute the incoming tangent direction for a self-loop edge.
///
/// For a self-loop (origin == end node), the position-based direction is zero.
/// Instead, evaluate the curve geometry at its end parameter to get the tangent.
/// The incoming direction is negated (pointing back towards where we came from).
pub(super) fn incoming_tangent_for_self_loop(graph: &IntersectionGraph, he_idx: usize) -> Vec3 {
    // Try the forward half-edge's curve first, then the twin's.
    if let Some(ref curve) = graph.half_edges[he_idx].curve {
        let tangent = curve.geometry.tangent_at(curve.t_end);
        // Negate: "incoming" direction at the end node points backwards.
        return tangent * -1.0;
    }
    let twin_idx = graph.half_edges[he_idx].twin;
    if let Some(ref curve) = graph.half_edges[twin_idx].curve {
        // Twin runs in the opposite direction; its tangent at t_start
        // corresponds to the same endpoint. Negate for incoming.
        let tangent = curve.geometry.tangent_at(curve.t_start);
        return tangent;
    }
    // No curve geometry — fall back to zero (pick_by_angle will use candidates[0]).
    Vec3::new(0.0, 0.0, 0.0)
}

/// Compute the outgoing tangent direction for a self-loop candidate edge.
///
/// Similar to [`incoming_tangent_for_self_loop`] but returns the tangent at
/// the start of the curve (the direction the edge departs from the node).
pub(super) fn outgoing_tangent_for_self_loop(graph: &IntersectionGraph, he_idx: usize) -> Vec3 {
    if let Some(ref curve) = graph.half_edges[he_idx].curve {
        return curve.geometry.tangent_at(curve.t_start);
    }
    let twin_idx = graph.half_edges[he_idx].twin;
    if let Some(ref curve) = graph.half_edges[twin_idx].curve {
        // Twin's end tangent (negated) gives the forward direction.
        return curve.geometry.tangent_at(curve.t_end) * -1.0;
    }
    Vec3::new(0.0, 0.0, 0.0)
}

/// Look up the surface normal at a graph node's position on the given face.
///
/// Uses the node's stored UV coordinates (if available) to evaluate the
/// surface normal from the face's geometry.  Returns `None` if the face
/// has no surface geometry or the node has no UV coordinates for this face.
pub(super) fn surface_normal_at_node(
    graph: &IntersectionGraph,
    node_idx: usize,
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<Vec3> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = match face.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let face_data = brep.faces.get(face.face)?;
    let surface = face_data.surface.as_ref()?;
    let geom = &surface.geometry;

    // Try UV coords first.
    if let Some(&(u, v)) = graph.nodes[node_idx].uv_coords.get(face) {
        let n = geom.normal_at(u, v);
        return n.normalized();
    }

    // Fallback: compute surface normal directly from the 3D position for
    // curved analytic surfaces where the cross-product fallback may fail
    // (e.g. self-loop edges produce zero incoming direction).
    // Planes are excluded: the cross-product fallback works reliably for
    // planar faces and changing the normal source could alter edge selection.
    let pos = graph.nodes[node_idx].position;
    match geom {
        SurfaceGeometry::Plane { .. } => None,
        SurfaceGeometry::Cylinder { origin, axis, .. } => {
            let to_pt = pos - *origin;
            let along_axis = *axis * to_pt.dot(axis);
            let radial = to_pt - along_axis;
            radial.normalized()
        }
        SurfaceGeometry::Sphere { center, .. } => {
            let radial = pos - *center;
            radial.normalized()
        }
        _ => None,
    }
}

/// Compute a signed area proxy for a loop using Newell's method.
///
/// Positive = CCW (outer), negative = CW (inner/hole).
pub(super) fn signed_area(graph: &IntersectionGraph, loop_hes: &[usize]) -> f64 {
    // Sum cross products of consecutive vertex positions projected
    // onto an approximate normal.
    let positions: Vec<_> = loop_hes
        .iter()
        .map(|&he_idx| graph.nodes[graph.half_edges[he_idx].origin].position)
        .collect();

    if positions.len() < 3 {
        // Self-loop edge (closed curve like a circle): estimate area from
        // the curve geometry by sampling points along the curve.
        if loop_hes.len() == 1 {
            let he = &graph.half_edges[loop_hes[0]];
            if let Some(ref curve) = he.curve {
                use popkern_core::geometry::curve::CurveEval;
                let n_samples = 16;
                let (t0, t1) = (curve.t_start, curve.t_end);
                let mut pts = Vec::with_capacity(n_samples);
                for i in 0..n_samples {
                    let t = t0 + (t1 - t0) * (i as f64) / (n_samples as f64);
                    pts.push(curve.geometry.point_at(t));
                }
                // Newell's method on sampled points
                let mut nx = 0.0;
                let mut ny = 0.0;
                let mut nz = 0.0;
                for i in 0..pts.len() {
                    let cur = &pts[i];
                    let next = &pts[(i + 1) % pts.len()];
                    nx += (cur.y - next.y) * (cur.z + next.z);
                    ny += (cur.z - next.z) * (cur.x + next.x);
                    nz += (cur.x - next.x) * (cur.y + next.y);
                }
                let normal = Vec3::new(nx, ny, nz);
                return normal.length() * 0.5 * normal.z.signum();
            }
        }
        return 0.0;
    }

    // Newell's method for 3D polygon normal (proportional to area).
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;

    let n = positions.len();
    for i in 0..n {
        let cur = &positions[i];
        let next = &positions[(i + 1) % n];
        nx += (cur.y - next.y) * (cur.z + next.z);
        ny += (cur.z - next.z) * (cur.x + next.x);
        nz += (cur.x - next.x) * (cur.y + next.y);
    }

    let normal = Vec3::new(nx, ny, nz);
    // The length of the normal vector is proportional to twice the area.
    // The sign indicates orientation.
    normal.length() * normal.z.signum()
}
