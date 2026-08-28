//! Target selection helpers for apex vertex coherence passes.

use super::*;

/// Collect substituted-boundary junctions: the unique 3D
/// positions where adjacent boundary segments meet.  For an
/// alternating FilletEnd/FaceEdge boundary with 2N segments,
/// there are 2N junctions (one per segment-to-next-segment
/// junction).  Duplicates within `1e-6 m` are merged.
pub(super) fn collect_boundary_junctions(boundary: &crate::blend::BlendBoundary) -> Vec<Point3> {
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
pub(super) fn nearest_junction(junctions: &[Point3], pos: Point3, max_dist: f64) -> Option<Point3> {
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

/// Reverse a NURBS curve in-place — knot reflection + CP/weight
/// reversal.  Same direction-flip recipe as
/// `popkern_core::topology::split_helpers::reverse_nurbs_curve`.
pub(super) fn reverse_nurbs(
    curve: &popkern_core::geometry::NurbsCurve,
) -> popkern_core::geometry::NurbsCurve {
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
pub(super) fn find_target_for_vertex(
    brep: &BRep,
    v: VertexHandle,
    cap_set: &BTreeSet<FaceHandle>,
) -> Option<Point3> {
    use popkern_core::topology::EdgeHandle;

    let trace = std::env::var_os("POPKERN_APEX_REBUILD_TRACE").is_some();
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
