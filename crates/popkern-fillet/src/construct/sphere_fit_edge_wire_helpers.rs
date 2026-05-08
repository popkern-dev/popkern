//! Curve and pcurve helpers for sphere-fit edge wiring.

use super::*;

/// Build an analytic Linear pcurve on `surface` by projecting the
/// curve endpoints to surface UV.  Exact when the surface contains
/// the curve and the curve traces a linear path in the surface's
/// UV space (cylinder seam-perpendicular arcs, plane edges).
/// Approximate but small-error for short arcs on cylinders where
/// the (theta, axial) trajectory is mildly curved.
///
/// Returns `None` when projection fails (point off surface, degenerate
/// UV, etc.) — caller falls through to the deg-3 NURBS refit.
pub(super) fn analytic_linear_pcurve_on_surface(
    curve: &NurbsCurve,
    surface: &SurfaceGeometry,
    t_start: f64,
    t_end: f64,
) -> Option<Pcurve> {
    let p_start = curve.point_at(t_start);
    let p_end = curve.point_at(t_end);
    let uv_start = project_point_to_uv(surface, &p_start)?;
    let uv_end = project_point_to_uv(surface, &p_end)?;
    Some(Pcurve::linear(uv_start, uv_end, t_start, t_end))
}

/// Phase 6c-D: reverse a NurbsCurve by flipping CPs, weights, and
/// reflecting the knot vector. Equivalent to
/// `popkern_core::geometry::edge_triple::analytic_helpers::reverse_nurbs`
/// (kept private there); duplicated locally for cohesion.
pub(super) fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut reversed_cps = curve.control_points.clone();
    reversed_cps.reverse();
    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cps,
        weights: reversed_weights,
    }
}

/// Phase 6c-D: snap a vertex to `new_pos`, syncing every incident
/// half-edge's edge curve CP[0] (primary side) or CP[last] (twin
/// side) to keep loop closure intact. No-op when already within
/// 1e-9 of `new_pos`.
///
/// This is a drop-in equivalent of
/// `var_radius::update_incident_edges_for_vertex_move` — kept
/// local to the wire pass for cohesion. Both functions produce
/// identical results on identical inputs.
pub(super) fn snap_vertex_to(brep: &mut BRep, v: VertexHandle, new_pos: Point3) {
    let Some(vd) = brep.vertices.get(v) else {
        return;
    };
    if vd.position.distance_to(&new_pos) < 1e-9 {
        return;
    }

    let halves: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(h, he)| (he.origin == v).then_some(h))
        .collect();

    for he_h in halves {
        let edge_h = brep.half_edges[he_h].edge;
        let primary = brep.edges.get(edge_h).map(|e| e.half_edge);
        let is_primary = primary == Some(he_h);
        let Some(edge) = brep.edges.get_mut(edge_h) else {
            continue;
        };
        let Some(coe) = edge.curve.as_mut() else {
            continue;
        };
        if let CurveGeometry::Nurbs(nc) = &mut coe.geometry
            && !nc.control_points.is_empty()
        {
            let n = nc.control_points.len();
            if is_primary {
                nc.control_points[0] = new_pos;
            } else {
                nc.control_points[n - 1] = new_pos;
            }
        }
    }

    if let Some(vd) = brep.vertices.get_mut(v) {
        vd.position = new_pos;
    }
}

/// Sample `curve_3d` at uniform parameter values, project each
/// sample to `surface`'s UV, and fit a degree-3 NURBS pcurve
/// through those (u, v, 0) points. Returns `None` if any
/// projection fails or the fit raises an error.
/// Phase 6c-D Stage 4 — analytic `Linear` cap-side pcurve for an
/// edge wired with a sphere-fit-substituted boundary segment.
///
/// Segment index parity tells us which cap boundary the edge is
/// on (per the layout `sphere_fit_boundary::apply_sphere_fit_substitution`
/// emits — `segments[2k] = FilletEnd[k]`,
/// `segments[2k+1] = FaceEdge[k]`):
///
/// * `seg_idx % 2 == 0` (FilletEnd) → cap **row 0** (`u = 0`):
///   pcurve traces `(0, v)` over `v ∈ [v_lo, v_hi]`.
/// * `seg_idx % 2 == 1` (FaceEdge) → cap **col 9** (`v = 1`):
///   pcurve traces `(u, 1)` over `u ∈ [u_lo, u_hi]`.
///
/// `reversed` flips the start/end UVs so the pcurve advances in
/// the same direction as the underlying 3D curve. Parameter
/// range `[t_start, t_end]` matches the edge curve.
///
/// Phase 6c-D Stage 4.1 — `(boundary_lo, boundary_hi)` is the
/// sub-edge's sub-range on the parent segment, normalised to
/// `[0, 1]`. For a full-edge match (single sub-edge spans whole
/// segment), `(boundary_lo, boundary_hi) = (0, 1)`. For multi-
/// sub-edge cases, each sub-edge gets a sub-range computed by
/// projecting its endpoint vertices onto the parent segment.
pub(super) fn cap_pcurve_for_segment(
    seg_idx: usize,
    reversed: bool,
    t_start: f64,
    t_end: f64,
    boundary_lo: f64,
    boundary_hi: f64,
) -> Option<Pcurve> {
    let (start, end) = if seg_idx.is_multiple_of(2) {
        // Row 0: u = 0, v varies over the sub-range.
        ((0.0, boundary_lo), (0.0, boundary_hi))
    } else {
        // Col 9: v = 1, u varies over the sub-range.
        ((boundary_lo, 1.0), (boundary_hi, 1.0))
    };
    let (start, end) = if reversed { (end, start) } else { (start, end) };
    Some(Pcurve::linear(start, end, t_start, t_end))
}

/// Phase 6c-D Stage 4.1 — closest-point projection of `point`
/// onto `curve`. Returns the parameter `t ∈ [t_lo, t_hi]` such
/// that `curve.point_at(t)` is closest to `point`.
///
/// Sampling-based: 100 uniform parameter samples, take the one
/// with minimum distance, refine with 5 bisection iterations.
/// Sufficient accuracy for sub-edge t-range computation in the
/// wire step (we don't need machine precision — micrometre-level
/// is plenty for cap-outer-G0 on a 0.1m cap).
pub(super) fn project_point_on_curve(curve: &NurbsCurve, point: Point3) -> f64 {
    let (t_lo, t_hi) = curve.parameter_range();
    if (t_hi - t_lo).abs() < 1e-15 {
        return t_lo;
    }
    let mut best_t = t_lo;
    let mut best_d = f64::INFINITY;
    const N: usize = 100;
    for i in 0..=N {
        let s = (i as f64) / (N as f64);
        let t = t_lo + s * (t_hi - t_lo);
        let p = curve.point_at(t);
        let d = p.distance_to(&point);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    // Refine via 5 bisection iterations around best_t.
    let mut lo = (best_t - (t_hi - t_lo) / (N as f64)).max(t_lo);
    let mut hi = (best_t + (t_hi - t_lo) / (N as f64)).min(t_hi);
    for _ in 0..5 {
        let mid_lo = lo + (hi - lo) * 0.25;
        let mid_hi = lo + (hi - lo) * 0.75;
        let d_lo = curve.point_at(mid_lo).distance_to(&point);
        let d_hi = curve.point_at(mid_hi).distance_to(&point);
        if d_lo < d_hi {
            hi = (lo + hi) * 0.5;
        } else {
            lo = (lo + hi) * 0.5;
        }
    }
    (lo + hi) * 0.5
}

pub(super) fn refit_pcurve_on_surface(
    curve_3d: &NurbsCurve,
    surface: &SurfaceGeometry,
    t0: f64,
    t1: f64,
) -> Option<NurbsCurve> {
    const N: usize = 12;
    let mut points: Vec<Point3> = Vec::with_capacity(N);
    let mut params: Vec<f64> = Vec::with_capacity(N);
    for i in 0..N {
        let s = i as f64 / (N - 1) as f64;
        let t = t0 + s * (t1 - t0);
        let p3 = curve_3d.point_at(t);
        let (u, v) = project_point_to_uv(surface, &p3)?;
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3usize;
    let n_cps = (degree + 2).min(points.len());
    fitting::fit_curve(&points, &params, degree, n_cps).ok()
}

pub(super) fn collect_loop_half_edges(brep: &BRep, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
    let mut hes = Vec::new();
    let mut cursor = start;
    let mut safety = 64usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    hes
}
