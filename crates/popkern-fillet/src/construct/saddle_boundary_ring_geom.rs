//! Geometric helpers for the saddle boundary-ring carve (PK-234):
//! junction location/splitting, ring-curve authoring, bridge support
//! selection, and the on-surface postcondition. Split from
//! `saddle_boundary_ring.rs` for the structural file budget.

use super::saddle_boundary_ring::JUNCTION_ON_EDGE_BAND;
use super::*;

/// Find or create the junction vertex where `point` lands on `face`'s
/// outer-loop boundary: reuse a loop vertex within the junction band,
/// else split the containing edge at the point.
pub(super) fn junction_vertex(
    brep: &mut BRep,
    face: FaceHandle,
    point: Point3,
    tolerance: &Tolerance,
) -> Result<Option<VertexHandle>, String> {
    // Reuse: nearest existing loop vertex.
    let mut nearest_vertex: Option<(f64, VertexHandle)> = None;
    for he in face_loop_half_edges(brep, face)? {
        let v = brep.half_edges[he].origin;
        if let Some(data) = brep.vertices.get(v) {
            let d = data.position.distance_to(&point);
            if nearest_vertex.map(|(bd, _)| d < bd).unwrap_or(true) {
                nearest_vertex = Some((d, v));
            }
        }
    }
    if let Some((d, v)) = nearest_vertex
        && d <= JUNCTION_ON_EDGE_BAND
    {
        return Ok(Some(v));
    }
    // Split: the loop edge whose curve passes nearest the point.
    // Coarse sample, then windowed parametric refinement — the
    // coarse minimum alone overstates the distance by up to half a
    // sample step (measured 1.56e-3 on a 0.3-long rail whose true
    // point-to-line distance is ~0).
    let mut best: Option<(f64, f64, EdgeHandle)> = None;
    for he in face_loop_half_edges(brep, face)? {
        let edge = brep.half_edges[he].edge;
        let Some(curve) = brep.edges[edge].curve.as_ref() else {
            continue;
        };
        let (e0, e1) = (curve.t_start, curve.t_end);
        if !(e0.is_finite() && e1.is_finite()) {
            continue;
        }
        let eval: &dyn CurveEval = &curve.geometry;
        let mut d_min = f64::INFINITY;
        let mut t_min = e0;
        for i in 0..=64 {
            let t = e0 + (e1 - e0) * i as f64 / 64.0;
            let d = eval.point_at(t).distance_to(&point);
            if d < d_min {
                d_min = d;
                t_min = t;
            }
        }
        let mut window = (e1 - e0).abs() / 64.0;
        for _ in 0..24 {
            let (lo, hi) = (
                (t_min - window).max(e0.min(e1)),
                (t_min + window).min(e0.max(e1)),
            );
            for i in 0..=16 {
                let t = lo + (hi - lo) * i as f64 / 16.0;
                let d = eval.point_at(t).distance_to(&point);
                if d < d_min {
                    d_min = d;
                    t_min = t;
                }
            }
            window *= 0.4;
        }
        if best.map(|(bd, _, _)| d_min < bd).unwrap_or(true) {
            best = Some((d_min, t_min, edge));
        }
    }
    let (d, t_foot, edge) = best.ok_or_else(|| "no curved boundary edge to split".to_string())?;
    if d > JUNCTION_ON_EDGE_BAND {
        if std::env::var_os("PK_CARVE_TRACE").is_some() {
            for he in face_loop_half_edges(brep, face)? {
                let hd = &brep.half_edges[he];
                let p = brep.vertices[hd.origin].position;
                let other = brep.loops[brep.half_edges[hd.twin].loop_].face;
                let mid = brep.edges[hd.edge].curve.as_ref().map(|c| {
                    let eval: &dyn CurveEval = &c.geometry;
                    eval.point_at(0.5 * (c.t_start + c.t_end))
                });
                eprintln!(
                    "    [saddle-ring] {face:?} loop v ({:.4},{:.4},{:.4}) e {:?} mid {:?} | other {other:?}",
                    p.x,
                    p.y,
                    p.z,
                    hd.edge,
                    mid.map(|m| (
                        (m.x * 1e4).round() / 1e4,
                        (m.y * 1e4).round() / 1e4,
                        (m.z * 1e4).round() / 1e4
                    ))
                );
            }
        }
        // The mitre-complex class: no rail edge under the junction.
        // The caller decides whether the rail-spur recut applies.
        if std::env::var_os("PK_CARVE_TRACE").is_some() {
            eprintln!(
                "    [saddle-ring] junction ({:.4},{:.4},{:.4}) off {face:?} boundary by {d:.3e} — recut candidate",
                point.x, point.y, point.z
            );
        }
        return Ok(None);
    }
    // Split at the FOOT on the edge curve (the junction target may
    // sit a hair off the carved edge; the vertex must be on it).
    let foot = {
        let curve = brep.edges[edge]
            .curve
            .as_ref()
            .ok_or_else(|| "split edge lost its curve".to_string())?;
        let eval: &dyn CurveEval = &curve.geometry;
        eval.point_at(t_foot)
    };
    brep.split_edge_at_point(edge, foot, tolerance)
        .map(|r| Some(r.vertex))
        .map_err(|e| format!("split at junction: {e}"))
}

/// Write a boundary segment curve onto a freshly-mef'd ring edge,
/// oriented `v_from → v_to` (the mef argument order) with endpoint
/// control points snapped to the topology vertices.
pub(super) fn write_ring_curve(
    brep: &mut BRep,
    edge: EdgeHandle,
    curve: &NurbsCurve,
    v_from: VertexHandle,
    v_to: VertexHandle,
) {
    let (t0, t1) = curve.parameter_range();
    let p_from = brep.vertices[v_from].position;
    let mut oriented = curve.clone();
    if oriented.point_at(t1).distance_to(&p_from) < oriented.point_at(t0).distance_to(&p_from) {
        oriented = reverse_nurbs(&oriented);
    }
    if let (Some(first), Some(_)) = (
        oriented.control_points.first().copied(),
        oriented.control_points.last().copied(),
    ) {
        let _ = first;
        let len = oriented.control_points.len();
        oriented.control_points[0] = p_from;
        oriented.control_points[len - 1] = brep.vertices[v_to].position;
    }
    let (t0, t1) = oriented.parameter_range();
    if let Some(edge_data) = brep.edges.get_mut(edge) {
        edge_data.curve = Some(popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::CurveGeometry::Nurbs(oriented),
            t_start: t0,
            t_end: t1,
        });
    }
}

/// Reverse a clamped NURBS curve (control points, weights, mirrored
/// knots) so `point_at(t0)` and `point_at(t1)` swap.
fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (lo, hi) = curve.parameter_range();
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    let mut knots: Vec<f64> = curve.knots.iter().map(|k| lo + hi - k).collect();
    knots.reverse();
    NurbsCurve {
        degree: curve.degree,
        control_points,
        weights,
        knots,
    }
}

/// Worst sampled distance from any ring edge curve to its OUTER
/// adjacent face surface (9 samples per edge, cold projection — the
/// outer faces are analytic strips and support planes). The blend
/// side is excluded: its surface is the merge leftover until the
/// downstream cap assignment replaces it.
pub(super) fn worst_ring_on_surface_deviation(
    brep: &BRep,
    ring_edges: &std::collections::BTreeSet<EdgeHandle>,
    blend_face: FaceHandle,
) -> Result<f64, String> {
    let mut worst = 0.0_f64;
    for &eh in ring_edges {
        let edge = brep
            .edges
            .get(eh)
            .ok_or_else(|| format!("ring edge {eh:?} vanished during merge"))?;
        let curve = edge
            .curve
            .as_ref()
            .ok_or_else(|| format!("ring edge {eh:?} lost its curve during merge"))?;
        let he = edge.half_edge;
        let he_data = brep
            .half_edges
            .get(he)
            .ok_or_else(|| format!("ring edge {eh:?} has no half-edge"))?;
        for h in [he, he_data.twin] {
            let Some(hd) = brep.half_edges.get(h) else {
                continue;
            };
            let face = brep.loops[hd.loop_].face;
            if face == blend_face {
                continue;
            }
            let Some(surface) = brep.faces[face].surface.as_ref() else {
                continue;
            };
            let eval: &dyn CurveEval = &curve.geometry;
            for i in 0..=8 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 8.0;
                let p = eval.point_at(t);
                let d = popkern_core::geometry::uv_projection::project_point_to_uv(
                    &surface.geometry,
                    &p,
                )
                .map(|(u, v)| {
                    (&surface.geometry as &dyn SurfaceEval)
                        .point_at(u, v)
                        .distance_to(&p)
                })
                .unwrap_or(f64::INFINITY);
                worst = worst.max(d);
            }
        }
    }
    Ok(worst)
}

/// The support face a bridge crosses: carries BOTH junction vertices
/// on its loop, is not a strip face, and its surface contains the
/// bridge curve's midpoint. The boundary's own face handle cannot
/// make this call — coplanar splits of one support plane are
/// indistinguishable by surface, and the recut's chord kefs move
/// junction ownership.
pub(super) fn select_bridge_support(
    brep: &BRep,
    v_a: VertexHandle,
    v_b: VertexHandle,
    strip_faces: &[FaceHandle],
    curve: &NurbsCurve,
) -> Result<FaceHandle, String> {
    let (t0, t1) = curve.parameter_range();
    let mid = curve.point_at(0.5 * (t0 + t1));
    let mut best: Option<(f64, FaceHandle)> = None;
    for (fh, face_data) in brep.faces.iter() {
        if strip_faces.contains(&fh) {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        let Ok(hes) = face_loop_half_edges(brep, fh) else {
            continue;
        };
        let mut has_a = false;
        let mut has_b = false;
        for he in hes {
            let v = brep.half_edges[he].origin;
            has_a |= v == v_a;
            has_b |= v == v_b;
        }
        if !(has_a && has_b) {
            continue;
        }
        let d = popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &mid)
            .map(|(u, v)| {
                (&surface.geometry as &dyn SurfaceEval)
                    .point_at(u, v)
                    .distance_to(&mid)
            })
            .unwrap_or(f64::INFINITY);
        if best.map(|(bd, _)| d < bd).unwrap_or(true) {
            best = Some((d, fh));
        }
    }
    match best {
        Some((d, fh)) if d <= 1.0e-6 => Ok(fh),
        Some((d, fh)) => Err(format!(
            "best common face {fh:?} is {d:.3e} off the bridge midpoint"
        )),
        None => Err("no common non-strip face carries both junction vertices".to_string()),
    }
}
