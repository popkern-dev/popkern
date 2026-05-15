//! Helper routines for edge/surface compatibility repair.

use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::{
    EdgeAnchors, EdgeHint, EdgeTriple, EdgeTripleOptions, SurfaceSide, compute_edge_geometry,
    compute_edge_geometry_unanchored,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::bounded_parameter_range;
use crate::healing::deviation::edge_surface_deviation;
use crate::topology::brep_edge_surface_compatibility_boundary::infer_surface_boundary_hint;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

pub(super) fn build_edge_triple(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    start: crate::geometry::Point3,
    end: crate::geometry::Point3,
    hint: Option<EdgeHint>,
    brep: &BRep,
) -> Result<EdgeTriple, String> {
    let mut edge_options = EdgeTripleOptions {
        marching_timeout_ms: Some(100),
        ..Default::default()
    };
    edge_options.verify_in_release = true;
    compute_edge_geometry(
        surf_a,
        surf_b,
        EdgeAnchors::Open { start, end },
        hint,
        &brep.tolerance,
        &edge_options,
    )
    .map_err(|err| format!("{err:?}"))
}

pub(super) fn build_reanchored_intersection_triple(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    start: crate::geometry::Point3,
    end: crate::geometry::Point3,
    brep: &BRep,
    samples: usize,
) -> Result<(EdgeTriple, f64), String> {
    let mut edge_options = EdgeTripleOptions {
        marching_timeout_ms: Some(100),
        ..Default::default()
    };
    edge_options.verify_in_release = true;
    let vertices = [start, end];
    let (u1_range, v1_range) = bounded_parameter_range(&vertices, surf_a);
    let (u2_range, v2_range) = bounded_parameter_range(&vertices, surf_b);
    let candidates = compute_edge_geometry_unanchored(
        surf_a,
        u1_range,
        v1_range,
        surf_b,
        u2_range,
        v2_range,
        &brep.tolerance,
        &edge_options,
    )
    .map_err(|err| format!("{err:?}"))?;
    candidates
        .into_iter()
        .filter_map(|triple| reanchor_candidate(triple, start, end, samples))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .ok_or_else(|| "no_unanchored_branch_near_edge_anchors".to_string())
}

fn reanchor_candidate(
    mut triple: EdgeTriple,
    start: crate::geometry::Point3,
    end: crate::geometry::Point3,
    samples: usize,
) -> Option<(EdgeTriple, f64)> {
    let samples = samples.max(17);
    let (range_start, range_end) = (triple.t_start, triple.t_end);
    let start_param = closest_curve_parameter(&triple, start, range_start, range_end, samples)?;
    let end_param = closest_curve_parameter(&triple, end, range_start, range_end, samples)?;
    if (end_param - start_param).abs() <= 1e-12 {
        return None;
    }
    let start_motion = triple.curve_3d.point_at(start_param).distance_to(&start);
    let end_motion = triple.curve_3d.point_at(end_param).distance_to(&end);
    triple.t_start = start_param;
    triple.t_end = end_param;
    Some((triple, start_motion.max(end_motion)))
}

fn closest_curve_parameter(
    triple: &EdgeTriple,
    target: crate::geometry::Point3,
    t_start: f64,
    t_end: f64,
    samples: usize,
) -> Option<f64> {
    let mut best_index = 0usize;
    let mut best_distance = f64::INFINITY;
    for index in 0..samples {
        let f = index as f64 / (samples - 1) as f64;
        let t = t_start + (t_end - t_start) * f;
        let distance = triple.curve_3d.point_at(t).distance_to(&target);
        if distance < best_distance {
            best_distance = distance;
            best_index = index;
        }
    }
    if !best_distance.is_finite() {
        return None;
    }
    let step = (t_end - t_start) / (samples - 1) as f64;
    let mut lo = t_start + step * best_index.saturating_sub(1) as f64;
    let mut hi = t_start + step * (best_index + 1).min(samples - 1) as f64;
    if lo > hi {
        std::mem::swap(&mut lo, &mut hi);
    }
    for _ in 0..32 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        let d1 = triple.curve_3d.point_at(m1).distance_to(&target);
        let d2 = triple.curve_3d.point_at(m2).distance_to(&target);
        if d1 <= d2 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    Some(0.5 * (lo + hi))
}

pub(super) fn infer_boundary_hint(
    brep: &mut BRep,
    he_a: HalfEdgeHandle,
    he_b: HalfEdgeHandle,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tolerance: f64,
) -> Option<EdgeHint> {
    infer_surface_boundary_hint(brep, he_a, surf_a, SurfaceSide::A, tolerance)
        .or_else(|| infer_surface_boundary_hint(brep, he_b, surf_b, SurfaceSide::B, tolerance))
}

pub(super) fn endpoint_gap(
    curve: &dyn CurveEval,
    t_start: f64,
    t_end: f64,
    start: crate::geometry::Point3,
    end: crate::geometry::Point3,
) -> f64 {
    curve
        .point_at(t_start)
        .distance_to(&start)
        .max(curve.point_at(t_end).distance_to(&end))
}

pub(super) fn residual_for_candidate(
    curve: &dyn CurveEval,
    t_range: (f64, f64),
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    samples: usize,
) -> f64 {
    let a = edge_surface_deviation(curve, t_range, surf_a, samples).max_deviation;
    let b = edge_surface_deviation(curve, t_range, surf_b, samples).max_deviation;
    a.max(b)
}

pub(super) fn edge_vertices(brep: &BRep, edge: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let half_edge = brep.edges.get(edge)?.half_edge;
    let start = brep.half_edges.get(half_edge)?.origin;
    let next = brep.half_edges.get(half_edge)?.next;
    let end = brep.half_edges.get(next)?.origin;
    Some((start, end))
}

pub(super) fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    let mut edges = Vec::new();
    for (_, half_edge) in &brep.half_edges {
        let next_origin = brep.half_edges.get(half_edge.next).map(|next| next.origin);
        if half_edge.origin == vertex || next_origin == Some(vertex) {
            edges.push(half_edge.edge);
        }
    }
    edges.sort_by_key(|edge| format!("{edge:?}"));
    edges.dedup();
    edges
}

pub(super) fn move_vertex_with_tolerance(
    brep: &mut BRep,
    vertex: VertexHandle,
    target: crate::geometry::Point3,
    motion: f64,
) {
    let current = brep.vertices[vertex]
        .tolerance
        .unwrap_or(brep.tolerance.linear());
    brep.vertices[vertex].position = target;
    brep.vertices[vertex].tolerance = Some(current);
    if motion > brep.max_imprecision {
        brep.max_imprecision = (motion * 1.01).max(motion + brep.tolerance.numerical_floor());
    }
}

pub(super) fn proposed_vertex_incident_endpoint_gap(
    brep: &BRep,
    target_edge: EdgeHandle,
    vertex: VertexHandle,
    proposed: crate::geometry::Point3,
    tolerance: f64,
) -> f64 {
    let mut max_gap = 0.0_f64;
    for (_, half_edge) in &brep.half_edges {
        let next_origin = brep.half_edges.get(half_edge.next).map(|next| next.origin);
        if half_edge.edge == target_edge
            || (half_edge.origin != vertex && next_origin != Some(vertex))
        {
            continue;
        }
        let Some(curve) = brep
            .edges
            .get(half_edge.edge)
            .and_then(|edge| edge.curve.as_ref())
        else {
            continue;
        };
        let start_gap = curve
            .geometry
            .point_at(curve.t_start)
            .distance_to(&proposed);
        let end_gap = curve.geometry.point_at(curve.t_end).distance_to(&proposed);
        let gap = start_gap.min(end_gap);
        if gap > tolerance {
            max_gap = max_gap.max(gap);
        }
    }
    max_gap
}

pub(super) fn unique_edge_half_edges(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Vec<HalfEdgeHandle> {
    let mut half_edges = Vec::with_capacity(2);
    if brep.half_edges.contains_key(half_edge) {
        half_edges.push(half_edge);
        let twin = brep.half_edges[half_edge].twin;
        if twin != half_edge && brep.half_edges.contains_key(twin) {
            half_edges.push(twin);
        }
    }
    half_edges
}

pub(super) fn face_surface_for_half_edge(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<(FaceHandle, &SurfaceGeometry)> {
    let loop_ = brep.half_edges.get(half_edge)?.loop_;
    let face = brep.loops.get(loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    Some((face, surface))
}

pub(super) fn adjacent_surfaced_pair(
    brep: &BRep,
    edge: EdgeHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle, FaceHandle, FaceHandle)> {
    let edge_data = brep.edges.get(edge)?;
    let half_edges = unique_edge_half_edges(brep, edge_data.half_edge);
    if half_edges.len() != 2 {
        return None;
    }
    let (face_a, _) = face_surface_for_half_edge(brep, half_edges[0])?;
    let (face_b, _) = face_surface_for_half_edge(brep, half_edges[1])?;
    if face_a == face_b {
        return None;
    }
    Some((half_edges[0], half_edges[1], face_a, face_b))
}
