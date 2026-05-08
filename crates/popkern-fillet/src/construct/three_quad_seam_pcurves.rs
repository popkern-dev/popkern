//! Boundary-owned seam pcurve authoring for three-quad cap faces.

use super::three_quad_face_matching::{
    face_has_cylinder_surface, face_surface_matches, half_edge_endpoints_match_curve,
};
use super::*;
use crate::construct::helpers;

pub(super) fn author_boundary_owned_seam_pcurves(
    brep: &mut BRep,
    sub_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    subface_index: usize,
) {
    if boundary.segments.len() < 2 * (subface_index + 1) {
        return;
    }
    for (segment_offset, segment_index) in [2 * subface_index, 2 * subface_index + 1]
        .into_iter()
        .enumerate()
    {
        let Some(segment) = boundary.segments.get(segment_index) else {
            continue;
        };
        let Some(he) =
            find_subface_half_edge_for_boundary_segment(brep, sub_face, segment, fillet_surfaces)
        else {
            if crate::a2_trace::cap_g0_enabled() {
                eprintln!(
                    "[a2-trace cap-g0] boundary-owned pcurve: subface[{subface_index}] \
                     segment[{segment_index}] no matching half-edge"
                );
            }
            continue;
        };
        let twin = brep.half_edges[he].twin;
        let Some(cap_surface) = brep.faces.get(sub_face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let twin_face = brep.loops[brep.half_edges[twin].loop_].face;
        let Some(twin_surface) = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let cap_pcurve = best_gregory_boundary_pcurve_for_half_edge(brep, he, cap_surface)
            .unwrap_or_else(|| {
                canonical_gregory_boundary_pcurve_for_half_edge(
                    brep,
                    he,
                    cap_surface,
                    segment_offset,
                )
            });
        let other_pcurve =
            project_boundary_segment_to_half_edge_pcurve(brep, twin, twin_surface, segment);
        brep.half_edges[he].pcurve = Some(cap_pcurve);
        if let Some(pcurve) = other_pcurve {
            brep.half_edges[twin].pcurve = Some(pcurve);
        }
        if crate::a2_trace::cap_g0_enabled() {
            eprintln!(
                "[a2-trace cap-g0] boundary-owned pcurve: subface[{subface_index}] \
                 segment[{segment_index}] wrote he={he:?} twin={twin:?} edge={:?}",
                brep.half_edges[he].edge
            );
        }
    }
}

fn find_subface_half_edge_for_boundary_segment(
    brep: &BRep,
    sub_face: FaceHandle,
    segment: &BoundarySegment,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<popkern_core::topology::HalfEdgeHandle> {
    let face_data = brep.faces.get(sub_face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..(brep.half_edges.len() + 1) {
        let he = brep.half_edges.get(cursor)?;
        let twin_face = brep.loops.get(brep.half_edges.get(he.twin)?.loop_)?.face;
        let matches = match segment {
            BoundarySegment::FaceEdge { face, curve } => {
                twin_face == *face || half_edge_endpoints_match_curve(brep, cursor, curve, 1.0e-5)
            }
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => {
                half_edge_endpoints_match_curve(brep, cursor, curve, 1.0e-5)
                    || fillet_surfaces.get(*fillet_index).is_some_and(|fillet| {
                        face_surface_matches(brep, twin_face, &fillet.surface)
                    })
                    || face_has_cylinder_surface(brep, twin_face)
            }
        };
        if matches {
            return Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    None
}

fn canonical_gregory_boundary_pcurve_for_half_edge(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    cap_surface: &SurfaceOnFace,
    segment_offset: usize,
) -> popkern_core::geometry::pcurve::Pcurve {
    use popkern_core::geometry::surface::SurfaceEval;

    let (start_uv, end_uv) = match segment_offset {
        // Boundary 0: FilletEnd[k], u = 0, v = 0 -> 1.
        0 => (Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        // Boundary 1: FaceEdge[k], v = 1, u = 0 -> 1.
        _ => (Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
    };
    let he_origin = brep.vertices[brep.half_edges[he].origin].position;
    let p_start = cap_surface.geometry.point_at(start_uv.x, start_uv.y);
    let p_end = cap_surface.geometry.point_at(end_uv.x, end_uv.y);
    let (start, end) = if he_origin.distance_to(&p_start) <= he_origin.distance_to(&p_end) {
        (start_uv, end_uv)
    } else {
        (end_uv, start_uv)
    };
    popkern_core::geometry::pcurve::Pcurve::linear((start.x, start.y), (end.x, end.y), 0.0, 1.0)
}

fn best_gregory_boundary_pcurve_for_half_edge(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    cap_surface: &SurfaceOnFace,
) -> Option<popkern_core::geometry::pcurve::Pcurve> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let edge = brep.half_edges.get(he)?.edge;
    let coe = brep.edges.get(edge)?.curve.as_ref()?;
    let he_origin = brep.vertices.get(brep.half_edges.get(he)?.origin)?.position;
    let p_t_start = coe.geometry.point_at(coe.t_start);
    let p_t_end = coe.geometry.point_at(coe.t_end);
    let edge_forward = he_origin.distance_to(&p_t_start) <= he_origin.distance_to(&p_t_end);

    type BoundaryUv = ((f64, f64), (f64, f64));
    let candidates = [
        ((0.0, 0.0), (0.0, 1.0)),
        ((1.0, 0.0), (1.0, 1.0)),
        ((0.0, 0.0), (1.0, 0.0)),
        ((0.0, 1.0), (1.0, 1.0)),
    ];
    let mut best: Option<(BoundaryUv, f64)> = None;
    for (start, end) in candidates {
        for reverse in [false, true] {
            let (uv0, uv1) = if reverse { (end, start) } else { (start, end) };
            let mut worst = 0.0_f64;
            for i in 0..16 {
                let s = i as f64 / 15.0;
                let edge_s = if edge_forward { s } else { 1.0 - s };
                let edge_t = coe.t_start + edge_s * (coe.t_end - coe.t_start);
                let p_edge = coe.geometry.point_at(edge_t);
                let u = uv0.0 + s * (uv1.0 - uv0.0);
                let v = uv0.1 + s * (uv1.1 - uv0.1);
                let p_cap = cap_surface.geometry.point_at(u, v);
                worst = worst.max(p_edge.distance_to(&p_cap));
            }
            match best {
                Some((_, best_worst)) if worst >= best_worst => {}
                _ => best = Some(((uv0, uv1), worst)),
            }
        }
    }
    let ((uv0, uv1), worst) = best?;
    if crate::a2_trace::cap_g0_enabled() {
        eprintln!(
            "[a2-trace cap-g0] best Gregory side for he={he:?}: ({:.1},{:.1})->({:.1},{:.1}) worst={worst:.6e}",
            uv0.0, uv0.1, uv1.0, uv1.1
        );
    }
    Some(popkern_core::geometry::pcurve::Pcurve::linear(
        uv0, uv1, 0.0, 1.0,
    ))
}

fn project_boundary_segment_to_half_edge_pcurve(
    brep: &BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    surface: &SurfaceOnFace,
    segment: &BoundarySegment,
) -> Option<popkern_core::geometry::pcurve::Pcurve> {
    use popkern_core::geometry::curve::CurveEval;

    let curve = match segment {
        BoundarySegment::FilletEnd { curve, .. } | BoundarySegment::FaceEdge { curve, .. } => curve,
    };
    let (t0, t1) = curve.parameter_range();
    let p0 = curve.point_at(t0);
    let p1 = curve.point_at(t1);
    let he_origin = brep.vertices[brep.half_edges[he].origin].position;
    let forward = he_origin.distance_to(&p0) <= he_origin.distance_to(&p1);

    let n_samples = 12usize;
    let mut uv_points = Vec::with_capacity(n_samples);
    let mut params = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let curve_s = if forward { s } else { 1.0 - s };
        let p = curve.point_at(t0 + curve_s * (t1 - t0));
        let (u, v) = helpers::project_to_uv(&surface.geometry, &p);
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &params, 5)
        .ok()
        .map(Into::into)
}
