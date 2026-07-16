//! Geometry and support-side helpers for exact-boundary commit.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3};
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::super::sector::{ExactBoundaryCapSectorInput, ExactBoundaryOuterSide};
use super::super::solve::{ExactBoundarySectorSurface, ExactBoundarySolvedSurfaces};
use super::commit_surface_equivalence::equivalent_surface;

#[path = "commit_trim_curve.rs"]
mod commit_trim_curve;
pub(super) use commit_trim_curve::{canonical_shared_trim_curve, pcurve_trace_curve_deviation};

pub(super) struct SharedTrimCurveSelection<'a> {
    pub(super) edge: EdgeHandle,
    pub(super) cap_half_edge: HalfEdgeHandle,
    pub(super) cap_pcurve: &'a Pcurve,
    pub(super) cap_surface: &'a SurfaceGeometry,
    pub(super) support_half_edge: HalfEdgeHandle,
    pub(super) support_pcurve: &'a Pcurve,
    pub(super) support_surface: &'a SurfaceGeometry,
}

pub(super) fn resolved_support_side_for_commit(
    brep: &BRep,
    side: ExactBoundaryOuterSide,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    prepared_support_face: FaceHandle,
    prepared_support_surface: &SurfaceGeometry,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    let support_face = face_for_half_edge(brep, support_half_edge)?;
    if support_face == prepared_support_face
        && !support_face_needs_semantic_owner(brep, support_face, prepared_support_surface)
    {
        return Some((support_half_edge, support_face));
    }
    if let Some(prepared_side) =
        coincident_prepared_support_side(brep, cap_half_edge, prepared_support_face)
    {
        return Some((prepared_side, prepared_support_face));
    }
    if support_face != cap_face
        && support_face_needs_semantic_owner(brep, support_face, prepared_support_surface)
        && let Some(surface_side) =
            coincident_surface_support_side(brep, cap_face, cap_half_edge, prepared_support_surface)
    {
        return Some(surface_side);
    }
    if side == ExactBoundaryOuterSide::FilletEnd
        && !half_edge_has_cylinder_surface(brep, support_half_edge)
        && let Some(cylinder_side) =
            same_oriented_cylinder_support_side(brep, cap_face, cap_half_edge)
    {
        return Some(cylinder_side);
    }
    (support_face != cap_face).then_some((support_half_edge, support_face))
}

fn support_face_needs_semantic_owner(
    brep: &BRep,
    support_face: FaceHandle,
    prepared_support_surface: &SurfaceGeometry,
) -> bool {
    brep.faces
        .get(support_face)
        .and_then(|face| face.surface.as_ref())
        .is_none_or(|surface| {
            !equivalent_surface(
                &surface.geometry,
                prepared_support_surface,
                brep.tolerance.linear(),
            )
        })
}

fn half_edge_has_cylinder_surface(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    face_for_half_edge(brep, half_edge)
        .and_then(|face| brep.faces.get(face))
        .and_then(|face| face.surface.as_ref())
        .is_some_and(|surface| matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }))
}

fn coincident_surface_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
    prepared_support_surface: &SurfaceGeometry,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    let cap_edge = brep.half_edges.get(cap_half_edge)?.edge;
    let edge_curve = brep.edges.get(cap_edge)?.curve.as_ref();
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !surface_matches_prepared_support(surface, prepared_support_surface, tolerance) {
            continue;
        }
        if let Some(curve) = edge_curve
            && edge_surface_deviation(
                &curve.geometry,
                (curve.t_start, curve.t_end),
                &surface.geometry,
                17,
            )
            .max_deviation
                > tolerance
        {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if coincident_side(brep, cap_half_edge, side) {
                return Some((side, face));
            }
        }
    }
    None
}

fn coincident_prepared_support_side(
    brep: &BRep,
    cap_half_edge: HalfEdgeHandle,
    prepared_support_face: FaceHandle,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(prepared_support_face)?;
    let sides = collect_loop_half_edges(brep, face_data.outer_loop)?;
    sides
        .into_iter()
        .find(|side| coincident_side(brep, cap_half_edge, *side))
}

fn same_oriented_cylinder_support_side(
    brep: &BRep,
    cap_face: FaceHandle,
    cap_half_edge: HalfEdgeHandle,
) -> Option<(HalfEdgeHandle, FaceHandle)> {
    for (face, face_data) in brep.faces.iter() {
        if face == cap_face {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let Some(sides) = collect_loop_half_edges(brep, face_data.outer_loop) else {
            continue;
        };
        for side in sides {
            if same_oriented_side(brep, cap_half_edge, side) {
                return Some((side, face));
            }
        }
    }
    None
}

fn surface_matches_prepared_support(
    candidate: &SurfaceOnFace,
    prepared: &SurfaceGeometry,
    tolerance: f64,
) -> bool {
    equivalent_surface(&candidate.geometry, prepared, tolerance)
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}

fn same_oriented_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    if cap_start != support_start {
        return false;
    }
    half_edge_destination_vertex(brep, cap_side)
        .zip(half_edge_destination_vertex(brep, support_side))
        .is_some_and(|(cap_end, support_end)| cap_end == support_end)
}

fn coincident_side(brep: &BRep, cap_side: HalfEdgeHandle, support_side: HalfEdgeHandle) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(support_start) = brep.half_edges.get(support_side).map(|data| data.origin) else {
        return false;
    };
    let Some(cap_end) = half_edge_destination_vertex(brep, cap_side) else {
        return false;
    };
    let Some(support_end) = half_edge_destination_vertex(brep, support_side) else {
        return false;
    };
    (cap_start == support_start && cap_end == support_end)
        || (cap_start == support_end && cap_end == support_start)
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    let loop_ = brep.half_edges.get(half_edge)?.loop_;
    brep.loops.get(loop_).map(|loop_| loop_.face)
}

fn half_edge_destination_vertex(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
) -> Option<popkern_core::topology::VertexHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    brep.half_edges.get(twin).map(|twin| twin.origin)
}

pub(super) fn project_curve_to_half_edge_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    surface: &SurfaceGeometry,
) -> Option<Pcurve> {
    let forward = half_edge_curve_forward(brep, half_edge, curve)?;
    let samples = 33usize;
    // Periodic surfaces (fillet cylinders) need branch unwrapping: a seam
    // crossing the parameter seam otherwise interpolates across the 2-pi
    // jump and the fitted pcurve leaves the seam mid-curve.
    let (period_u, period_v) = surface.periodicity();
    let mut previous_uv = None;
    let mut uv_points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let curve_fraction = if forward { fraction } else { 1.0 - fraction };
        let curve_t = curve.t_start + (curve.t_end - curve.t_start) * curve_fraction;
        let point = curve.geometry.point_at(curve_t);
        let (mut u, mut v) = surface.closest_point(&point)?;
        crate::apex_construction::exact_boundary_cap::ribbon_projection::unwrap_periodic_uv(
            &mut u,
            &mut v,
            previous_uv,
            period_u,
            period_v,
        );
        previous_uv = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
        parameters.push(fraction);
    }
    let degree = 5.min(uv_points.len().saturating_sub(1));
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &parameters, degree)
        .ok()
        .map(Pcurve::from_nurbs)
}

fn half_edge_curve_forward(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
) -> Option<bool> {
    let origin = brep
        .vertices
        .get(brep.half_edges.get(half_edge)?.origin)?
        .position;
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    Some(origin.distance_to(&curve_start) <= origin.distance_to(&curve_end))
}

pub(super) fn reversed_pcurve(pcurve: &Pcurve) -> Pcurve {
    let (t_start, t_end) = pcurve.parameter_range();
    let samples = 17usize;
    let mut uv_points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let t = t_end + (t_start - t_end) * fraction;
        uv_points.push(pcurve.point_at(t));
        parameters.push(fraction);
    }
    let degree = 5.min(uv_points.len().saturating_sub(1));
    popkern_core::geometry::fitting::interpolate_curve(&uv_points, &parameters, degree)
        .ok()
        .map(Pcurve::from_nurbs)
        .unwrap_or_else(|| pcurve.clone())
}

fn orient_curve_to_edge_vertices(
    brep: &BRep,
    edge: EdgeHandle,
    curve: CurveOnEdge,
) -> Option<CurveOnEdge> {
    let direct = edge_endpoint_gap(brep, edge, &curve)?;
    let reversed = reversed_curve_on_edge(&curve)?;
    let reverse = edge_endpoint_gap(brep, edge, &reversed)?;
    Some(if reverse < direct { reversed } else { curve })
}

fn reversed_curve_on_edge(curve: &CurveOnEdge) -> Option<CurveOnEdge> {
    let samples = 17usize;
    let mut points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let t = curve.t_end + (curve.t_start - curve.t_end) * fraction;
        points.push(curve.geometry.point_at(t));
        parameters.push(fraction);
    }
    let degree = 5.min(points.len().saturating_sub(1));
    let geometry = popkern_core::geometry::fitting::interpolate_curve(&points, &parameters, degree)
        .ok()
        .map(CurveGeometry::Nurbs)
        .or_else(|| line_curve_from_endpoints(&points))?;
    Some(CurveOnEdge {
        geometry,
        t_start: 0.0,
        t_end: 1.0,
    })
}

fn edge_endpoint_gap(brep: &BRep, edge: EdgeHandle, curve: &CurveOnEdge) -> Option<f64> {
    let canonical = brep.edges.get(edge)?.half_edge;
    let start = brep
        .vertices
        .get(brep.half_edges.get(canonical)?.origin)?
        .position;
    let end = brep
        .vertices
        .get(
            brep.half_edges
                .get(brep.half_edges.get(canonical)?.twin)?
                .origin,
        )?
        .position;
    Some(
        curve
            .geometry
            .point_at(curve.t_start)
            .distance_to(&start)
            .max(curve.geometry.point_at(curve.t_end).distance_to(&end)),
    )
}

pub(super) fn curve_from_pcurve_in_edge_direction(
    brep: &BRep,
    edge: EdgeHandle,
    half_edge: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
) -> Option<CurveOnEdge> {
    let edge_data = brep.edges.get(edge)?;
    let forward = half_edge == edge_data.half_edge;
    let (pc_start, pc_end) = pcurve.parameter_range();
    let samples = 17usize;
    let mut points = Vec::with_capacity(samples);
    let mut parameters = Vec::with_capacity(samples);
    for index in 0..samples {
        let fraction = index as f64 / (samples - 1) as f64;
        let pcurve_fraction = if forward { fraction } else { 1.0 - fraction };
        let pcurve_t = pc_start + (pc_end - pc_start) * pcurve_fraction;
        let uv = pcurve.point_at(pcurve_t);
        points.push(surface.point_at(uv.x, uv.y));
        parameters.push(fraction);
    }
    let degree = 5.min(points.len().saturating_sub(1));
    let curve = popkern_core::geometry::fitting::interpolate_curve(&points, &parameters, degree)
        .ok()
        .map(CurveGeometry::Nurbs)
        .or_else(|| line_curve_from_endpoints(&points))?;
    Some(CurveOnEdge {
        geometry: curve,
        t_start: 0.0,
        t_end: 1.0,
    })
}

fn line_curve_from_endpoints(points: &[Point3]) -> Option<CurveGeometry> {
    let start = *points.first()?;
    let end = *points.last()?;
    CurveOnEdge::line_segment(start, end, 1.0e-12).map(|curve| curve.geometry)
}

pub(super) fn normalize_semantic_line_edge_curve(brep: &mut BRep, edge: EdgeHandle) {
    let Some(curve) = brep
        .edges
        .get_mut(edge)
        .and_then(|edge| edge.curve.as_mut())
    else {
        return;
    };
    if let Some(line) = curve
        .geometry
        .canonical_line_in_range((curve.t_start, curve.t_end), 1.0e-8)
    {
        curve.geometry = line;
    }
}

pub(super) fn solved_surface<'a>(
    solved: &'a ExactBoundarySolvedSurfaces,
    sector: &ExactBoundaryCapSectorInput,
) -> Option<&'a ExactBoundarySectorSurface> {
    solved
        .surfaces
        .iter()
        .find(|surface| surface.sector_index == sector.sector_index)
}
