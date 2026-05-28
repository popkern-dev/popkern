//! Support-pcurve authority selection for exact-boundary curves.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::BRep;

use super::super::sector::ExactBoundaryOuterSeam;

const SUPPORT_SOURCE_SAMPLES: usize = 17;
const SUPPORT_SOURCE_TOLERANCE: f64 = 1.0e-6;

pub(super) fn validated_support_pcurve_source_curve(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    edge_curve: &NurbsCurve,
) -> Option<NurbsCurve> {
    let support_curve = support_pcurve_source_curve(brep, seam)
        .or_else(|| projected_support_source_curve(brep, seam, edge_curve))?;
    let (origin, destination) = super::cap_half_edge_vertex_span(brep, seam.cap_half_edge)?;
    let direct_endpoint_gap = super::endpoint_orientation_gap(&support_curve, origin, destination);
    let reverse_endpoint_gap = super::endpoint_orientation_gap(
        &super::reverse_nurbs_curve(&support_curve),
        origin,
        destination,
    );
    if direct_endpoint_gap > super::TOPOLOGY_ENDPOINT_TOLERANCE
        && reverse_endpoint_gap > super::TOPOLOGY_ENDPOINT_TOLERANCE
    {
        return None;
    }
    let support_surface = support_surface_for_seam(brep, seam)?;
    let support_range = support_curve.parameter_range();
    let support_deviation = edge_surface_deviation(
        &support_curve,
        support_range,
        support_surface,
        SUPPORT_SOURCE_SAMPLES,
    )
    .max_deviation;
    if support_deviation > SUPPORT_SOURCE_TOLERANCE {
        return None;
    }
    Some(support_curve)
}

fn support_pcurve_source_curve(brep: &BRep, seam: &ExactBoundaryOuterSeam) -> Option<NurbsCurve> {
    let support_half_edge = brep.half_edges.get(seam.support_half_edge)?;
    let pcurve = support_half_edge.pcurve.as_ref()?;
    let support_surface = support_surface_for_seam(brep, seam)?;
    let (t0, t1) = pcurve.parameter_range();
    let mut points = Vec::with_capacity(SUPPORT_SOURCE_SAMPLES);
    let mut params = Vec::with_capacity(SUPPORT_SOURCE_SAMPLES);
    for index in 0..SUPPORT_SOURCE_SAMPLES {
        let parameter = index as f64 / (SUPPORT_SOURCE_SAMPLES - 1) as f64;
        let pcurve_t = t0 + (t1 - t0) * parameter;
        let uv = pcurve.point_at(pcurve_t);
        points.push(support_surface.point_at(uv.x, uv.y));
        params.push(parameter);
    }
    let degree = 9.min(points.len().saturating_sub(1));
    popkern_core::geometry::fitting::fit_curve(&points, &params, degree, points.len()).ok()
}

fn projected_support_source_curve(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    edge_curve: &NurbsCurve,
) -> Option<NurbsCurve> {
    let support_surface = support_surface_for_seam(brep, seam)?;
    let edge_range = edge_curve.parameter_range();
    let mut points = Vec::with_capacity(SUPPORT_SOURCE_SAMPLES);
    let mut params = Vec::with_capacity(SUPPORT_SOURCE_SAMPLES);
    for index in 0..SUPPORT_SOURCE_SAMPLES {
        let parameter = index as f64 / (SUPPORT_SOURCE_SAMPLES - 1) as f64;
        let t = edge_range.0 + (edge_range.1 - edge_range.0) * parameter;
        let point = edge_curve.point_at(t);
        let (u, v) = support_surface.closest_point(&point)?;
        points.push(support_surface.point_at(u, v));
        params.push(parameter);
    }
    let degree = 9.min(points.len().saturating_sub(1));
    popkern_core::geometry::fitting::fit_curve(&points, &params, degree, points.len()).ok()
}

fn support_surface_for_seam<'a>(
    brep: &'a BRep,
    seam: &ExactBoundaryOuterSeam,
) -> Option<&'a dyn SurfaceEval> {
    let support_half_edge = brep.half_edges.get(seam.support_half_edge)?;
    let support_face = brep.loops.get(support_half_edge.loop_)?.face;
    Some(&brep.faces.get(support_face)?.surface.as_ref()?.geometry)
}
