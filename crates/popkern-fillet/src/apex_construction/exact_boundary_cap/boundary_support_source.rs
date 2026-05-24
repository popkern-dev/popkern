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
    let support_curve = support_pcurve_source_curve(brep, seam)?;
    let (origin, destination) = super::cap_half_edge_vertex_span(brep, seam.cap_half_edge)?;
    if super::endpoint_orientation_gap(&support_curve, origin, destination)
        > super::TOPOLOGY_ENDPOINT_TOLERANCE
        && super::endpoint_orientation_gap(
            &super::reverse_nurbs_curve(&support_curve),
            origin,
            destination,
        ) > super::TOPOLOGY_ENDPOINT_TOLERANCE
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
    let edge_range = edge_curve.parameter_range();
    let edge_deviation = edge_surface_deviation(
        edge_curve,
        edge_range,
        support_surface,
        SUPPORT_SOURCE_SAMPLES,
    )
    .max_deviation;
    let cap_pcurve_missing = brep
        .half_edges
        .get(seam.cap_half_edge)
        .is_some_and(|half_edge| half_edge.pcurve.is_none());
    (edge_deviation > support_deviation.max(SUPPORT_SOURCE_TOLERANCE) || cap_pcurve_missing)
        .then_some(support_curve)
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

fn support_surface_for_seam<'a>(
    brep: &'a BRep,
    seam: &ExactBoundaryOuterSeam,
) -> Option<&'a dyn SurfaceEval> {
    let support_half_edge = brep.half_edges.get(seam.support_half_edge)?;
    let support_face = brep.loops.get(support_half_edge.loop_)?.face;
    Some(&brep.faces.get(support_face)?.surface.as_ref()?.geometry)
}
