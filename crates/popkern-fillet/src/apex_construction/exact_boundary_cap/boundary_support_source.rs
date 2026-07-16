//! Support-pcurve authority selection for exact-boundary curves.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::BRep;

use super::super::sector::ExactBoundaryOuterSeam;

const SUPPORT_SOURCE_SAMPLES: usize = 17;
const SUPPORT_SOURCE_TOLERANCE: f64 = 1.0e-6;
/// Hard ceiling for replacing an off-support boundary curve with the
/// support-pcurve authority. Boundary-conditioned fillet marches move
/// band curves by up to ~5e-4 (see corner.rs anchor reconciliation);
/// anything farther indicates a topology mismatch, not march error.
const SUPPORT_SOURCE_MAX_REPLACEMENT_DEVIATION: f64 = 5.0e-4;

pub(super) fn validated_support_pcurve_source_curve(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    edge_curve: &NurbsCurve,
) -> Option<NurbsCurve> {
    let support_curve = support_pcurve_source_curve(brep, seam)
        .or_else(|| projected_support_source_curve(brep, seam, edge_curve))?;
    let (origin, destination) =
        super::boundary_curve_orientation::cap_half_edge_vertex_span(brep, seam.cap_half_edge)?;
    let direct_endpoint_gap = super::boundary_curve_orientation::endpoint_orientation_gap(
        &support_curve,
        origin,
        destination,
    );
    let reverse_endpoint_gap = super::boundary_curve_orientation::endpoint_orientation_gap(
        &super::boundary_curve_orientation::reverse_nurbs_curve(&support_curve),
        origin,
        destination,
    );
    if direct_endpoint_gap > super::TOPOLOGY_ENDPOINT_TOLERANCE
        && reverse_endpoint_gap > super::TOPOLOGY_ENDPOINT_TOLERANCE
    {
        return None;
    }
    let support_range = support_curve.parameter_range();
    let edge_range = edge_curve.parameter_range();
    let edge_deviation =
        curve_pair_deviation(&support_curve, support_range, edge_curve, edge_range);
    // The authority may legitimately differ from the raw curve by about
    // the raw curve's own measured off-support error (that error is what
    // it corrects); beyond that scale the candidate is rejected as a
    // topology mismatch rather than a repair.
    let raw_off_support = raw_curve_support_deviation(brep, seam, edge_curve)?;
    let replacement_allowance = raw_off_support
        .mul_add(2.0, SUPPORT_SOURCE_TOLERANCE)
        .min(SUPPORT_SOURCE_MAX_REPLACEMENT_DEVIATION);
    if edge_deviation > replacement_allowance {
        return None;
    }
    validate_support_curve_on_seam(brep, seam, &support_curve)?;
    Some(support_curve)
}

/// Measured deviation of a candidate boundary curve from the seam's
/// support surface.
pub(super) fn raw_curve_support_deviation(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    curve: &NurbsCurve,
) -> Option<f64> {
    Some(
        edge_surface_deviation(
            curve,
            curve.parameter_range(),
            support_surface_for_seam(brep, seam)?,
            SUPPORT_SOURCE_SAMPLES,
        )
        .max_deviation,
    )
}

/// Whether a prepared boundary curve is measurably off its support
/// surface and should defer to the support-pcurve authority.
pub(super) fn raw_curve_needs_support_authority(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    curve: &NurbsCurve,
) -> bool {
    raw_curve_support_deviation(brep, seam, curve)
        .is_some_and(|deviation| deviation > SUPPORT_SOURCE_TOLERANCE)
}

pub(super) fn validated_projected_support_source_curve(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    topology_curve: &NurbsCurve,
) -> Option<NurbsCurve> {
    let support_curve = projected_support_source_curve(brep, seam, topology_curve)?;
    validate_support_curve_on_seam(brep, seam, &support_curve)?;
    Some(support_curve)
}

fn validate_support_curve_on_seam(
    brep: &BRep,
    seam: &ExactBoundaryOuterSeam,
    support_curve: &NurbsCurve,
) -> Option<()> {
    let (origin, destination) =
        super::boundary_curve_orientation::cap_half_edge_vertex_span(brep, seam.cap_half_edge)?;
    let direct_endpoint_gap = super::boundary_curve_orientation::endpoint_orientation_gap(
        support_curve,
        origin,
        destination,
    );
    let reverse_endpoint_gap = super::boundary_curve_orientation::endpoint_orientation_gap(
        &super::boundary_curve_orientation::reverse_nurbs_curve(support_curve),
        origin,
        destination,
    );
    if direct_endpoint_gap > super::TOPOLOGY_ENDPOINT_TOLERANCE
        && reverse_endpoint_gap > super::TOPOLOGY_ENDPOINT_TOLERANCE
    {
        return None;
    }
    let support_deviation = edge_surface_deviation(
        support_curve,
        support_curve.parameter_range(),
        support_surface_for_seam(brep, seam)?,
        SUPPORT_SOURCE_SAMPLES,
    )
    .max_deviation;
    if support_deviation > SUPPORT_SOURCE_TOLERANCE {
        return None;
    }
    Some(())
}

fn curve_pair_deviation(
    candidate: &NurbsCurve,
    candidate_range: (f64, f64),
    authority: &NurbsCurve,
    authority_range: (f64, f64),
) -> f64 {
    let mut worst = 0.0_f64;
    for index in 0..SUPPORT_SOURCE_SAMPLES {
        let parameter = index as f64 / (SUPPORT_SOURCE_SAMPLES - 1) as f64;
        let candidate_t = candidate_range.0 + (candidate_range.1 - candidate_range.0) * parameter;
        let authority_t = authority_range.0 + (authority_range.1 - authority_range.0) * parameter;
        worst = worst.max(
            candidate
                .point_at(candidate_t)
                .distance_to(&authority.point_at(authority_t)),
        );
    }
    worst
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
