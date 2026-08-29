//! Support pcurve projection and sample construction for exact-boundary ribbons.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use popkern_core::geometry::{Point3, Vec3};

use super::boundary::PreparedExactBoundaryCurve;
use super::ribbon::{
    ExactBoundaryRibbonBudget, ExactBoundaryRibbonDeclineReason, ExactBoundaryRibbonSample,
    SupportRibbonSampleRequest,
};

pub(super) struct SupportPcurveOrientation {
    pub(super) reversed: bool,
    pub(super) max_residual: f64,
}

pub(super) fn projected_support_pcurve(
    boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    budget: ExactBoundaryRibbonBudget,
) -> Option<Pcurve> {
    let count = sample_count(budget).max(5);
    let (period_u, period_v) = surface.periodicity();
    let mut previous_uv = None;
    let mut uv_points = Vec::with_capacity(count);
    let mut parameters = Vec::with_capacity(count);
    for index in 0..count {
        let parameter = index as f64 / (count - 1) as f64;
        let point = boundary.curve.point_at(parameter);
        let (mut u, mut v) = if let Some((previous_u, previous_v)) = previous_uv {
            project_point_to_uv_near(surface, &point, previous_u, previous_v)?
        } else {
            project_point_to_uv(surface, &point)?
        };
        unwrap_periodic_uv(&mut u, &mut v, previous_uv, period_u, period_v);
        previous_uv = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
        parameters.push(parameter);
    }
    let degree = 3.min(uv_points.len().saturating_sub(1));
    popkern_core::geometry::fitting::fit_curve(&uv_points, &parameters, degree, uv_points.len())
        .ok()
        .map(Pcurve::from_nurbs)
}

pub(in crate::apex_construction::exact_boundary_cap) fn unwrap_periodic_uv(
    u: &mut f64,
    v: &mut f64,
    previous_uv: Option<(f64, f64)>,
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    let Some((previous_u, previous_v)) = previous_uv else {
        return;
    };
    if let Some(period) = period_u {
        while *u - previous_u > period * 0.5 {
            *u -= period;
        }
        while previous_u - *u > period * 0.5 {
            *u += period;
        }
    }
    if let Some(period) = period_v {
        while *v - previous_v > period * 0.5 {
            *v -= period;
        }
        while previous_v - *v > period * 0.5 {
            *v += period;
        }
    }
}

pub(super) fn choose_support_pcurve_orientation(
    boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    budget: ExactBoundaryRibbonBudget,
) -> SupportPcurveOrientation {
    let forward = max_support_residual(boundary, surface, pcurve, false, budget);
    let reversed = max_support_residual(boundary, surface, pcurve, true, budget);
    if reversed <= forward {
        SupportPcurveOrientation {
            reversed: true,
            max_residual: reversed,
        }
    } else {
        SupportPcurveOrientation {
            reversed: false,
            max_residual: forward,
        }
    }
}

fn max_support_residual(
    boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    reversed: bool,
    budget: ExactBoundaryRibbonBudget,
) -> f64 {
    sample_parameters(budget)
        .map(|parameter| {
            let boundary_point = boundary.curve.point_at(parameter);
            let uv = pcurve.point_at(support_pcurve_parameter(pcurve, parameter, reversed));
            let support_point = surface.point_at(uv.x, uv.y);
            boundary_point.distance_to(&support_point)
        })
        .fold(0.0, f64::max)
}

#[cfg(test)]
pub(super) fn compute_support_ribbon_sample(
    request: SupportRibbonSampleRequest<'_>,
) -> Result<ExactBoundaryRibbonSample, ExactBoundaryRibbonDeclineReason> {
    let parameter = request.parameter;
    compute_support_ribbon_samples(request)?
        .into_iter()
        .min_by(|a, b| {
            (a.parameter - parameter)
                .abs()
                .partial_cmp(&(b.parameter - parameter).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .ok_or(ExactBoundaryRibbonDeclineReason::TangentAmbiguity)
}

pub(super) fn compute_support_ribbon_samples(
    request: SupportRibbonSampleRequest<'_>,
) -> Result<Vec<ExactBoundaryRibbonSample>, ExactBoundaryRibbonDeclineReason> {
    let SupportRibbonSampleRequest {
        boundary,
        other_boundary,
        surface,
        pcurve,
        parameter: _requested_parameter,
        pcurve_reversed,
        derivative_scale,
        budget,
    } = request;
    let trace = support_trace_points(boundary, surface, pcurve, pcurve_reversed, budget)?;
    let mut samples = Vec::with_capacity(trace.len());
    // Inward disambiguation is decided ONCE per ribbon and then
    // propagated by field continuity. The old per-sample test
    // (`dot(other_boundary.point_at(0.5) - sample_point)`) flipped
    // the field's sign mid-seam whenever the reference point sat on
    // the wrong side of a sample — measured on pyramid4_uniform,
    // whose seam field is the constant cylinder axis: the reference
    // (the point-like face boundary's midpoint = the corner itself)
    // flipped D at t ≈ 0.3, and the solved surface smoothed the 180°
    // discontinuity into the symmetric cap-vs-cylinder normal error
    // that grew to 4.8e-1 at the seam ends (PK-125).
    let mut previous_cross: Option<Vec3> = None;
    let inward_reference = other_boundary.curve.point_at(0.5) - boundary.curve.point_at(0.5);
    for index in 0..trace.len() {
        let point = &trace[index];
        let boundary_tangent =
            support_trace_tangent(&trace, index, boundary, other_boundary, surface, budget)?
                .normalized()
                .ok_or(ExactBoundaryRibbonDeclineReason::TangentAmbiguity)?;
        let cross = point.support_normal.cross(&boundary_tangent);
        if cross.length() < budget.min_cross_length {
            return Err(ExactBoundaryRibbonDeclineReason::CrossDerivativeDegenerate);
        }
        let Some(mut cross_unit) = cross.normalized() else {
            return Err(ExactBoundaryRibbonDeclineReason::CrossDerivativeDegenerate);
        };
        let orientation_sign = match previous_cross {
            Some(previous) => cross_unit.dot(&previous),
            None => cross_unit.dot(&inward_reference),
        };
        if orientation_sign < 0.0 {
            cross_unit = cross_unit * -1.0;
        }
        previous_cross = Some(cross_unit);
        samples.push(ExactBoundaryRibbonSample {
            parameter: point.parameter,
            support_uv: point.support_uv,
            boundary_point: point.boundary_point,
            support_point: point.support_point,
            boundary_tangent,
            support_normal: point.support_normal,
            cross_derivative: cross_unit * derivative_scale,
        });
    }
    Ok(samples)
}

pub(super) fn local_sector_width(
    boundary: &PreparedExactBoundaryCurve,
    other_boundary: &PreparedExactBoundaryCurve,
) -> f64 {
    boundary
        .curve
        .point_at(0.5)
        .distance_to(&other_boundary.curve.point_at(0.5))
}

pub(super) fn sample_count(budget: ExactBoundaryRibbonBudget) -> usize {
    budget.sample_count.max(2)
}

pub(super) fn sample_parameters(
    budget: ExactBoundaryRibbonBudget,
) -> impl Clone + ExactSizeIterator<Item = f64> {
    let count = sample_count(budget);
    (0..count).map(move |index| index as f64 / (count - 1) as f64)
}

fn support_pcurve_parameter(pcurve: &Pcurve, parameter: f64, reversed: bool) -> f64 {
    let mapped = if reversed { 1.0 - parameter } else { parameter };
    pcurve.t_start + (pcurve.t_end - pcurve.t_start) * mapped
}

struct SupportTracePoint {
    parameter: f64,
    support_uv: (f64, f64),
    boundary_point: Point3,
    support_point: Point3,
    support_normal: Vec3,
}

fn support_trace_points(
    boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    pcurve_reversed: bool,
    budget: ExactBoundaryRibbonBudget,
) -> Result<Vec<SupportTracePoint>, ExactBoundaryRibbonDeclineReason> {
    sample_parameters(budget)
        .map(|parameter| {
            let pcurve_parameter = support_pcurve_parameter(pcurve, parameter, pcurve_reversed);
            let uv = pcurve.point_at(pcurve_parameter);
            let support_point = surface.point_at(uv.x, uv.y);
            let boundary_point = boundary.curve.point_at(parameter);
            let normal = surface.normal_at(uv.x, uv.y);
            if normal.length() < budget.min_normal_length {
                return Err(ExactBoundaryRibbonDeclineReason::DegenerateSupportNormal);
            }
            let Some(support_normal) = normal.normalized() else {
                return Err(ExactBoundaryRibbonDeclineReason::DegenerateSupportNormal);
            };
            Ok(SupportTracePoint {
                parameter,
                support_uv: (uv.x, uv.y),
                boundary_point,
                support_point,
                support_normal,
            })
        })
        .collect()
}

fn support_trace_tangent(
    trace: &[SupportTracePoint],
    index: usize,
    boundary: &PreparedExactBoundaryCurve,
    other_boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    budget: ExactBoundaryRibbonBudget,
) -> Result<Vec3, ExactBoundaryRibbonDeclineReason> {
    let point = &trace[index];
    let tangent = boundary.curve.derivative_at(point.parameter, 1);
    let tangent_in_plane = tangent - point.support_normal * tangent.dot(&point.support_normal);
    if tangent_in_plane.length() >= budget.min_tangent_length {
        return Ok(tangent_in_plane);
    }
    let Some(chord) = support_trace_chord(trace, index) else {
        return support_frame_tangent(surface, point, other_boundary, budget)
            .ok_or(ExactBoundaryRibbonDeclineReason::TangentAmbiguity);
    };
    let chord_in_plane = chord - point.support_normal * chord.dot(&point.support_normal);
    if chord_in_plane.length() < budget.min_tangent_length {
        return support_frame_tangent(surface, point, other_boundary, budget)
            .ok_or(ExactBoundaryRibbonDeclineReason::TangentAmbiguity);
    }
    Ok(chord_in_plane)
}

fn support_frame_tangent(
    surface: &SurfaceGeometry,
    point: &SupportTracePoint,
    other_boundary: &PreparedExactBoundaryCurve,
    budget: ExactBoundaryRibbonBudget,
) -> Option<Vec3> {
    let inward = other_boundary.curve.point_at(0.5) - point.boundary_point;
    let inward_in_plane = inward - point.support_normal * inward.dot(&point.support_normal);
    if inward_in_plane.length() >= budget.min_tangent_length {
        return Some(inward_in_plane.cross(&point.support_normal));
    }
    support_surface_tangent(surface, point, budget)
}

fn support_surface_tangent(
    surface: &SurfaceGeometry,
    point: &SupportTracePoint,
    budget: ExactBoundaryRibbonBudget,
) -> Option<Vec3> {
    let u = surface.partial_u(point.support_uv.0, point.support_uv.1);
    let v = surface.partial_v(point.support_uv.0, point.support_uv.1);
    let u = u - point.support_normal * u.dot(&point.support_normal);
    let v = v - point.support_normal * v.dot(&point.support_normal);
    let tangent = if u.length() >= v.length() { u } else { v };
    (tangent.length() >= budget.min_tangent_length).then_some(tangent)
}

fn support_trace_chord(trace: &[SupportTracePoint], index: usize) -> Option<Vec3> {
    let current = trace.get(index)?;
    if let Some(next) = trace.get(index + 1) {
        let chord = next.support_point - current.support_point;
        if chord.length() > 0.0 {
            return Some(chord);
        }
    }
    if index > 0
        && let Some(previous) = trace.get(index - 1)
    {
        let chord = current.support_point - previous.support_point;
        if chord.length() > 0.0 {
            return Some(chord);
        }
    }
    None
}
