//! Support pcurve projection and sample construction for exact-boundary ribbons.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};

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

fn unwrap_periodic_uv(
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

pub(super) fn compute_support_ribbon_sample(
    request: SupportRibbonSampleRequest<'_>,
) -> Result<ExactBoundaryRibbonSample, ExactBoundaryRibbonDeclineReason> {
    let SupportRibbonSampleRequest {
        boundary,
        other_boundary,
        surface,
        pcurve,
        parameter,
        pcurve_reversed,
        derivative_scale,
        budget,
    } = request;
    let pcurve_parameter = support_pcurve_parameter(pcurve, parameter, pcurve_reversed);
    let uv = pcurve.point_at(pcurve_parameter);
    let support_point = surface.point_at(uv.x, uv.y);
    let boundary_point = boundary.curve.point_at(parameter);
    let normal = surface.normal_at(uv.x, uv.y);
    let Some(support_normal) = normal.normalized() else {
        return Err(ExactBoundaryRibbonDeclineReason::DegenerateSupportNormal);
    };
    if normal.length() < budget.min_normal_length {
        return Err(ExactBoundaryRibbonDeclineReason::DegenerateSupportNormal);
    }

    let tangent = boundary.curve.derivative_at(parameter, 1);
    let tangent_in_plane = tangent - support_normal * tangent.dot(&support_normal);
    if tangent_in_plane.length() < budget.min_tangent_length {
        return Err(ExactBoundaryRibbonDeclineReason::TangentAmbiguity);
    }
    let Some(boundary_tangent) = tangent_in_plane.normalized() else {
        return Err(ExactBoundaryRibbonDeclineReason::TangentAmbiguity);
    };

    let cross = support_normal.cross(&boundary_tangent);
    if cross.length() < budget.min_cross_length {
        return Err(ExactBoundaryRibbonDeclineReason::CrossDerivativeDegenerate);
    }
    let Some(mut cross_unit) = cross.normalized() else {
        return Err(ExactBoundaryRibbonDeclineReason::CrossDerivativeDegenerate);
    };
    let inward = other_boundary.curve.point_at(0.5) - boundary_point;
    if cross_unit.dot(&inward) < 0.0 {
        cross_unit = cross_unit * -1.0;
    }
    let cross_derivative = cross_unit * derivative_scale;

    Ok(ExactBoundaryRibbonSample {
        parameter,
        support_uv: (uv.x, uv.y),
        boundary_point,
        support_point,
        boundary_tangent,
        support_normal,
        cross_derivative,
    })
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
