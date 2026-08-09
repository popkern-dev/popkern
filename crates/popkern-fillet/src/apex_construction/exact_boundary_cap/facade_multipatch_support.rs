//! Support-authority evidence for exact-boundary multipatch route attempts.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::super::super::multipatch::ExactBoundaryMultiPatchCandidate;
use super::super::super::multipatch::budget::ExactBoundaryMultiPatchBudget;
use super::super::super::multipatch::report::ExactBoundaryMultiPatchDiagnostic;
use super::super::super::pcurves::PreparedExactBoundaryPcurvePair;
use super::super::super::ribbon::{ExactBoundaryRibbonSample, PreparedExactBoundaryRibbon};
use super::super::super::sector::ExactBoundaryCapSectorInput;

#[path = "facade_multipatch_support_optimizer.rs"]
mod optimizer;
#[path = "facade_multipatch_support_reauthor.rs"]
mod reauthor;

pub(super) use optimizer::optimize_selected_strip_live_support_seam;
pub(super) use reauthor::support_authority_boundary_attempt_with_control_limit;
#[cfg(test)]
use reauthor::{support_authority_boundary_curve, tangent_support_normal_dot};

pub(super) fn dense_support_samples(
    boundary: &NurbsCurve,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> Vec<ExactBoundaryRibbonSample> {
    let mut samples = Vec::new();
    let support_forward = pcurve_follows_boundary(
        &pair.support_surface,
        &pair.support_pcurve,
        boundary.point_at(0.0),
        boundary.point_at(1.0),
    );
    for index in 0..64 {
        let parameter = index as f64 / 63.0;
        if let Some(sample) =
            dense_support_sample(boundary, pair, ribbon, parameter, support_forward)
        {
            samples.push(sample);
        }
    }
    samples.sort_by(|a, b| {
        a.parameter
            .partial_cmp(&b.parameter)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    samples
}

pub(super) fn selected_strip_live_equivalent_diagnostic(
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
    boundary: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
    budget: ExactBoundaryMultiPatchBudget,
) -> Option<ExactBoundaryMultiPatchDiagnostic> {
    let report = selected_strip_support_seam_report(candidate, boundary, support_pair)?;
    (report.max_g1_deficit > budget.support_normal_deficit).then(|| {
        ExactBoundaryMultiPatchDiagnostic::support_strip_band(
            sector,
            &candidate.strip_inventory,
            report.max_g1_deficit,
            report.sample_parameter,
            budget.support_normal_deficit,
            Vec::new(),
        )
    })
}

fn dense_support_sample(
    boundary: &NurbsCurve,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
    parameter: f64,
    support_forward: bool,
) -> Option<ExactBoundaryRibbonSample> {
    let support_parameter = if support_forward {
        parameter
    } else {
        1.0 - parameter
    };
    let support_t = pair.support_pcurve.t_start
        + (pair.support_pcurve.t_end - pair.support_pcurve.t_start) * support_parameter;
    let support_uv = pair.support_pcurve.point_at(support_t);
    let support_point = pair.support_surface.point_at(support_uv.x, support_uv.y);
    let support_normal = pair
        .support_surface
        .partial_u(support_uv.x, support_uv.y)
        .cross(&pair.support_surface.partial_v(support_uv.x, support_uv.y))
        .normalized()?;
    let boundary_point = boundary.point_at(parameter);
    let boundary_tangent =
        support_boundary_tangent(pair, support_t, support_uv.x, support_uv.y, support_forward)?;
    let nearest = nearest_support_sample(&ribbon.samples, parameter)?;
    Some(ExactBoundaryRibbonSample {
        parameter,
        support_uv: (support_uv.x, support_uv.y),
        boundary_point,
        support_point,
        boundary_tangent,
        support_normal,
        cross_derivative: nearest.cross_derivative,
    })
}

fn support_boundary_tangent(
    pair: &PreparedExactBoundaryPcurvePair,
    support_t: f64,
    support_u: f64,
    support_v: f64,
    support_forward: bool,
) -> Option<Vec3> {
    let direction_sign = if support_forward { 1.0 } else { -1.0 };
    let support_tangent = pair.support_pcurve.tangent_at(support_t);
    (pair.support_surface.partial_u(support_u, support_v) * support_tangent.x
        + pair.support_surface.partial_v(support_u, support_v) * support_tangent.y)
        .normalized()
        .map(|tangent| tangent * direction_sign)
}

fn nearest_support_sample(
    samples: &[ExactBoundaryRibbonSample],
    parameter: f64,
) -> Option<ExactBoundaryRibbonSample> {
    samples.iter().copied().min_by(|a, b| {
        (a.parameter - parameter)
            .abs()
            .partial_cmp(&(b.parameter - parameter).abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    })
}

#[derive(Debug, Clone, Copy)]
pub(super) struct CandidateSupportSeamReport {
    pub(super) max_g1_deficit: f64,
    pub(super) sample_parameter: f64,
}

pub(super) fn selected_strip_support_seam_report(
    candidate: &ExactBoundaryMultiPatchCandidate,
    boundary: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> Option<CandidateSupportSeamReport> {
    let cap_pcurve = &candidate.strip_candidate.outer_pcurve;
    let support_pcurve = &support_pair.support_pcurve;
    let (cap_lo, cap_hi) = cap_pcurve.parameter_range();
    let (support_lo, support_hi) = support_pcurve.parameter_range();
    let boundary_start = boundary.point_at(0.0);
    let boundary_end = boundary.point_at(1.0);
    let cap_forward = pcurve_follows_boundary(
        &candidate.strip_candidate.surface,
        cap_pcurve,
        boundary_start,
        boundary_end,
    );
    let support_forward = pcurve_follows_boundary(
        &support_pair.support_surface,
        support_pcurve,
        boundary_start,
        boundary_end,
    );

    let mut report = CandidateSupportSeamReport {
        max_g1_deficit: 0.0,
        sample_parameter: 0.0,
    };
    for index in 0..100 {
        let parameter = 0.05 + 0.90 * index as f64 / 99.0;
        let cap_s = if cap_forward {
            parameter
        } else {
            1.0 - parameter
        };
        let support_s = if support_forward {
            parameter
        } else {
            1.0 - parameter
        };
        let cap_uv = cap_pcurve.point_at(cap_lo + cap_s * (cap_hi - cap_lo));
        let support_uv =
            support_pcurve.point_at(support_lo + support_s * (support_hi - support_lo));
        let Some(deficit) = interior_normal_deficit_between_surfaces(
            &candidate.strip_candidate.surface,
            cap_uv.x,
            cap_uv.y,
            &support_pair.support_surface,
            support_uv.x,
            support_uv.y,
        ) else {
            continue;
        };
        if deficit > report.max_g1_deficit {
            report = CandidateSupportSeamReport {
                max_g1_deficit: deficit,
                sample_parameter: parameter,
            };
        }
    }
    Some(report)
}

pub(super) fn normal_deficit_between_surfaces(
    cap_surface: &impl SurfaceEval,
    cap_u: f64,
    cap_v: f64,
    support_surface: &impl SurfaceEval,
    support_u: f64,
    support_v: f64,
) -> Option<f64> {
    let cap_normal = cap_surface
        .partial_u(cap_u, cap_v)
        .cross(&cap_surface.partial_v(cap_u, cap_v));
    let support_normal = support_surface
        .partial_u(support_u, support_v)
        .cross(&support_surface.partial_v(support_u, support_v));
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length <= 1.0e-300
        || support_length <= 1.0e-300
        || !cap_length.is_finite()
        || !support_length.is_finite()
    {
        return None;
    }
    Some(1.0 - (cap_normal.dot(&support_normal) / (cap_length * support_length)).abs())
}

pub(super) fn interior_normal_deficit_between_surfaces(
    cap_surface: &impl SurfaceEval,
    cap_u: f64,
    cap_v: f64,
    support_surface: &impl SurfaceEval,
    support_u: f64,
    support_v: f64,
) -> Option<f64> {
    let ((cap_u_min, cap_u_max), (cap_v_min, cap_v_max)) = cap_surface.parameter_range();
    let ((support_u_min, support_u_max), (support_v_min, support_v_max)) =
        support_surface.parameter_range();
    let (cap_period_u, cap_period_v) = cap_surface.parameter_periodicity();
    let (support_period_u, support_period_v) = support_surface.parameter_periodicity();
    normal_deficit_between_surfaces(
        cap_surface,
        interior_parameter(
            wrap_periodic(cap_u, cap_u_min, cap_u_max, cap_period_u),
            cap_u_min,
            cap_u_max,
        ),
        interior_parameter(
            wrap_periodic(cap_v, cap_v_min, cap_v_max, cap_period_v),
            cap_v_min,
            cap_v_max,
        ),
        support_surface,
        interior_parameter(
            wrap_periodic(support_u, support_u_min, support_u_max, support_period_u),
            support_u_min,
            support_u_max,
        ),
        interior_parameter(
            wrap_periodic(support_v, support_v_min, support_v_max, support_period_v),
            support_v_min,
            support_v_max,
        ),
    )
}

/// Periodic directions WRAP into the domain instead of clamping. A
/// support pcurve legally chained onto an adjacent period branch
/// (u − 2π on a cylinder chart) evaluates exactly, but clamping such
/// a parameter to the domain edge samples the normal up to a full
/// period away (the PK-125 lesson, recurred here: this check clamped
/// u = 7.78 to 2π on an r = 0.2 strip — 1.50 rad ≈ 86° off — and
/// declined an e-4-tangent candidate at 0.888).
fn wrap_periodic(value: f64, min: f64, max: f64, period: Option<f64>) -> f64 {
    let Some(period) = period else {
        return value;
    };
    if !(value.is_finite() && min.is_finite() && max.is_finite()) || period <= 0.0 {
        return value;
    }
    let mut wrapped = value;
    while wrapped < min {
        wrapped += period;
    }
    while wrapped > max {
        wrapped -= period;
    }
    wrapped
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * 1.0e-5).max(1.0e-5);
    value.clamp(min + margin, max - margin)
}

pub(super) fn pcurve_follows_boundary(
    surface: &impl SurfaceEval,
    pcurve: &Pcurve,
    boundary_start: Point3,
    boundary_end: Point3,
) -> bool {
    let (lo, hi) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(lo);
    let end_uv = pcurve.point_at(hi);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    boundary_start.distance_to(&start) + boundary_end.distance_to(&end)
        <= boundary_start.distance_to(&end) + boundary_end.distance_to(&start)
}

#[cfg(test)]
#[path = "facade_multipatch_support_tests.rs"]
mod tests;
