//! Support-authority boundary reauthoring for multipatch route attempts.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use crate::apex_construction::exact_boundary_cap::boundary::PreparedExactBoundaryCurve;
use crate::apex_construction::exact_boundary_cap::multipatch::report::ExactBoundaryMultiPatchDeclineReason;
use crate::apex_construction::exact_boundary_cap::pcurves::PreparedExactBoundaryPcurvePair;

const SUPPORT_BOUNDARY_HERMITE_SEGMENT_BUDGETS: [usize; 5] = [4, 8, 16, 32, 64];
const SUPPORT_BOUNDARY_FIT_SAMPLE_BUDGETS: [usize; 3] = [17, 33, 65];
const SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE: f64 = 5.0e-3;
const SUPPORT_TANGENT_PLANE_DEFICIT: f64 = 5.0e-2;

#[cfg(test)]
pub(crate) fn support_authority_boundary_curve(
    boundary: &PreparedExactBoundaryCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> Option<NurbsCurve> {
    support_authority_boundary_attempt(boundary, support_pair).candidate
}

pub(crate) fn support_authority_boundary_attempt(
    boundary: &PreparedExactBoundaryCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> SupportAuthorityBoundaryAttempt {
    match support_hermite_boundary_curve(boundary, support_pair) {
        ReauthoredBoundaryResult::Accepted(candidate) => SupportAuthorityBoundaryAttempt {
            candidate: Some(candidate),
            decline: None,
        },
        ReauthoredBoundaryResult::Declined(_) => {
            match sampled_support_authority_boundary_curve(boundary, support_pair) {
                ReauthoredBoundaryResult::Accepted(candidate) => SupportAuthorityBoundaryAttempt {
                    candidate: Some(candidate),
                    decline: None,
                },
                ReauthoredBoundaryResult::Declined(decline) => SupportAuthorityBoundaryAttempt {
                    candidate: None,
                    decline: Some(decline),
                },
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct SupportAuthorityBoundaryAttempt {
    pub(crate) candidate: Option<NurbsCurve>,
    pub(crate) decline: Option<SupportBoundaryReauthorDecline>,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct SupportBoundaryReauthorDecline {
    pub(crate) reason: ExactBoundaryMultiPatchDeclineReason,
    pub(crate) detail_code: &'static str,
    pub(crate) residual: f64,
    pub(crate) sample_parameter: f64,
    pub(crate) tolerance: f64,
}

impl SupportBoundaryReauthorDecline {
    fn candidate_failed(
        detail_code: &'static str,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self {
            reason: ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorCandidateFailed,
            detail_code,
            residual,
            sample_parameter,
            tolerance,
        }
    }
}

enum ReauthoredBoundaryResult {
    Accepted(NurbsCurve),
    Declined(SupportBoundaryReauthorDecline),
}

fn sampled_support_authority_boundary_curve(
    boundary: &PreparedExactBoundaryCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> ReauthoredBoundaryResult {
    if let ReauthoredBoundaryResult::Accepted(candidate) =
        support_hermite_boundary_curve(boundary, support_pair)
    {
        return ReauthoredBoundaryResult::Accepted(candidate);
    }
    let support_forward = pcurve_follows_boundary(
        &support_pair.support_surface,
        &support_pair.support_pcurve,
        boundary.curve.point_at(0.0),
        boundary.curve.point_at(1.0),
    );
    let (start, end) = support_pair.support_pcurve.parameter_range();
    let mut last_decline = reauthor_decline(
        "support_boundary_reauthor.candidate_uninitialized",
        f64::INFINITY,
        0.0,
        SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE,
    );
    for sample_count in SUPPORT_BOUNDARY_FIT_SAMPLE_BUDGETS {
        let mut points = Vec::with_capacity(sample_count);
        let mut params = Vec::with_capacity(sample_count);
        for index in 0..sample_count {
            let parameter = index as f64 / (sample_count - 1) as f64;
            let support_parameter = if support_forward {
                parameter
            } else {
                1.0 - parameter
            };
            let t = start + (end - start) * support_parameter;
            let uv = support_pair.support_pcurve.point_at(t);
            points.push(support_pair.support_surface.point_at(uv.x, uv.y));
            params.push(parameter);
        }
        let candidate = match derivative_aware_fit_curve(&points, &params, support_pair) {
            Ok(candidate) => candidate,
            Err(decline) => {
                last_decline = ReauthoredBoundaryResult::Declined(decline);
                continue;
            }
        };
        match certify_support_boundary_candidate(
            boundary,
            support_pair,
            &candidate,
            support_forward,
        ) {
            None => return ReauthoredBoundaryResult::Accepted(candidate),
            Some(decline) => last_decline = ReauthoredBoundaryResult::Declined(decline),
        }
    }
    last_decline
}

fn derivative_aware_fit_curve(
    points: &[Point3],
    params: &[f64],
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> Result<NurbsCurve, SupportBoundaryReauthorDecline> {
    let degree = 5.min(points.len().saturating_sub(1));
    let candidate = popkern_core::geometry::fitting::interpolate_curve(points, params, degree)
        .ok()
        .or_else(|| {
            popkern_core::geometry::fitting::fit_curve(points, params, degree, points.len()).ok()
        })
        .ok_or_else(|| {
            SupportBoundaryReauthorDecline::candidate_failed(
                "support_boundary_reauthor.sampled_fit_failed",
                f64::INFINITY,
                0.0,
                SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE,
            )
        })?;
    let tangent = tangent_support_normal_report(&candidate, support_pair);
    let tolerance = max_tangent_dot_for_deficit(SUPPORT_TANGENT_PLANE_DEFICIT);
    if tangent.residual > tolerance {
        return Err(SupportBoundaryReauthorDecline::candidate_failed(
            "support_boundary_reauthor.sampled_fit_tangent_unsatisfied",
            tangent.residual,
            tangent.sample_parameter,
            tolerance,
        ));
    }
    Ok(candidate)
}

fn support_hermite_boundary_curve(
    boundary: &PreparedExactBoundaryCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> ReauthoredBoundaryResult {
    let support_forward = pcurve_follows_boundary(
        &support_pair.support_surface,
        &support_pair.support_pcurve,
        boundary.curve.point_at(0.0),
        boundary.curve.point_at(1.0),
    );
    let mut last_decline = reauthor_decline(
        "support_boundary_reauthor.candidate_uninitialized",
        f64::INFINITY,
        0.0,
        SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE,
    );
    for segment_count in SUPPORT_BOUNDARY_HERMITE_SEGMENT_BUDGETS {
        let candidate = match support_hermite_boundary_candidate(
            support_pair,
            support_forward,
            segment_count,
        ) {
            Ok(candidate) => candidate,
            Err(decline) => {
                last_decline = ReauthoredBoundaryResult::Declined(decline);
                continue;
            }
        };
        match certify_support_boundary_candidate(
            boundary,
            support_pair,
            &candidate,
            support_forward,
        ) {
            None => return ReauthoredBoundaryResult::Accepted(candidate),
            Some(decline) => last_decline = ReauthoredBoundaryResult::Declined(decline),
        }
    }
    last_decline
}

fn support_hermite_boundary_candidate(
    support_pair: &PreparedExactBoundaryPcurvePair,
    support_forward: bool,
    segment_count: usize,
) -> Result<NurbsCurve, SupportBoundaryReauthorDecline> {
    let mut controls = Vec::with_capacity(segment_count * 3 + 1);
    for segment in 0..segment_count {
        let s0 = segment as f64 / segment_count as f64;
        let s1 = (segment + 1) as f64 / segment_count as f64;
        let sample0 =
            support_boundary_sample(support_pair, s0, support_forward).ok_or_else(|| {
                SupportBoundaryReauthorDecline::candidate_failed(
                    "support_boundary_reauthor.support_sample_unavailable",
                    f64::INFINITY,
                    s0,
                    SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE,
                )
            })?;
        let sample1 =
            support_boundary_sample(support_pair, s1, support_forward).ok_or_else(|| {
                SupportBoundaryReauthorDecline::candidate_failed(
                    "support_boundary_reauthor.support_sample_unavailable",
                    f64::INFINITY,
                    s1,
                    SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE,
                )
            })?;
        if segment == 0 {
            controls.push(sample0.point);
        }
        let scale = (s1 - s0) / 3.0;
        controls.push(sample0.point + sample0.tangent * scale);
        controls.push(sample1.point - sample1.tangent * scale);
        controls.push(sample1.point);
    }
    let mut knots = Vec::with_capacity(segment_count * 3 + 5);
    knots.extend(std::iter::repeat_n(0.0, 4));
    for segment in 1..segment_count {
        let knot = segment as f64 / segment_count as f64;
        knots.extend(std::iter::repeat_n(knot, 3));
    }
    knots.extend(std::iter::repeat_n(1.0, 4));
    Ok(NurbsCurve {
        degree: 3,
        knots,
        weights: vec![1.0; controls.len()],
        control_points: controls,
    })
}

fn reauthor_decline(
    detail_code: &'static str,
    residual: f64,
    sample_parameter: f64,
    tolerance: f64,
) -> ReauthoredBoundaryResult {
    ReauthoredBoundaryResult::Declined(SupportBoundaryReauthorDecline {
        reason: ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorCandidateFailed,
        detail_code,
        residual,
        sample_parameter,
        tolerance,
    })
}

#[derive(Debug, Clone, Copy)]
struct SupportBoundarySample {
    point: Point3,
    tangent: Vec3,
}

fn support_boundary_sample(
    support_pair: &PreparedExactBoundaryPcurvePair,
    parameter: f64,
    support_forward: bool,
) -> Option<SupportBoundarySample> {
    let support_parameter = if support_forward {
        parameter
    } else {
        1.0 - parameter
    };
    let (start, end) = support_pair.support_pcurve.parameter_range();
    let pcurve_t = start + (end - start) * support_parameter;
    let uv = support_pair.support_pcurve.point_at(pcurve_t);
    let uv_tangent = support_pair.support_pcurve.tangent_at(pcurve_t);
    let direction_sign = if support_forward { 1.0 } else { -1.0 };
    let tangent = support_pair.support_surface.partial_u(uv.x, uv.y) * uv_tangent.x
        + support_pair.support_surface.partial_v(uv.x, uv.y) * uv_tangent.y;
    Some(SupportBoundarySample {
        point: support_pair.support_surface.point_at(uv.x, uv.y),
        tangent: tangent * direction_sign * (end - start).abs(),
    })
}

fn support_boundary_endpoint_gap(
    boundary: &PreparedExactBoundaryCurve,
    candidate: &NurbsCurve,
) -> Option<f64> {
    let (start, end) = candidate.parameter_range();
    Some(
        candidate
            .point_at(start)
            .distance_to(&boundary.curve.point_at(0.0))
            .max(
                candidate
                    .point_at(end)
                    .distance_to(&boundary.curve.point_at(1.0)),
            ),
    )
}

#[derive(Debug, Clone, Copy)]
struct ResidualReport {
    residual: f64,
    sample_parameter: f64,
}

fn endpoint_gap_sample(boundary: &PreparedExactBoundaryCurve, candidate: &NurbsCurve) -> f64 {
    let (start, end) = candidate.parameter_range();
    let start_gap = candidate
        .point_at(start)
        .distance_to(&boundary.curve.point_at(0.0));
    let end_gap = candidate
        .point_at(end)
        .distance_to(&boundary.curve.point_at(1.0));
    if start_gap >= end_gap { 0.0 } else { 1.0 }
}

fn certify_support_boundary_candidate(
    boundary: &PreparedExactBoundaryCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
    candidate: &NurbsCurve,
    support_forward: bool,
) -> Option<SupportBoundaryReauthorDecline> {
    let endpoint_gap = support_boundary_endpoint_gap(boundary, candidate).unwrap_or(f64::INFINITY);
    if endpoint_gap > SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE {
        return Some(SupportBoundaryReauthorDecline {
            reason: ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorEndpointMismatch,
            detail_code:
                ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorEndpointMismatch.code(),
            residual: endpoint_gap,
            sample_parameter: endpoint_gap_sample(boundary, candidate),
            tolerance: SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE,
        });
    }
    let support_residual = support_boundary_residual(candidate, support_pair, support_forward);
    if support_residual.residual > SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE {
        return Some(SupportBoundaryReauthorDecline {
            reason: ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorSupportResidual,
            detail_code:
                ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorSupportResidual.code(),
            residual: support_residual.residual,
            sample_parameter: support_residual.sample_parameter,
            tolerance: SUPPORT_BOUNDARY_RESIDUAL_TOLERANCE,
        });
    }
    let tangent_residual = tangent_support_normal_report(candidate, support_pair);
    let tangent_tolerance = max_tangent_dot_for_deficit(SUPPORT_TANGENT_PLANE_DEFICIT);
    if tangent_residual.residual > tangent_tolerance {
        return Some(SupportBoundaryReauthorDecline {
            reason: ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorTangentResidual,
            detail_code:
                ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorTangentResidual.code(),
            residual: tangent_residual.residual,
            sample_parameter: tangent_residual.sample_parameter,
            tolerance: tangent_tolerance,
        });
    }
    None
}

fn support_boundary_residual(
    candidate: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
    support_forward: bool,
) -> ResidualReport {
    let (start, end) = candidate.parameter_range();
    let mut report = ResidualReport {
        residual: 0.0,
        sample_parameter: 0.0,
    };
    for index in 0..33 {
        let parameter = index as f64 / 32.0;
        let candidate_t = start + (end - start) * parameter;
        if let Some(sample) = support_boundary_sample(support_pair, parameter, support_forward) {
            let residual = candidate.point_at(candidate_t).distance_to(&sample.point);
            if residual > report.residual {
                report = ResidualReport {
                    residual,
                    sample_parameter: parameter,
                };
            }
        }
    }
    report
}

fn max_tangent_dot_for_deficit(deficit: f64) -> f64 {
    let alignment = 1.0 - deficit.clamp(0.0, 1.0);
    (1.0 - alignment * alignment).max(0.0).sqrt()
}

#[cfg(test)]
pub(crate) fn tangent_support_normal_dot(
    boundary: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> f64 {
    tangent_support_normal_report(boundary, support_pair).residual
}

fn tangent_support_normal_report(
    boundary: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
) -> ResidualReport {
    let support_forward = pcurve_follows_boundary(
        &support_pair.support_surface,
        &support_pair.support_pcurve,
        boundary.point_at(0.0),
        boundary.point_at(1.0),
    );
    let mut report = ResidualReport {
        residual: 0.0,
        sample_parameter: 0.0,
    };
    for index in 0..32 {
        let parameter = 0.05 + 0.90 * index as f64 / 31.0;
        let Some(boundary_tangent) = boundary.derivative_at(parameter, 1).normalized() else {
            continue;
        };
        let support_parameter = if support_forward {
            parameter
        } else {
            1.0 - parameter
        };
        let t = support_pair.support_pcurve.t_start
            + (support_pair.support_pcurve.t_end - support_pair.support_pcurve.t_start)
                * support_parameter;
        let uv = support_pair.support_pcurve.point_at(t);
        let Some(support_normal) = support_pair
            .support_surface
            .partial_u(uv.x, uv.y)
            .cross(&support_pair.support_surface.partial_v(uv.x, uv.y))
            .normalized()
        else {
            continue;
        };
        let residual = boundary_tangent.dot(&support_normal).abs();
        if residual > report.residual {
            report = ResidualReport {
                residual,
                sample_parameter: parameter,
            };
        }
    }
    report
}

fn pcurve_follows_boundary(
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
