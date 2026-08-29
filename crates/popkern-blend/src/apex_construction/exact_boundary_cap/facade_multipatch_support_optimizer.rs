//! Bounded live support-seam optimization for selected multipatch strips.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use super::{
    CandidateSupportSeamReport, interior_normal_deficit_between_surfaces, pcurve_follows_boundary,
    selected_strip_support_seam_report,
};
use crate::apex_construction::exact_boundary_cap::multipatch::ExactBoundaryMultiPatchCandidate;
use crate::apex_construction::exact_boundary_cap::multipatch::budget::ExactBoundaryMultiPatchBudget;
use crate::apex_construction::exact_boundary_cap::pcurves::PreparedExactBoundaryPcurvePair;

pub(crate) fn optimize_selected_strip_live_support_seam(
    candidate: &mut ExactBoundaryMultiPatchCandidate,
    boundary: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
    budget: ExactBoundaryMultiPatchBudget,
) {
    let Some(mut best_report) =
        selected_strip_support_seam_report(candidate, boundary, support_pair)
    else {
        return;
    };
    if best_report.max_g1_deficit <= budget.support_normal_deficit {
        return;
    }
    let baseline_interface_limit = candidate
        .inner_patch
        .interface_g1_report
        .map(|report| report.tolerance)
        .unwrap_or(budget.interface_g1_deficit);
    for _ in 0..budget.max_interface_iterations {
        let Some(update) = support_row_update(candidate, boundary, support_pair, best_report)
        else {
            break;
        };
        let mut accepted = false;
        for damping in [0.75, 0.5, 0.25, 0.125, 0.0625] {
            let mut trial = candidate.clone();
            apply_support_row_update(&mut trial.strip_candidate.surface, &update, damping);
            if !strip_outer_and_interface_rows_preserved(candidate, &trial)
                || !strip_inner_interface_g1_preserved(&trial, baseline_interface_limit)
            {
                continue;
            }
            let Some(report) = selected_strip_support_seam_report(&trial, boundary, support_pair)
            else {
                continue;
            };
            if report.max_g1_deficit + 1.0e-12 < best_report.max_g1_deficit {
                best_report = report;
                candidate.strip_candidate = trial.strip_candidate;
                accepted = true;
                break;
            }
        }
        if best_report.max_g1_deficit <= budget.support_normal_deficit || !accepted {
            break;
        }
    }
}

#[derive(Debug, Clone)]
struct SupportRowUpdate {
    deltas: Vec<Vec3>,
}

fn support_row_update(
    candidate: &ExactBoundaryMultiPatchCandidate,
    boundary: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
    report: CandidateSupportSeamReport,
) -> Option<SupportRowUpdate> {
    let sample = support_seam_sample(candidate, boundary, support_pair, report.sample_parameter)?;
    let cap_u = interior_parameter_for_surface(
        &candidate.strip_candidate.surface,
        sample.cap_u,
        SurfaceDirection::U,
    );
    let cap_v = interior_parameter_for_surface(
        &candidate.strip_candidate.surface,
        sample.cap_v,
        SurfaceDirection::V,
    );
    let cap_tangent = candidate
        .strip_candidate
        .surface
        .partial_u(cap_u, cap_v)
        .normalized()?;
    let support_u = interior_parameter_for_surface(
        &support_pair.support_surface,
        sample.support_u,
        SurfaceDirection::U,
    );
    let support_v = interior_parameter_for_surface(
        &support_pair.support_surface,
        sample.support_v,
        SurfaceDirection::V,
    );
    let support_normal = support_pair
        .support_surface
        .normal_at(support_u, support_v)
        .normalized()?;
    let direction = support_normal.cross(&cap_tangent).normalized()?;
    let boundary_point = candidate
        .strip_candidate
        .surface
        .point_at(sample.cap_u, sample.cap_v);
    let current = support_row_point_at(&candidate.strip_candidate.surface, sample.cap_u)?;
    let width = current.distance_to(&boundary_point).max(1.0e-6);
    let target = boundary_point + direction * width;
    let delta = target - current;
    if delta.length() <= 1.0e-14 {
        return None;
    }
    let basis = support_row_rational_basis(&candidate.strip_candidate.surface, sample.cap_u)?;
    let scale = basis.iter().map(|value| value * value).sum::<f64>();
    if scale <= 1.0e-14 {
        return None;
    }
    Some(SupportRowUpdate {
        deltas: basis
            .into_iter()
            .map(|basis| delta * (basis / scale))
            .collect(),
    })
}

#[derive(Debug, Clone, Copy)]
struct SupportSeamSample {
    cap_u: f64,
    cap_v: f64,
    support_u: f64,
    support_v: f64,
}

fn support_seam_sample(
    candidate: &ExactBoundaryMultiPatchCandidate,
    boundary: &NurbsCurve,
    support_pair: &PreparedExactBoundaryPcurvePair,
    parameter: f64,
) -> Option<SupportSeamSample> {
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
    let support_uv = support_pcurve.point_at(support_lo + support_s * (support_hi - support_lo));
    Some(SupportSeamSample {
        cap_u: cap_uv.x,
        cap_v: cap_uv.y,
        support_u: support_uv.x,
        support_v: support_uv.y,
    })
}

fn apply_support_row_update(surface: &mut NurbsSurface, update: &SupportRowUpdate, damping: f64) {
    for (index, delta) in update.deltas.iter().enumerate() {
        let point_index = index * surface.num_v + 1;
        surface.control_points[point_index] =
            surface.control_points[point_index] + *delta * damping;
    }
}

fn strip_outer_and_interface_rows_preserved(
    baseline: &ExactBoundaryMultiPatchCandidate,
    trial: &ExactBoundaryMultiPatchCandidate,
) -> bool {
    let base = &baseline.strip_candidate.surface;
    let surface = &trial.strip_candidate.surface;
    if base.num_u != surface.num_u || base.num_v != surface.num_v {
        return false;
    }
    for index in 0..surface.num_u {
        for row in [0, surface.num_v - 1] {
            let point_index = index * surface.num_v + row;
            if surface.control_points[point_index].distance_to(&base.control_points[point_index])
                > baseline.strip_candidate.tolerance
            {
                return false;
            }
        }
    }
    true
}

fn strip_inner_interface_g1_preserved(
    candidate: &ExactBoundaryMultiPatchCandidate,
    tolerance: f64,
) -> bool {
    let Some(report) = strip_inner_interface_g1_report(candidate) else {
        return false;
    };
    report.max_g1_deficit <= tolerance
}

fn strip_inner_interface_g1_report(
    candidate: &ExactBoundaryMultiPatchCandidate,
) -> Option<CandidateSupportSeamReport> {
    let (t_start, t_end) = candidate.interface.curve.parameter_range();
    let mut report = CandidateSupportSeamReport {
        max_g1_deficit: 0.0,
        sample_parameter: t_start,
    };
    for index in 1..16 {
        let parameter = t_start + (t_end - t_start) * index as f64 / 16.0;
        let strip_uv = candidate.interface.strip_pcurve.point_at(parameter);
        let inner_uv = candidate.inner_patch.interface_pcurve.point_at(parameter);
        let deficit = interior_normal_deficit_between_surfaces(
            &candidate.strip_candidate.surface,
            strip_uv.x,
            strip_uv.y,
            &candidate.inner_patch.surface,
            inner_uv.x,
            inner_uv.y,
        )?;
        if deficit > report.max_g1_deficit {
            report = CandidateSupportSeamReport {
                max_g1_deficit: deficit,
                sample_parameter: parameter,
            };
        }
    }
    Some(report)
}

#[derive(Debug, Clone, Copy)]
enum SurfaceDirection {
    U,
    V,
}

fn interior_parameter_for_surface(
    surface: &impl SurfaceEval,
    value: f64,
    direction: SurfaceDirection,
) -> f64 {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let (period_u, period_v) = surface.parameter_periodicity();
    let (min, max, period) = match direction {
        SurfaceDirection::U => (u_min, u_max, period_u),
        SurfaceDirection::V => (v_min, v_max, period_v),
    };
    crate::chart_sampling::interior_parameter(
        crate::chart_sampling::wrap_periodic(value, min, max, period),
        min,
        max,
    )
}

fn support_row_point_at(surface: &NurbsSurface, parameter: f64) -> Option<Point3> {
    let basis = support_row_rational_basis(surface, parameter)?;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (index, rational_basis) in basis.into_iter().enumerate() {
        let point = surface.control_points[index * surface.num_v + 1];
        x += point.x * rational_basis;
        y += point.y * rational_basis;
        z += point.z * rational_basis;
    }
    Some(Point3::new(x, y, z))
}

fn support_row_rational_basis(surface: &NurbsSurface, parameter: f64) -> Option<Vec<f64>> {
    let mut weighted = Vec::with_capacity(surface.num_u);
    let mut denominator = 0.0;
    for index in 0..surface.num_u {
        let basis = basis_function(index, surface.degree_u, parameter, &surface.knots_u);
        let weight = surface.weights[index * surface.num_v + 1];
        let value = basis * weight;
        weighted.push(value);
        denominator += value;
    }
    if denominator.abs() <= 1.0e-14 || !denominator.is_finite() {
        return None;
    }
    Some(
        weighted
            .into_iter()
            .map(|value| value / denominator)
            .collect(),
    )
}

fn basis_function(index: usize, degree: usize, parameter: f64, knots: &[f64]) -> f64 {
    if degree == 0 {
        let is_last_span = index + 2 == knots.len();
        if (knots[index] <= parameter && parameter < knots[index + 1])
            || (is_last_span && (parameter - knots[index + 1]).abs() <= 1.0e-14)
        {
            return 1.0;
        }
        return 0.0;
    }
    let mut value = 0.0;
    let left = knots[index + degree] - knots[index];
    if left.abs() > 1.0e-14 {
        value +=
            (parameter - knots[index]) / left * basis_function(index, degree - 1, parameter, knots);
    }
    let right = knots[index + degree + 1] - knots[index + 1];
    if right.abs() > 1.0e-14 {
        value += (knots[index + degree + 1] - parameter) / right
            * basis_function(index + 1, degree - 1, parameter, knots);
    }
    value
}
