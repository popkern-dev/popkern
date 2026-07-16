//! Post-pcurve support-normal strip refinement for exact-boundary solves.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};

use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::ribbon::{
    ExactBoundaryRibbonSample, PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons,
};
use super::super::sector::ExactBoundaryOuterSide;
use super::ExactBoundarySolvedSurfaces;
use super::solve_support_pcurve_optimize::optimize_pair_controls;
use super::solve_support_pcurve_residual::worst_cap_support_position_residual;
use super::solve_support_strip::{
    MAX_ITERATIONS, MAX_SUPPORT_NORMAL_DEFICIT, StripSide, correct_strip_point_least_norm,
    interior_rational_strip_basis, scaled_target_derivatives, strip_curve,
    support_normal_deficit_at,
};

const MAX_REFINED_BOUNDARY_RESIDUAL: f64 = 1.0e-6;

pub(crate) fn refine_support_normal_strips_with_pcurves(
    solved: &mut ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    ribbons: &PreparedExactBoundaryRibbons,
) {
    let baseline_deficit = first_pcurve_decline_deficit(solved, pcurves, ribbons)
        .unwrap_or_else(|| worst_pcurve_support_normal_deficit(solved, pcurves, ribbons));
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace exact-boundary/pcurve-strip] baseline deficit {baseline_deficit:.4e}"
        );
    }
    if baseline_deficit <= MAX_SUPPORT_NORMAL_DEFICIT {
        return;
    }
    let mut candidate = solved.clone();
    for _ in 0..MAX_ITERATIONS {
        let mut corrected = false;
        for pair in &pcurves.pairs {
            let Some(surface) = candidate.surfaces.iter_mut().find(|surface| {
                surface.sector_index == pair.sector_index && surface.sector_name == pair.sector_name
            }) else {
                continue;
            };
            let Some(ribbon) = ribbon_for_pair(ribbons, pair) else {
                continue;
            };
            corrected |= refine_pair_samples(&mut surface.surface, pair, ribbon);
        }
        if !corrected {
            break;
        }
    }
    let candidate_deficit = first_pcurve_decline_deficit(&candidate, pcurves, ribbons)
        .unwrap_or_else(|| worst_pcurve_support_normal_deficit(&candidate, pcurves, ribbons));
    let candidate_boundary_residual =
        worst_cap_support_position_residual(&candidate, pcurves, ribbons);
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace exact-boundary/pcurve-strip] candidate deficit {candidate_deficit:.4e} boundary_residual {candidate_boundary_residual:.4e}"
        );
    }
    if candidate_deficit.is_finite()
        && candidate_deficit < baseline_deficit
        && candidate_boundary_residual <= MAX_REFINED_BOUNDARY_RESIDUAL
    {
        *solved = candidate;
    }
}

fn first_pcurve_decline_deficit(
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    ribbons: &PreparedExactBoundaryRibbons,
) -> Option<f64> {
    for pair in &pcurves.pairs {
        let surface = solved.surfaces.iter().find(|surface| {
            surface.sector_index == pair.sector_index && surface.sector_name == pair.sector_name
        })?;
        let ribbon = ribbon_for_pair(ribbons, pair)?;
        let deficit = worst_pair_support_normal_deficit(&surface.surface, pair, ribbon);
        if deficit > MAX_SUPPORT_NORMAL_DEFICIT {
            return Some(deficit);
        }
    }
    None
}

fn worst_pcurve_support_normal_deficit(
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    ribbons: &PreparedExactBoundaryRibbons,
) -> f64 {
    pcurves
        .pairs
        .iter()
        .filter_map(|pair| {
            let surface = solved.surfaces.iter().find(|surface| {
                surface.sector_index == pair.sector_index && surface.sector_name == pair.sector_name
            })?;
            let ribbon = ribbon_for_pair(ribbons, pair)?;
            Some(worst_pair_support_normal_deficit(
                &surface.surface,
                pair,
                ribbon,
            ))
        })
        .fold(0.0, f64::max)
}

pub(super) fn worst_pair_support_normal_deficit(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> f64 {
    ribbon
        .samples
        .iter()
        .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
        .filter_map(|sample| pcurve_support_normal_deficit(surface, pair, sample))
        .fold(0.0, f64::max)
}

fn pcurve_support_normal_deficit(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    sample: &ExactBoundaryRibbonSample,
) -> Option<f64> {
    let cap_uv = pair.cap_pcurve.point_at(
        pair.cap_pcurve.t_start
            + (pair.cap_pcurve.t_end - pair.cap_pcurve.t_start) * sample.parameter,
    );
    support_normal_deficit_at(surface, cap_uv.x, cap_uv.y, sample.support_normal)
}

fn refine_pair_samples(
    surface: &mut NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> bool {
    let mut corrected = refine_pair_least_squares(surface, pair, ribbon);
    corrected |= optimize_pair_controls(surface, pair, ribbon);
    for sample in ribbon
        .samples
        .iter()
        .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
    {
        corrected |= refine_pair_sample(surface, pair, sample);
    }
    corrected
}

fn refine_pair_least_squares(
    surface: &mut NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> bool {
    let targets = collect_strip_targets(surface, pair, ribbon);
    if targets.is_empty() {
        return false;
    }
    let side = targets[0].side;
    let control_indices = interior_rational_strip_basis(surface, side, 0.5)
        .into_iter()
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if control_indices.is_empty() {
        return false;
    }
    let mut normal = vec![vec![0.0; control_indices.len()]; control_indices.len()];
    let mut rhs = vec![Vec3::zero(); control_indices.len()];
    for target in &targets {
        let basis = interior_rational_strip_basis(surface, side, target.parameter);
        let basis_values = control_indices
            .iter()
            .map(|control_index| {
                basis
                    .iter()
                    .find(|(index, _)| index == control_index)
                    .map(|(_, value)| *value)
                    .unwrap_or(0.0)
            })
            .collect::<Vec<_>>();
        let delta = target.target_point - strip_curve(surface, side).point_at(target.parameter);
        for row in 0..control_indices.len() {
            rhs[row] = rhs[row] + delta * (basis_values[row] * target.weight);
            for col in 0..control_indices.len() {
                normal[row][col] += basis_values[row] * basis_values[col] * target.weight;
            }
        }
    }
    for (index, row) in normal.iter_mut().enumerate() {
        row[index] += 1.0e-8;
    }
    let Some(deltas) = solve_vector_system(normal, rhs) else {
        return false;
    };
    let baseline = targets
        .iter()
        .filter_map(|target| {
            support_normal_deficit_at(surface, target.cap_u, target.cap_v, target.support_normal)
        })
        .fold(0.0, f64::max);
    let mut candidate = surface.clone();
    for (control_index, delta) in control_indices.iter().zip(deltas) {
        candidate.control_points[*control_index] = candidate.control_points[*control_index] + delta;
    }
    let candidate_deficit = targets
        .iter()
        .filter_map(|target| {
            support_normal_deficit_at(
                &candidate,
                target.cap_u,
                target.cap_v,
                target.support_normal,
            )
        })
        .fold(0.0, f64::max);
    if candidate_deficit.is_finite() && candidate_deficit < baseline {
        *surface = candidate;
        return true;
    }
    false
}

fn collect_strip_targets(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> Vec<StripTarget> {
    ribbon
        .samples
        .iter()
        .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
        .filter_map(|sample| strip_target(surface, pair, sample))
        .collect()
}

fn strip_target(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    sample: &ExactBoundaryRibbonSample,
) -> Option<StripTarget> {
    let cap_uv = pair.cap_pcurve.point_at(
        pair.cap_pcurve.t_start
            + (pair.cap_pcurve.t_end - pair.cap_pcurve.t_start) * sample.parameter,
    );
    let (side, parameter) = strip_side_from_cap_uv(pair.side, cap_uv.x, cap_uv.y)?;
    let baseline = support_normal_deficit_at(surface, cap_uv.x, cap_uv.y, sample.support_normal)?;
    let current_strip_point = strip_curve(surface, side).point_at(parameter);
    let (target_point, weight) = if baseline <= MAX_SUPPORT_NORMAL_DEFICIT {
        (current_strip_point, 0.5)
    } else {
        let boundary_point = surface.point_at(cap_uv.x, cap_uv.y);
        let step = surface_derivative_control_step(surface, side).max(1.0e-9);
        let target_point = best_local_target_point(
            surface,
            side,
            parameter,
            boundary_point,
            step,
            cap_uv.x,
            cap_uv.y,
            sample.support_normal,
            sample.cross_derivative,
        )?;
        (target_point, 1.0)
    };
    Some(StripTarget {
        side,
        parameter,
        target_point,
        weight,
        cap_u: cap_uv.x,
        cap_v: cap_uv.y,
        support_normal: sample.support_normal,
    })
}

#[allow(clippy::too_many_arguments)]
fn best_local_target_point(
    surface: &NurbsSurface,
    side: StripSide,
    parameter: f64,
    boundary_point: Point3,
    step: f64,
    cap_u: f64,
    cap_v: f64,
    support_normal: Vec3,
    cross_derivative: Vec3,
) -> Option<Point3> {
    let mut best_target = None;
    let mut best_deficit =
        support_normal_deficit_at(surface, cap_u, cap_v, support_normal).unwrap_or(f64::INFINITY);
    for derivative in
        scaled_target_derivatives(surface, side, parameter, support_normal, cross_derivative)
    {
        let mut candidate = surface.clone();
        let target = boundary_point + derivative * step;
        if !correct_strip_point_least_norm(&mut candidate, side, parameter, target) {
            continue;
        }
        let deficit = support_normal_deficit_at(&candidate, cap_u, cap_v, support_normal)
            .unwrap_or(f64::INFINITY);
        if deficit.is_finite() && deficit < best_deficit {
            best_deficit = deficit;
            best_target = Some(target);
        }
    }
    best_target
}

fn refine_pair_sample(
    surface: &mut NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    sample: &ExactBoundaryRibbonSample,
) -> bool {
    let cap_uv = pair.cap_pcurve.point_at(
        pair.cap_pcurve.t_start
            + (pair.cap_pcurve.t_end - pair.cap_pcurve.t_start) * sample.parameter,
    );
    let Some(baseline) =
        support_normal_deficit_at(surface, cap_uv.x, cap_uv.y, sample.support_normal)
    else {
        return false;
    };
    if baseline <= MAX_SUPPORT_NORMAL_DEFICIT {
        return false;
    }
    let Some((side, parameter)) = strip_side_from_cap_uv(pair.side, cap_uv.x, cap_uv.y) else {
        return false;
    };
    let boundary_point = surface.point_at(cap_uv.x, cap_uv.y);
    let step = surface_derivative_control_step(surface, side).max(1.0e-9);
    let targets = scaled_target_derivatives(
        surface,
        side,
        parameter,
        sample.support_normal,
        sample.cross_derivative,
    );
    apply_best_pcurve_strip_target(
        surface,
        PcurveStripCorrection {
            side,
            parameter,
            boundary_point,
            step,
            cap_u: cap_uv.x,
            cap_v: cap_uv.y,
            support_normal: sample.support_normal,
        },
        &targets,
    )
}

struct PcurveStripCorrection {
    side: StripSide,
    parameter: f64,
    boundary_point: Point3,
    step: f64,
    cap_u: f64,
    cap_v: f64,
    support_normal: Vec3,
}

struct StripTarget {
    side: StripSide,
    parameter: f64,
    target_point: Point3,
    weight: f64,
    cap_u: f64,
    cap_v: f64,
    support_normal: Vec3,
}

fn apply_best_pcurve_strip_target(
    surface: &mut NurbsSurface,
    correction: PcurveStripCorrection,
    target_derivatives: &[Vec3],
) -> bool {
    let baseline = support_normal_deficit_at(
        surface,
        correction.cap_u,
        correction.cap_v,
        correction.support_normal,
    )
    .unwrap_or(f64::INFINITY);
    let mut best_surface = None;
    let mut best_deficit = baseline;
    for derivative in target_derivatives {
        let mut candidate = surface.clone();
        let target = correction.boundary_point + *derivative * correction.step;
        if !correct_strip_point_least_norm(
            &mut candidate,
            correction.side,
            correction.parameter,
            target,
        ) {
            continue;
        }
        let deficit = support_normal_deficit_at(
            &candidate,
            correction.cap_u,
            correction.cap_v,
            correction.support_normal,
        )
        .unwrap_or(f64::INFINITY);
        if deficit.is_finite() && deficit < best_deficit {
            best_deficit = deficit;
            best_surface = Some(candidate);
        }
    }
    if let Some(best_surface) = best_surface {
        *surface = best_surface;
        return true;
    }
    false
}

#[allow(clippy::needless_range_loop)]
fn solve_vector_system(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<Vec3>) -> Option<Vec<Vec3>> {
    let n = rhs.len();
    for pivot in 0..n {
        let pivot_row = (pivot..n).max_by(|a, b| {
            matrix[*a][pivot]
                .abs()
                .partial_cmp(&matrix[*b][pivot].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if matrix[pivot_row][pivot].abs() <= 1.0e-12 {
            return None;
        }
        matrix.swap(pivot, pivot_row);
        rhs.swap(pivot, pivot_row);
        let pivot_value = matrix[pivot][pivot];
        for col in pivot..n {
            matrix[pivot][col] /= pivot_value;
        }
        rhs[pivot] = rhs[pivot] * (1.0 / pivot_value);
        for row in 0..n {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            if factor.abs() <= 1.0e-15 {
                continue;
            }
            for col in pivot..n {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] = rhs[row] - rhs[pivot] * factor;
        }
    }
    Some(rhs)
}

fn ribbon_for_pair<'a>(
    ribbons: &'a PreparedExactBoundaryRibbons,
    pair: &PreparedExactBoundaryPcurvePair,
) -> Option<&'a PreparedExactBoundaryRibbon> {
    ribbons.ribbons.iter().find(|ribbon| {
        ribbon.sector_index == pair.sector_index
            && ribbon.side == pair.side
            && ribbon.edge == pair.edge
            && ribbon.support_half_edge == pair.support_half_edge
    })
}

pub(super) fn strip_side_from_cap_uv(
    side: ExactBoundaryOuterSide,
    u: f64,
    v: f64,
) -> Option<(StripSide, f64)> {
    match side {
        ExactBoundaryOuterSide::FilletEnd => Some((StripSide::FilletEnd, u)),
        ExactBoundaryOuterSide::FaceEdge => Some((StripSide::FaceEdge, v)),
    }
}

fn surface_derivative_control_step(surface: &NurbsSurface, side: StripSide) -> f64 {
    let (degree, knots) = match side {
        StripSide::FilletEnd => (surface.degree_v, surface.knots_v.as_slice()),
        StripSide::FaceEdge => (surface.degree_u, surface.knots_u.as_slice()),
    };
    if degree == 0 {
        return 0.0;
    }
    let lower = knots.get(1).copied().unwrap_or(0.0);
    let upper = knots.get(degree + 1).copied().unwrap_or(1.0);
    ((upper - lower) / degree as f64).abs()
}
