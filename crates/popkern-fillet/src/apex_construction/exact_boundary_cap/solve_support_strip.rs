//! Local support-normal strip correction for exact-boundary surface solves.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use super::super::corner::PreparedExactBoundarySectorBoundary;
use super::super::ribbon::{
    ExactBoundaryRibbonSample, PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons,
};
use super::super::sector::ExactBoundaryOuterSide;
use super::solve_derivatives::{
    BoundaryOrientation, boundary_derivative_control_step, ribbon_derivative,
};

pub(super) const MAX_SUPPORT_NORMAL_DEFICIT: f64 = 5.0e-2;
const MIN_NORMAL_LENGTH: f64 = 1.0e-12;
const INTERIOR_MARGIN: f64 = 1.0e-5;
pub(super) const MAX_ITERATIONS: usize = 12;
const CORRECTION_DAMPING: f64 = 1.0;
const MAX_CONTROL_CORRECTION_LENGTH: f64 = 0.5;

pub(super) fn apply_support_normal_strip_correction(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbons: &PreparedExactBoundaryRibbons,
    orientation: BoundaryOrientation,
    surface: &mut NurbsSurface,
) {
    let baseline_deficit = worst_support_normal_deficit(sector, ribbons, orientation, surface);
    if baseline_deficit <= MAX_SUPPORT_NORMAL_DEFICIT {
        return;
    }
    let mut candidate = surface.clone();
    for _ in 0..MAX_ITERATIONS {
        let mut corrected = false;
        for ribbon in ribbons
            .ribbons
            .iter()
            .filter(|ribbon| ribbon.sector_index == sector.sector_index)
        {
            corrected |= correct_ribbon_samples(sector, ribbon, orientation, &mut candidate);
        }
        if !corrected {
            break;
        }
    }
    let candidate_deficit = worst_support_normal_deficit(sector, ribbons, orientation, &candidate);
    if candidate_deficit.is_finite() && candidate_deficit < baseline_deficit {
        *surface = candidate;
    }
}

fn worst_support_normal_deficit(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbons: &PreparedExactBoundaryRibbons,
    orientation: BoundaryOrientation,
    surface: &NurbsSurface,
) -> f64 {
    ribbons
        .ribbons
        .iter()
        .filter(|ribbon| ribbon.sector_index == sector.sector_index)
        .flat_map(|ribbon| {
            ribbon
                .samples
                .iter()
                .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
                .filter_map(move |sample| {
                    support_normal_deficit(surface, ribbon.side, orientation, sample)
                })
        })
        .fold(0.0, f64::max)
}

fn correct_ribbon_samples(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbon: &PreparedExactBoundaryRibbon,
    orientation: BoundaryOrientation,
    surface: &mut NurbsSurface,
) -> bool {
    let mut corrected = false;
    for sample in ribbon
        .samples
        .iter()
        .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
    {
        corrected |= correct_sample_if_needed(sector, ribbon, orientation, sample, surface);
    }
    corrected
}

fn correct_sample_if_needed(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbon: &PreparedExactBoundaryRibbon,
    orientation: BoundaryOrientation,
    sample: &ExactBoundaryRibbonSample,
    surface: &mut NurbsSurface,
) -> bool {
    let Some(deficit) = support_normal_deficit(surface, ribbon.side, orientation, sample) else {
        return false;
    };
    if deficit <= MAX_SUPPORT_NORMAL_DEFICIT {
        return false;
    }
    match ribbon.side {
        ExactBoundaryOuterSide::FilletEnd => {
            correct_fillet_end_strip(sector, ribbon, orientation, sample, surface)
        }
        ExactBoundaryOuterSide::FaceEdge => {
            correct_face_edge_strip(sector, ribbon, orientation, sample, surface)
        }
    }
}

fn correct_fillet_end_strip(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbon: &PreparedExactBoundaryRibbon,
    orientation: BoundaryOrientation,
    sample: &ExactBoundaryRibbonSample,
    surface: &mut NurbsSurface,
) -> bool {
    let parameter = oriented_parameter(sample.parameter, orientation.fillet_reversed);
    let ribbon_derivative = ribbon_derivative(ribbon, parameter, orientation.fillet_reversed);
    let step = boundary_derivative_control_step(&sector.face_boundary.curve);
    let boundary_point = surface.point_at(parameter, 0.0);
    let targets = target_derivatives(
        surface,
        StripSide::FilletEnd,
        parameter,
        sample.support_normal,
        ribbon_derivative,
    );
    apply_best_strip_target(
        surface,
        StripCorrection {
            side: StripSide::FilletEnd,
            outer_side: ribbon.side,
            orientation,
            parameter,
            boundary_point,
            step,
            sample,
        },
        &targets,
    )
}

fn correct_face_edge_strip(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbon: &PreparedExactBoundaryRibbon,
    orientation: BoundaryOrientation,
    sample: &ExactBoundaryRibbonSample,
    surface: &mut NurbsSurface,
) -> bool {
    let parameter = oriented_parameter(sample.parameter, orientation.face_reversed);
    let ribbon_derivative = ribbon_derivative(ribbon, parameter, orientation.face_reversed);
    let step = boundary_derivative_control_step(&sector.fillet_boundary.curve);
    let boundary_point = surface.point_at(0.0, parameter);
    let targets = target_derivatives(
        surface,
        StripSide::FaceEdge,
        parameter,
        sample.support_normal,
        ribbon_derivative,
    );
    apply_best_strip_target(
        surface,
        StripCorrection {
            side: StripSide::FaceEdge,
            outer_side: ribbon.side,
            orientation,
            parameter,
            boundary_point,
            step,
            sample,
        },
        &targets,
    )
}

struct StripCorrection<'a> {
    side: StripSide,
    outer_side: ExactBoundaryOuterSide,
    orientation: BoundaryOrientation,
    parameter: f64,
    boundary_point: Point3,
    step: f64,
    sample: &'a ExactBoundaryRibbonSample,
}

fn apply_best_strip_target(
    surface: &mut NurbsSurface,
    correction: StripCorrection<'_>,
    target_derivatives: &[Vec3],
) -> bool {
    let baseline = support_normal_deficit(
        surface,
        correction.outer_side,
        correction.orientation,
        correction.sample,
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
        let deficit = support_normal_deficit(
            &candidate,
            correction.outer_side,
            correction.orientation,
            correction.sample,
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

fn target_derivatives(
    surface: &NurbsSurface,
    side: StripSide,
    parameter: f64,
    support_normal: Vec3,
    ribbon_derivative: Vec3,
) -> Vec<Vec3> {
    let scale = ribbon_derivative.length().max(1.0e-9);
    let mut derivatives = vec![ribbon_derivative, -ribbon_derivative];
    if let Some(normal_derivative) =
        normal_constrained_derivative(surface, side, parameter, support_normal, scale)
    {
        derivatives.push(normal_derivative);
        derivatives.push(-normal_derivative);
    }
    derivatives
}

pub(super) fn scaled_target_derivatives(
    surface: &NurbsSurface,
    side: StripSide,
    parameter: f64,
    support_normal: Vec3,
    ribbon_derivative: Vec3,
) -> Vec<Vec3> {
    target_derivatives(surface, side, parameter, support_normal, ribbon_derivative)
        .into_iter()
        .flat_map(|derivative| {
            [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
                .into_iter()
                .map(move |scale| derivative * scale)
        })
        .collect()
}

fn normal_constrained_derivative(
    surface: &NurbsSurface,
    side: StripSide,
    parameter: f64,
    support_normal: Vec3,
    scale: f64,
) -> Option<Vec3> {
    let support_normal = support_normal.normalized()?;
    let derivative = match side {
        StripSide::FilletEnd => {
            let tangent = surface.partial_u(parameter, 0.0).normalized()?;
            support_normal.cross(&tangent)
        }
        StripSide::FaceEdge => {
            let tangent = surface.partial_v(0.0, parameter).normalized()?;
            tangent.cross(&support_normal)
        }
    };
    derivative.normalized().map(|unit| unit * scale)
}

pub(super) fn correct_strip_point_least_norm(
    surface: &mut NurbsSurface,
    side: StripSide,
    parameter: f64,
    target: Point3,
) -> bool {
    let basis = interior_rational_strip_basis(surface, side, parameter);
    let basis_norm = basis.iter().map(|(_, value)| value * value).sum::<f64>();
    if basis_norm <= 1.0e-12 {
        return false;
    }
    let current = strip_curve(surface, side).point_at(parameter);
    let delta = (target - current) * CORRECTION_DAMPING;
    for (control_index, basis_value) in basis {
        let control_delta = limit_vector(delta * (basis_value / basis_norm));
        surface.control_points[control_index] =
            surface.control_points[control_index] + control_delta;
    }
    true
}

pub(super) fn interior_rational_strip_basis(
    surface: &NurbsSurface,
    side: StripSide,
    parameter: f64,
) -> Vec<(usize, f64)> {
    (1..strip_control_count(surface, side).saturating_sub(1))
        .filter_map(|local_index| {
            let value = rational_strip_basis(surface, side, local_index, parameter);
            (value.is_finite() && value.abs() > 1.0e-9)
                .then_some((strip_surface_index(surface, side, local_index), value))
        })
        .collect()
}

fn rational_strip_basis(
    surface: &NurbsSurface,
    side: StripSide,
    local_index: usize,
    parameter: f64,
) -> f64 {
    let controls = (0..strip_control_count(surface, side))
        .map(|index| {
            if index == local_index {
                Point3::new(1.0, 0.0, 0.0)
            } else {
                Point3::new(0.0, 0.0, 0.0)
            }
        })
        .collect();
    let basis_curve = NurbsCurve {
        degree: strip_degree(surface, side),
        knots: strip_knots(surface, side).to_vec(),
        control_points: controls,
        weights: strip_weights(surface, side),
    };
    basis_curve.point_at(parameter).x
}

pub(super) fn strip_curve(surface: &NurbsSurface, side: StripSide) -> NurbsCurve {
    NurbsCurve {
        degree: strip_degree(surface, side),
        knots: strip_knots(surface, side).to_vec(),
        control_points: (0..strip_control_count(surface, side))
            .map(|local_index| {
                surface.control_points[strip_surface_index(surface, side, local_index)]
            })
            .collect(),
        weights: strip_weights(surface, side),
    }
}

fn support_normal_deficit(
    surface: &NurbsSurface,
    side: ExactBoundaryOuterSide,
    orientation: BoundaryOrientation,
    sample: &ExactBoundaryRibbonSample,
) -> Option<f64> {
    let (u, v) = match side {
        ExactBoundaryOuterSide::FilletEnd => (
            oriented_parameter(sample.parameter, orientation.fillet_reversed),
            0.0,
        ),
        ExactBoundaryOuterSide::FaceEdge => (
            0.0,
            oriented_parameter(sample.parameter, orientation.face_reversed),
        ),
    };
    support_normal_deficit_at(surface, u, v, sample.support_normal)
}

pub(super) fn support_normal_deficit_at(
    surface: &NurbsSurface,
    u: f64,
    v: f64,
    support_normal: Vec3,
) -> Option<f64> {
    let (u, v) = interior_normal_sample_parameters(surface, u, v);
    let cap_normal = surface.partial_u(u, v).cross(&surface.partial_v(u, v));
    let cap_length = cap_normal.length();
    let support_length = support_normal.length();
    if cap_length < MIN_NORMAL_LENGTH
        || support_length < MIN_NORMAL_LENGTH
        || !cap_length.is_finite()
        || !support_length.is_finite()
    {
        return Some(f64::INFINITY);
    }
    Some(1.0 - (cap_normal.dot(&support_normal) / (cap_length * support_length)).abs())
}

fn interior_normal_sample_parameters(surface: &NurbsSurface, u: f64, v: f64) -> (f64, f64) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    (
        interior_parameter(u, u_min, u_max),
        interior_parameter(v, v_min, v_max),
    )
}

fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * INTERIOR_MARGIN).max(INTERIOR_MARGIN);
    value.clamp(min + margin, max - margin)
}

fn oriented_parameter(parameter: f64, reversed: bool) -> f64 {
    if reversed { 1.0 - parameter } else { parameter }
}

fn limit_vector(vector: Vec3) -> Vec3 {
    let length = vector.length();
    if length <= MAX_CONTROL_CORRECTION_LENGTH || length <= f64::EPSILON {
        vector
    } else {
        vector * (MAX_CONTROL_CORRECTION_LENGTH / length)
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) enum StripSide {
    FilletEnd,
    FaceEdge,
}

fn strip_control_count(surface: &NurbsSurface, side: StripSide) -> usize {
    match side {
        StripSide::FilletEnd => surface.num_u,
        StripSide::FaceEdge => surface.num_v,
    }
}

fn strip_degree(surface: &NurbsSurface, side: StripSide) -> usize {
    match side {
        StripSide::FilletEnd => surface.degree_u,
        StripSide::FaceEdge => surface.degree_v,
    }
}

fn strip_knots(surface: &NurbsSurface, side: StripSide) -> &[f64] {
    match side {
        StripSide::FilletEnd => &surface.knots_u,
        StripSide::FaceEdge => &surface.knots_v,
    }
}

fn strip_weights(surface: &NurbsSurface, side: StripSide) -> Vec<f64> {
    (0..strip_control_count(surface, side))
        .map(|local_index| surface.weights[strip_surface_index(surface, side, local_index)])
        .collect()
}

fn strip_surface_index(surface: &NurbsSurface, side: StripSide, local_index: usize) -> usize {
    match side {
        StripSide::FilletEnd => local_index * surface.num_v + 1,
        StripSide::FaceEdge => surface.num_v + local_index,
    }
}
