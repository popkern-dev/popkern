//! Support-ribbon derivative rows for exact-boundary surface solves.

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::super::corner::PreparedExactBoundarySectorBoundary;
use super::super::ribbon::{PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons};
use super::super::sector::ExactBoundaryOuterSide;

#[derive(Debug, Clone, Copy)]
pub(super) struct BoundaryOrientation {
    pub(super) fillet_reversed: bool,
    pub(super) face_reversed: bool,
}

pub(super) fn apply_derivative_rows(
    sector: &PreparedExactBoundarySectorBoundary,
    ribbons: &PreparedExactBoundaryRibbons,
    orientation: BoundaryOrientation,
    controls: &mut [Point3],
    num_u: usize,
    num_v: usize,
) {
    let mut derivative_controls = vec![None; controls.len()];
    if let Some(ribbon) = ribbon_for(
        ribbons,
        sector.sector_index,
        ExactBoundaryOuterSide::FilletEnd,
    ) {
        for i in 1..num_u.saturating_sub(1) {
            let parameter = i as f64 / (num_u - 1) as f64;
            let derivative = ribbon_derivative(ribbon, parameter, orientation.fillet_reversed);
            let base = controls[index(i, 0, num_v)];
            merge_derivative_control(
                &mut derivative_controls,
                index(i, 1, num_v),
                base + derivative * boundary_derivative_control_step(&sector.face_boundary.curve),
            );
        }
    }
    if let Some(ribbon) = ribbon_for(
        ribbons,
        sector.sector_index,
        ExactBoundaryOuterSide::FaceEdge,
    ) {
        for j in 1..num_v.saturating_sub(1) {
            let parameter = j as f64 / (num_v - 1) as f64;
            let derivative = ribbon_derivative(ribbon, parameter, orientation.face_reversed);
            let base = controls[index(0, j, num_v)];
            merge_derivative_control(
                &mut derivative_controls,
                index(1, j, num_v),
                base + derivative * boundary_derivative_control_step(&sector.fillet_boundary.curve),
            );
        }
    }
    for (control, derivative_control) in controls.iter_mut().zip(derivative_controls) {
        if let Some(derivative_control) = derivative_control {
            *control = derivative_control;
        }
    }
}

fn merge_derivative_control(controls: &mut [Option<Point3>], index: usize, candidate: Point3) {
    controls[index] = Some(match controls[index] {
        Some(existing) => existing + (candidate - existing) * 0.5,
        None => candidate,
    });
}

fn ribbon_for(
    ribbons: &PreparedExactBoundaryRibbons,
    sector_index: usize,
    side: ExactBoundaryOuterSide,
) -> Option<&PreparedExactBoundaryRibbon> {
    ribbons
        .ribbons
        .iter()
        .find(|ribbon| ribbon.sector_index == sector_index && ribbon.side == side)
}

fn ribbon_derivative(
    ribbon: &PreparedExactBoundaryRibbon,
    parameter: f64,
    boundary_reversed: bool,
) -> Vec3 {
    let ribbon_parameter = if boundary_reversed {
        1.0 - parameter
    } else {
        parameter
    };
    ribbon
        .samples
        .windows(2)
        .find_map(|window| {
            let left = window[0];
            let right = window[1];
            if ribbon_parameter < left.parameter || ribbon_parameter > right.parameter {
                return None;
            }
            let span = right.parameter - left.parameter;
            if span.abs() <= f64::EPSILON {
                return Some(left.cross_derivative);
            }
            let fraction = ((ribbon_parameter - left.parameter) / span).clamp(0.0, 1.0);
            Some(
                left.cross_derivative + (right.cross_derivative - left.cross_derivative) * fraction,
            )
        })
        .or_else(|| {
            ribbon
                .samples
                .iter()
                .min_by(|a, b| {
                    (a.parameter - ribbon_parameter)
                        .abs()
                        .partial_cmp(&(b.parameter - ribbon_parameter).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|sample| sample.cross_derivative)
        })
        .unwrap_or(Vec3::zero())
}

fn boundary_derivative_control_step(boundary: &NurbsCurve) -> f64 {
    if boundary.degree == 0 {
        return 0.0;
    }
    let lower = boundary.knots.get(1).copied().unwrap_or(0.0);
    let upper = boundary
        .knots
        .get(boundary.degree + 1)
        .copied()
        .unwrap_or(1.0);
    ((upper - lower) / boundary.degree as f64).abs()
}

fn index(i: usize, j: usize, num_v: usize) -> usize {
    i * num_v + j
}
