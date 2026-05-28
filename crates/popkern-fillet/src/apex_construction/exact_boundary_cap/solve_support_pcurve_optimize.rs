//! Bounded coordinate optimization for post-pcurve support-normal strips.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsSurface, Vec3};

use super::super::pcurves::PreparedExactBoundaryPcurvePair;
use super::super::ribbon::PreparedExactBoundaryRibbon;
use super::solve_support_pcurve_strip::{
    strip_side_from_cap_uv, worst_pair_support_normal_deficit,
};
use super::solve_support_strip::{
    MAX_SUPPORT_NORMAL_DEFICIT, StripSide, interior_rational_strip_basis, support_normal_deficit_at,
};

pub(super) fn optimize_pair_controls(
    surface: &mut NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> bool {
    let Some((side, _)) = first_pair_strip_side(surface, pair, ribbon) else {
        return false;
    };
    let mut corrected = false;
    let mut best_deficit = worst_pair_support_normal_deficit(surface, pair, ribbon);
    if best_deficit <= MAX_SUPPORT_NORMAL_DEFICIT {
        return false;
    }
    let control_indices = strip_band_control_indices(surface, side);
    let directions = pair_search_directions(surface, pair, ribbon, side);
    if control_indices.is_empty() || directions.is_empty() {
        return false;
    }
    for step in [0.25, 0.125, 0.0625, 0.03125, 0.015625] {
        let mut improved_at_step = true;
        while improved_at_step {
            improved_at_step = false;
            for control_index in &control_indices {
                for direction in &directions {
                    for sign in [-1.0, 1.0] {
                        let mut candidate = surface.clone();
                        candidate.control_points[*control_index] =
                            candidate.control_points[*control_index] + *direction * (step * sign);
                        let deficit = worst_pair_support_normal_deficit(&candidate, pair, ribbon);
                        if deficit.is_finite() && deficit < best_deficit {
                            *surface = candidate;
                            best_deficit = deficit;
                            corrected = true;
                            improved_at_step = true;
                        }
                    }
                }
            }
        }
    }
    corrected |= optimize_pair_gradient(surface, pair, ribbon, &control_indices, best_deficit);
    corrected
}

fn optimize_pair_gradient(
    surface: &mut NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
    control_indices: &[usize],
    initial_deficit: f64,
) -> bool {
    let mut best_deficit = initial_deficit;
    let mut corrected = false;
    for _ in 0..20 {
        let gradient = objective_gradient(surface, pair, ribbon, control_indices, 1.0e-4);
        let gradient_norm = gradient
            .iter()
            .map(Vec3::length_squared)
            .sum::<f64>()
            .sqrt();
        if gradient_norm <= 1.0e-12 || !gradient_norm.is_finite() {
            break;
        }
        let mut improved = false;
        for scale in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625] {
            let mut candidate = surface.clone();
            for (control_index, gradient) in control_indices.iter().zip(&gradient) {
                candidate.control_points[*control_index] =
                    candidate.control_points[*control_index] - *gradient * (scale / gradient_norm);
            }
            let deficit = worst_pair_support_normal_deficit(&candidate, pair, ribbon);
            if deficit.is_finite() && deficit < best_deficit {
                *surface = candidate;
                best_deficit = deficit;
                corrected = true;
                improved = true;
                break;
            }
        }
        if !improved {
            break;
        }
    }
    corrected
}

fn objective_gradient(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
    control_indices: &[usize],
    epsilon: f64,
) -> Vec<Vec3> {
    control_indices
        .iter()
        .map(|control_index| {
            Vec3::new(
                partial_objective(surface, pair, ribbon, *control_index, Axis::X, epsilon),
                partial_objective(surface, pair, ribbon, *control_index, Axis::Y, epsilon),
                partial_objective(surface, pair, ribbon, *control_index, Axis::Z, epsilon),
            )
        })
        .collect()
}

fn partial_objective(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
    control_index: usize,
    axis: Axis,
    epsilon: f64,
) -> f64 {
    let mut plus = surface.clone();
    adjust_control(&mut plus, control_index, axis, epsilon);
    let plus_objective = worst_pair_support_normal_deficit(&plus, pair, ribbon);
    let mut minus = surface.clone();
    adjust_control(&mut minus, control_index, axis, -epsilon);
    let minus_objective = worst_pair_support_normal_deficit(&minus, pair, ribbon);
    (plus_objective - minus_objective) / (2.0 * epsilon)
}

fn adjust_control(surface: &mut NurbsSurface, control_index: usize, axis: Axis, delta: f64) {
    match axis {
        Axis::X => surface.control_points[control_index].x += delta,
        Axis::Y => surface.control_points[control_index].y += delta,
        Axis::Z => surface.control_points[control_index].z += delta,
    }
}

#[derive(Debug, Clone, Copy)]
enum Axis {
    X,
    Y,
    Z,
}

fn first_pair_strip_side(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> Option<(StripSide, f64)> {
    ribbon
        .samples
        .iter()
        .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
        .find_map(|sample| {
            let cap_uv = pair.cap_pcurve.point_at(
                pair.cap_pcurve.t_start
                    + (pair.cap_pcurve.t_end - pair.cap_pcurve.t_start) * sample.parameter,
            );
            let deficit =
                support_normal_deficit_at(surface, cap_uv.x, cap_uv.y, sample.support_normal)?;
            (deficit > MAX_SUPPORT_NORMAL_DEFICIT)
                .then(|| strip_side_from_cap_uv(pair.side, cap_uv.x, cap_uv.y))
                .flatten()
        })
}

fn pair_search_directions(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
    side: StripSide,
) -> Vec<Vec3> {
    let mut directions = Vec::new();
    for sample in ribbon
        .samples
        .iter()
        .filter(|sample| sample.parameter > 0.0 && sample.parameter < 1.0)
    {
        let cap_uv = pair.cap_pcurve.point_at(
            pair.cap_pcurve.t_start
                + (pair.cap_pcurve.t_end - pair.cap_pcurve.t_start) * sample.parameter,
        );
        let Some((_, parameter)) = strip_side_from_cap_uv(pair.side, cap_uv.x, cap_uv.y) else {
            continue;
        };
        if let Some(unit) = sample.cross_derivative.normalized() {
            directions.push(unit);
        }
        if let Some(unit) = sample.support_normal.normalized() {
            directions.push(unit);
        }
        let tangent = match side {
            StripSide::FilletEnd => surface.partial_u(parameter, 0.0),
            StripSide::FaceEdge => surface.partial_v(0.0, parameter),
        };
        if let Some(unit) = tangent.normalized() {
            directions.push(unit);
        }
    }
    directions.extend([
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]);
    compact_directions(directions)
}

fn strip_band_control_indices(surface: &NurbsSurface, side: StripSide) -> Vec<usize> {
    let mut indices = Vec::new();
    for layer in [1, 2] {
        indices.extend(
            interior_rational_strip_basis(surface, side, 0.5)
                .into_iter()
                .filter_map(|(index, _)| strip_band_index(surface, side, index, layer)),
        );
    }
    indices
}

fn strip_band_index(
    surface: &NurbsSurface,
    side: StripSide,
    first_strip_index: usize,
    layer: usize,
) -> Option<usize> {
    match side {
        StripSide::FilletEnd => {
            let local = first_strip_index / surface.num_v;
            (local > 0 && local + 1 < surface.num_u && layer < surface.num_v)
                .then_some(local * surface.num_v + layer)
        }
        StripSide::FaceEdge => {
            let local = first_strip_index.saturating_sub(surface.num_v);
            (local > 0 && local + 1 < surface.num_v && layer < surface.num_u)
                .then_some(layer * surface.num_v + local)
        }
    }
}

fn compact_directions(directions: Vec<Vec3>) -> Vec<Vec3> {
    let mut compact = Vec::new();
    for direction in directions {
        if compact
            .iter()
            .all(|existing: &Vec3| existing.dot(&direction).abs() < 0.98)
        {
            compact.push(direction);
        }
    }
    compact
}
