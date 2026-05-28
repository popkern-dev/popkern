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
    let control_indices = interior_rational_strip_basis(surface, side, 0.5)
        .into_iter()
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
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
    corrected
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
    compact_directions(directions)
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
