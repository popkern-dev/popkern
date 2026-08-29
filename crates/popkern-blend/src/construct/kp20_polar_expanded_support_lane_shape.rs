//! Support-lane shape primitives for expanded KP20 topology.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

use super::kp20_polar_bi6_split::split_bi6_patch_u;
use super::kp20_polar_expanded_support_topology::Kp20ExpandedSupportMode;
use super::kp20_polar_fold_summary::Kp20FoldSite;
use super::kp20_polar_support_surface_split::evaluate_bi6_normal;

pub(super) fn reauthored_support_lanes(
    base: &Kp20Bi6Patch,
    target: Kp20FoldSite,
    split_parameter: f64,
    mode: Kp20ExpandedSupportMode,
    scale: f64,
) -> (Kp20Bi6Patch, Kp20Bi6Patch) {
    let (mut boundary_lane, mut interior_lane) = split_bi6_patch_u(base, split_parameter);
    let normal = evaluate_bi6_normal(base, target.u, target.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let span = support_span(base);
    for column in 0..7 {
        let v = column as f64 / 6.0;
        let local = (1.0 - ((v - target.v).abs() / 0.42)).max(0.0);
        let direction = expanded_support_direction(base, normal, mode, v);
        let sign = if matches!(mode, Kp20ExpandedSupportMode::InterfaceContract) {
            -1.0
        } else {
            1.0
        };
        let relief = direction * (sign * span * scale * 0.10 * local * local);
        boundary_lane.cps[6][column] = boundary_lane.cps[6][column] + relief;
        interior_lane.cps[0][column] = boundary_lane.cps[6][column];
        boundary_lane.cps[4][column] = boundary_lane.cps[4][column] + relief * 0.30;
        boundary_lane.cps[5][column] = boundary_lane.cps[5][column] + relief * 0.65;
        interior_lane.cps[1][column] = interior_lane.cps[1][column] + relief * 0.65;
        interior_lane.cps[2][column] = interior_lane.cps[2][column] + relief * 0.30;
        if matches!(mode, Kp20ExpandedSupportMode::CoupledLaneSpread) {
            spread_adjacent_lane_rows(&mut boundary_lane, &mut interior_lane, column, relief);
        }
    }
    (boundary_lane, interior_lane)
}

fn spread_adjacent_lane_rows(
    boundary_lane: &mut Kp20Bi6Patch,
    interior_lane: &mut Kp20Bi6Patch,
    column: usize,
    relief: Vec3,
) {
    boundary_lane.cps[2][column] = boundary_lane.cps[2][column] - relief * 0.18;
    boundary_lane.cps[3][column] = boundary_lane.cps[3][column] - relief * 0.10;
    interior_lane.cps[3][column] = interior_lane.cps[3][column] + relief * 0.22;
    interior_lane.cps[4][column] = interior_lane.cps[4][column] + relief * 0.14;
}

fn expanded_support_direction(
    base: &Kp20Bi6Patch,
    normal: Vec3,
    mode: Kp20ExpandedSupportMode,
    v: f64,
) -> Vec3 {
    match mode {
        Kp20ExpandedSupportMode::SupportNormalRelief => normal,
        Kp20ExpandedSupportMode::ChordRelief
        | Kp20ExpandedSupportMode::InterfaceContract
        | Kp20ExpandedSupportMode::CoupledLaneSpread => {
            let column = ((v * 6.0).round() as usize).clamp(0, 6);
            let chord = base.cps[1][6] - base.cps[1][0];
            let depth = base.cps[6][column] - base.cps[1][column];
            (chord + depth).normalized().unwrap_or(normal)
        }
    }
}

fn support_span(surface: &Kp20Bi6Patch) -> f64 {
    let support_width = surface.cps[1][0].distance_to(&surface.cps[1][6]);
    let interface_width = surface.cps[6][0].distance_to(&surface.cps[6][6]);
    let cross_depth = surface.cps[1][3].distance_to(&surface.cps[6][3]);
    ((support_width + interface_width + cross_depth) / 3.0).max(1.0e-9)
}
