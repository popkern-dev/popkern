//! Candidate-only KP20 inner-pole patch prototype.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

use super::kp20_polar_multipatch_interface::{
    Kp20MultiPatchInterfaceBuildReport, Kp20MultiPatchInterfaceCurve,
};
use super::kp20_polar_multipatch_inventory::Kp20MultiPatchInventoryReport;

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchInnerBuildReport {
    pub(super) attempted_inner_patch_count: usize,
    pub(super) built_inner_patch_count: usize,
    pub(super) finite_control_point_count: usize,
    pub(super) invalid_control_point_count: usize,
    pub(super) boundary_constraint_count: usize,
    pub(super) adjacent_row_constraint_count: usize,
    pub(super) interior_constraint_count: usize,
    pub(super) exact_boundary_owned_endpoint_count: usize,
    pub(super) max_boundary_g0_residual: f64,
    pub(super) max_adjacent_row_residual: f64,
    pub(super) max_interior_residual: f64,
    pub(super) min_adjacent_row_step: f64,
    pub(super) max_coupled_core_displacement: f64,
    pub(super) decline_reason: Option<&'static str>,
    pub(super) inner_patch_names: Vec<String>,
    pub(super) inner_patches: Vec<Kp20MultiPatchInnerPatch>,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchInnerPatch {
    pub(super) patch_name: String,
    pub(super) sector_index: usize,
    pub(super) boundary_interface_names: Vec<String>,
    pub(super) surface: Kp20Bi6Patch,
}

pub(super) fn build_multipatch_inner_patches(
    inventory: &Kp20MultiPatchInventoryReport,
    interfaces: &Kp20MultiPatchInterfaceBuildReport,
) -> Kp20MultiPatchInnerBuildReport {
    let mut inner_patches = Vec::new();
    let mut finite_control_point_count = 0_usize;
    let mut invalid_control_point_count = 0_usize;
    let mut boundary_constraint_count = 0_usize;
    let mut adjacent_row_constraint_count = 0_usize;
    let mut interior_constraint_count = 0_usize;
    let mut exact_boundary_owned_endpoint_count = 0_usize;
    let mut max_boundary_g0_residual = 0.0_f64;
    let mut max_adjacent_row_residual = 0.0_f64;
    let mut max_interior_residual = 0.0_f64;
    let mut min_adjacent_row_step = f64::INFINITY;
    let mut max_coupled_core_displacement = 0.0_f64;

    for sector in &inventory.sector_reports {
        let sector_interfaces = sector
            .interface_names
            .iter()
            .filter_map(|name| interface_by_name(interfaces, name))
            .collect::<Vec<_>>();
        if sector_interfaces.is_empty() {
            continue;
        }

        let centroid = interface_centroid(&sector_interfaces);
        let mut surface = Kp20Bi6Patch {
            cps: [[centroid; 7]; 7],
        };
        for (boundary_index, interface) in sector_interfaces.iter().enumerate() {
            assign_boundary_controls(&mut surface, boundary_index, interface.controls);
            boundary_constraint_count += interface.controls.len();
            assign_boundary_adjacent_controls(
                &mut surface,
                boundary_index,
                interface.controls,
                interface.support_adjacent_controls,
                sector_interfaces.len(),
            );
            adjacent_row_constraint_count += ADJACENT_CONTROL_COUNT;
            exact_boundary_owned_endpoint_count +=
                exact_boundary_owned_endpoint_count_for(boundary_index, sector_interfaces.len());
        }
        let core_report = assign_inner_core_controls(&mut surface);
        interior_constraint_count += INNER_CORE_CONTROL_COUNT;
        for (boundary_index, interface) in sector_interfaces.iter().enumerate() {
            max_boundary_g0_residual = max_boundary_g0_residual.max(boundary_residual(
                &surface,
                boundary_index,
                interface,
            ));
            let adjacent = adjacent_row_residual(
                &surface,
                boundary_index,
                interface.controls,
                interface.support_adjacent_controls,
                sector_interfaces.len(),
            );
            max_adjacent_row_residual =
                max_adjacent_row_residual.max(adjacent.max_construction_residual);
            min_adjacent_row_step = min_adjacent_row_step.min(adjacent.min_step);
        }
        max_interior_residual = max_interior_residual.max(inner_core_residual(&surface));
        max_coupled_core_displacement =
            max_coupled_core_displacement.max(core_report.max_displacement_from_centroid);
        let (finite, invalid) = count_control_points(&surface);
        finite_control_point_count += finite;
        invalid_control_point_count += invalid;
        inner_patches.push(Kp20MultiPatchInnerPatch {
            patch_name: sector.inner_patch_name.clone(),
            sector_index: sector.sector_index,
            boundary_interface_names: sector_interfaces
                .iter()
                .map(|interface| interface.interface_name.clone())
                .collect(),
            surface,
        });
    }

    let inner_patch_names = inner_patches
        .iter()
        .map(|patch| patch.patch_name.clone())
        .collect::<Vec<_>>();
    let decline_reason = if inner_patches.len() != inventory.attempted_inner_patch_count {
        Some("missing_inner_patch")
    } else if invalid_control_point_count > 0 {
        Some("non_finite_inner_control_point")
    } else if max_boundary_g0_residual > EXACT_BOUNDARY_G0_TOLERANCE {
        Some("inner_boundary_g0_residual")
    } else if max_adjacent_row_residual > EXACT_BOUNDARY_G0_TOLERANCE {
        Some("inner_adjacent_row_residual")
    } else if max_interior_residual > EXACT_BOUNDARY_G0_TOLERANCE {
        Some("inner_core_residual")
    } else if !min_adjacent_row_step.is_finite() {
        Some("missing_inner_adjacent_row")
    } else {
        None
    };
    Kp20MultiPatchInnerBuildReport {
        attempted_inner_patch_count: inventory.attempted_inner_patch_count,
        built_inner_patch_count: inner_patches.len(),
        finite_control_point_count,
        invalid_control_point_count,
        boundary_constraint_count,
        adjacent_row_constraint_count,
        interior_constraint_count,
        exact_boundary_owned_endpoint_count,
        max_boundary_g0_residual,
        max_adjacent_row_residual,
        max_interior_residual,
        min_adjacent_row_step,
        max_coupled_core_displacement,
        decline_reason,
        inner_patch_names,
        inner_patches,
    }
}

fn interface_by_name<'a>(
    interfaces: &'a Kp20MultiPatchInterfaceBuildReport,
    name: &str,
) -> Option<&'a Kp20MultiPatchInterfaceCurve> {
    interfaces
        .interface_curves
        .iter()
        .find(|interface| interface.interface_name == name)
}

fn interface_centroid(interfaces: &[&Kp20MultiPatchInterfaceCurve]) -> Point3 {
    let mut x = 0.0_f64;
    let mut y = 0.0_f64;
    let mut z = 0.0_f64;
    let mut count = 0_usize;
    for interface in interfaces {
        for point in interface.controls {
            x += point.x;
            y += point.y;
            z += point.z;
            count += 1;
        }
    }
    if count == 0 {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let scale = 1.0 / count as f64;
    Point3::new(x * scale, y * scale, z * scale)
}

fn assign_boundary_controls(
    surface: &mut Kp20Bi6Patch,
    boundary_index: usize,
    controls: [Point3; 7],
) {
    for (idx, control) in controls.into_iter().enumerate() {
        match boundary_index {
            0 => surface.cps[0][idx] = control,
            1 => surface.cps[6][idx] = control,
            2 => surface.cps[idx][0] = control,
            _ => surface.cps[idx][6] = control,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Kp20InnerAdjacentRowReport {
    max_construction_residual: f64,
    min_step: f64,
}

fn assign_boundary_adjacent_controls(
    surface: &mut Kp20Bi6Patch,
    boundary_index: usize,
    controls: [Point3; 7],
    support_adjacent_controls: [Point3; 7],
    boundary_count: usize,
) {
    for idx in 0..7 {
        if endpoint_exact_boundary_owner(boundary_index, idx)
            .is_some_and(|owner| owner < boundary_count)
        {
            continue;
        }
        let expected =
            support_tangent_adjacent_control(controls[idx], support_adjacent_controls[idx]);
        assign_adjacent_control(surface, boundary_index, idx, expected);
    }
}

fn adjacent_row_residual(
    surface: &Kp20Bi6Patch,
    boundary_index: usize,
    controls: [Point3; 7],
    support_adjacent_controls: [Point3; 7],
    boundary_count: usize,
) -> Kp20InnerAdjacentRowReport {
    let mut max_construction_residual = 0.0_f64;
    let mut min_step = f64::INFINITY;
    for idx in 0..7 {
        if endpoint_exact_boundary_owner(boundary_index, idx)
            .is_some_and(|owner| owner < boundary_count)
        {
            continue;
        }
        let boundary = controls[idx];
        let expected = support_tangent_adjacent_control(boundary, support_adjacent_controls[idx]);
        let assigned = adjacent_control(surface, boundary_index, idx);
        max_construction_residual = max_construction_residual.max(assigned.distance_to(&expected));
        min_step = min_step.min(assigned.distance_to(&boundary));
    }
    Kp20InnerAdjacentRowReport {
        max_construction_residual,
        min_step,
    }
}

fn exact_boundary_owned_endpoint_count_for(boundary_index: usize, boundary_count: usize) -> usize {
    [0_usize, 6]
        .into_iter()
        .filter(|idx| {
            endpoint_exact_boundary_owner(boundary_index, *idx)
                .is_some_and(|owner| owner < boundary_count)
        })
        .count()
}

fn adjacent_frame_centroid(surface: &Kp20Bi6Patch) -> Point3 {
    let mut x = 0.0_f64;
    let mut y = 0.0_f64;
    let mut z = 0.0_f64;
    let mut count = 0_usize;
    for u in 0..7 {
        for v in 0..7 {
            if is_inner_core_control(u, v) {
                continue;
            }
            let point = surface.cps[u][v];
            x += point.x;
            y += point.y;
            z += point.z;
            count += 1;
        }
    }
    let scale = 1.0 / count as f64;
    Point3::new(x * scale, y * scale, z * scale)
}

#[derive(Debug, Clone, Copy)]
struct Kp20InnerCoreConstructionReport {
    max_displacement_from_centroid: f64,
}

fn assign_inner_core_controls(surface: &mut Kp20Bi6Patch) -> Kp20InnerCoreConstructionReport {
    let centroid = adjacent_frame_centroid(surface);
    let mut max_displacement_from_centroid = 0.0_f64;
    for u in 2..=4 {
        for v in 2..=4 {
            let target = coupled_inner_core_target(surface, u, v);
            surface.cps[u][v] = target;
            max_displacement_from_centroid =
                max_displacement_from_centroid.max(target.distance_to(&centroid));
        }
    }
    Kp20InnerCoreConstructionReport {
        max_displacement_from_centroid,
    }
}

fn inner_core_residual(surface: &Kp20Bi6Patch) -> f64 {
    let mut residual = 0.0_f64;
    for u in 2..=4 {
        for v in 2..=4 {
            let target = coupled_inner_core_target(surface, u, v);
            residual = residual.max(surface.cps[u][v].distance_to(&target));
        }
    }
    residual
}

fn coupled_inner_core_target(surface: &Kp20Bi6Patch, u: usize, v: usize) -> Point3 {
    let su = (u - 2) as f64 / 2.0;
    let sv = (v - 2) as f64 / 2.0;
    let horizontal = surface.cps[u][1] + (surface.cps[u][5] - surface.cps[u][1]) * sv;
    let vertical = surface.cps[1][v] + (surface.cps[5][v] - surface.cps[1][v]) * su;
    horizontal + (vertical - horizontal) * 0.5
}

fn is_inner_core_control(u: usize, v: usize) -> bool {
    (2..=4).contains(&u) && (2..=4).contains(&v)
}

fn endpoint_exact_boundary_owner(boundary_index: usize, idx: usize) -> Option<usize> {
    match (boundary_index, idx) {
        (0, 0) | (1, 0) => Some(2),
        (0, 6) | (1, 6) => Some(3),
        (2, 0) | (3, 0) => Some(0),
        (2, 6) | (3, 6) => Some(1),
        _ => None,
    }
}

fn support_tangent_adjacent_control(boundary: Point3, support_adjacent: Point3) -> Point3 {
    boundary + (boundary - support_adjacent)
}

fn assign_adjacent_control(
    surface: &mut Kp20Bi6Patch,
    boundary_index: usize,
    idx: usize,
    control: Point3,
) {
    match boundary_index {
        0 => surface.cps[1][idx] = control,
        1 => surface.cps[5][idx] = control,
        2 => surface.cps[idx][1] = control,
        _ => surface.cps[idx][5] = control,
    }
}

fn adjacent_control(surface: &Kp20Bi6Patch, boundary_index: usize, idx: usize) -> Point3 {
    match boundary_index {
        0 => surface.cps[1][idx],
        1 => surface.cps[5][idx],
        2 => surface.cps[idx][1],
        _ => surface.cps[idx][5],
    }
}

fn boundary_residual(
    surface: &Kp20Bi6Patch,
    boundary_index: usize,
    interface: &Kp20MultiPatchInterfaceCurve,
) -> f64 {
    let mut residual = 0.0_f64;
    for (idx, control) in interface.controls.iter().enumerate() {
        let assigned = match boundary_index {
            0 => surface.cps[0][idx],
            1 => surface.cps[6][idx],
            2 => surface.cps[idx][0],
            _ => surface.cps[idx][6],
        };
        residual = residual.max(assigned.distance_to(control));
    }
    residual
}

fn count_control_points(surface: &Kp20Bi6Patch) -> (usize, usize) {
    let mut finite = 0_usize;
    let mut invalid = 0_usize;
    for row in surface.cps {
        for point in row {
            if point.x.is_finite() && point.y.is_finite() && point.z.is_finite() {
                finite += 1;
            } else {
                invalid += 1;
            }
        }
    }
    (finite, invalid)
}

const EXACT_BOUNDARY_G0_TOLERANCE: f64 = 1.0e-12;
const ADJACENT_CONTROL_COUNT: usize = 7;
const INNER_CORE_CONTROL_COUNT: usize = 9;

#[cfg(test)]
#[path = "kp20_polar_multipatch_inner_tests.rs"]
mod tests;
