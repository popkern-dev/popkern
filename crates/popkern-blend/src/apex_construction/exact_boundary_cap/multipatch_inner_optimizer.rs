//! Bounded normal-row optimization for multi-patch inner pole candidates.
#![allow(dead_code)]

use popkern_core::geometry::{Point3, Vec3};

use super::inner::{ExactBoundaryInnerPolePatchPrototype, ExactBoundaryStripInnerInterfaceCurve};
use super::inner_geometry::{
    basis_function, inner_interface_g0, interface_g1_report, interface_normal_row_target, pole_g0,
};
use super::strip::ExactBoundaryOuterStripBoundaryScaffold;

const DAMPING_STEPS: [f64; 5] = [0.75, 0.5, 0.25, 0.125, 0.0625];
const IMPROVEMENT_EPSILON: f64 = 1.0e-12;

pub(super) fn optimize_inner_interface_g1(
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    tolerance: f64,
) -> ExactBoundaryInnerPolePatchPrototype {
    let mut best = prototype.clone();
    let mut best_report = interface_g1_report(strip_candidate, interface, &best, tolerance);
    if best_report.max_interface_g1_deficit <= tolerance {
        best.interface_g1_report = Some(best_report);
        best.interface_g1_validated = true;
        return best;
    }

    for _ in 0..16 {
        let Some(update) = normal_row_update(
            strip_candidate,
            interface,
            &best,
            best_report.sample_parameter,
        ) else {
            break;
        };
        let mut accepted = false;
        for damping in DAMPING_STEPS {
            let mut candidate = best.clone();
            apply_normal_row_update(&mut candidate, &update, damping);
            if !preserves_g0_invariants(interface, prototype, &candidate) {
                continue;
            }
            let report = interface_g1_report(strip_candidate, interface, &candidate, tolerance);
            if report.max_interface_g1_deficit + IMPROVEMENT_EPSILON
                < best_report.max_interface_g1_deficit
            {
                best = candidate;
                best_report = report;
                accepted = true;
                break;
            }
        }
        if best_report.max_interface_g1_deficit <= tolerance || !accepted {
            break;
        }
    }

    best.interface_g1_report = Some(best_report);
    best.interface_g1_validated = best_report.max_interface_g1_deficit <= tolerance;
    best
}

fn normal_row_update(
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    parameter: f64,
) -> Option<Vec<Vec3>> {
    let target = interface_normal_row_target(
        strip_candidate,
        interface,
        parameter,
        prototype.surface.control_points[prototype.surface.num_v - 1],
    );
    let current = normal_row_point_at(prototype, parameter)?;
    let delta = target - current;
    if delta.length() <= 1.0e-14 {
        return None;
    }
    let rational_basis = normal_row_rational_basis(prototype, parameter)?;
    let scale = rational_basis
        .iter()
        .map(|basis| basis * basis)
        .sum::<f64>();
    if scale <= 1.0e-14 {
        return None;
    }
    Some(
        rational_basis
            .into_iter()
            .map(|basis| delta * (basis / scale))
            .collect(),
    )
}

fn apply_normal_row_update(
    prototype: &mut ExactBoundaryInnerPolePatchPrototype,
    update: &[Vec3],
    damping: f64,
) {
    for (index, delta) in update.iter().enumerate() {
        let point_index = index * prototype.surface.num_v + 1;
        prototype.surface.control_points[point_index] =
            prototype.surface.control_points[point_index] + *delta * damping;
    }
    prototype.interface_g1_report = None;
    prototype.interface_g1_validated = false;
}

fn preserves_g0_invariants(
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    baseline: &ExactBoundaryInnerPolePatchPrototype,
    candidate: &ExactBoundaryInnerPolePatchPrototype,
) -> bool {
    if candidate.surface.num_u != baseline.surface.num_u
        || candidate.surface.num_v != baseline.surface.num_v
    {
        return false;
    }
    for index in 0..candidate.surface.num_u {
        let interface_index = index * candidate.surface.num_v;
        let pole_index = interface_index + candidate.surface.num_v - 1;
        if candidate.surface.control_points[interface_index]
            .distance_to(&baseline.surface.control_points[interface_index])
            > baseline.tolerance
        {
            return false;
        }
        if candidate.surface.control_points[pole_index]
            .distance_to(&baseline.surface.control_points[pole_index])
            > baseline.tolerance
        {
            return false;
        }
    }
    let interface_g0 =
        inner_interface_g0(&candidate.surface, interface, &candidate.interface_pcurve);
    let pole_g0 = pole_g0(&candidate.surface, candidate.surface.control_points[2]);
    interface_g0.max_residual <= baseline.tolerance && pole_g0.max_residual <= baseline.tolerance
}

fn normal_row_point_at(
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    parameter: f64,
) -> Option<Point3> {
    let basis = normal_row_rational_basis(prototype, parameter)?;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;
    for (index, rational_basis) in basis.into_iter().enumerate() {
        let point_index = index * prototype.surface.num_v + 1;
        let point = prototype.surface.control_points[point_index];
        x += point.x * rational_basis;
        y += point.y * rational_basis;
        z += point.z * rational_basis;
    }
    Some(Point3::new(x, y, z))
}

fn normal_row_rational_basis(
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    parameter: f64,
) -> Option<Vec<f64>> {
    let mut weighted = Vec::with_capacity(prototype.surface.num_u);
    let mut denominator = 0.0;
    for index in 0..prototype.surface.num_u {
        let basis = basis_function(
            index,
            prototype.surface.degree_u,
            parameter,
            &prototype.surface.knots_u,
        );
        let weight = prototype.surface.weights[index * prototype.surface.num_v + 1];
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

#[cfg(test)]
mod tests {
    use popkern_core::geometry::pcurve::Pcurve;
    use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};

    use super::*;
    use crate::apex_construction::exact_boundary_cap::multipatch::strip::ExactBoundaryStripBoundaryDirection;

    #[test]
    fn optimizer_moves_only_normal_row_until_interface_g1_passes() {
        let strip = flat_strip();
        let interface = straight_interface();
        let prototype = tilted_inner_patch();
        let original_interface = row_points(&prototype, 0);
        let original_pole = row_points(&prototype, 2);

        let optimized = optimize_inner_interface_g1(&strip, &interface, &prototype, 5.0e-2);
        let report = interface_g1_report(&strip, &interface, &optimized, 5.0e-2);

        assert!(optimized.interface_g1_validated);
        assert!(report.max_interface_g1_deficit <= 5.0e-2);
        assert_eq!(row_points(&optimized, 0), original_interface);
        assert_eq!(row_points(&optimized, 2), original_pole);
        assert!(
            inner_interface_g0(&optimized.surface, &interface, &optimized.interface_pcurve)
                .max_residual
                <= 1.0e-8
        );
        assert!(pole_g0(&optimized.surface, Point3::new(0.5, 1.0, 0.0)).max_residual <= 1.0e-8);
    }

    fn flat_strip() -> ExactBoundaryOuterStripBoundaryScaffold {
        ExactBoundaryOuterStripBoundaryScaffold {
            patch_name: "strip".to_string(),
            surface: NurbsSurface {
                degree_u: 1,
                degree_v: 1,
                knots_u: vec![0.0, 0.0, 1.0, 1.0],
                knots_v: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![
                    Point3::new(0.0, -1.0, 0.0),
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, -1.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                ],
                num_u: 2,
                num_v: 2,
                weights: vec![1.0; 4],
            },
            outer_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
            boundary_direction: ExactBoundaryStripBoundaryDirection::Forward,
            max_outer_boundary_residual: 0.0,
            outer_boundary_sample_count: 2,
            tolerance: 1.0e-8,
            rational_boundary: false,
            support_normal_row_pending: false,
            support_normal_report: None,
        }
    }

    fn straight_interface() -> ExactBoundaryStripInnerInterfaceCurve {
        ExactBoundaryStripInnerInterfaceCurve {
            interface_name: "interface".to_string(),
            strip_patch_name: "strip".to_string(),
            inner_patch_name: "inner".to_string(),
            curve: NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                weights: vec![1.0; 2],
            },
            strip_pcurve: Pcurve::linear((0.0, 1.0), (1.0, 1.0), 0.0, 1.0),
            strip_parameter_reversed: false,
            max_strip_boundary_residual: 0.0,
            sample_parameter: 0.0,
            sample_count: 17,
            tolerance: 1.0e-8,
        }
    }

    fn tilted_inner_patch() -> ExactBoundaryInnerPolePatchPrototype {
        ExactBoundaryInnerPolePatchPrototype {
            patch_name: "inner".to_string(),
            interface_name: "interface".to_string(),
            surface: NurbsSurface {
                degree_u: 1,
                degree_v: 2,
                knots_u: vec![0.0, 0.0, 1.0, 1.0],
                knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                control_points: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.0, 0.0, 0.2),
                    Point3::new(0.5, 1.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.2),
                    Point3::new(0.5, 1.0, 0.0),
                ],
                num_u: 2,
                num_v: 3,
                weights: vec![1.0; 6],
            },
            interface_pcurve: Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0),
            max_interface_g0: 0.0,
            max_pole_g0: 0.0,
            tolerance: 1.0e-8,
            interface_g1_report: None,
            interface_g1_validated: false,
        }
    }

    fn row_points(prototype: &ExactBoundaryInnerPolePatchPrototype, row: usize) -> Vec<Point3> {
        (0..prototype.surface.num_u)
            .map(|index| prototype.surface.control_points[index * prototype.surface.num_v + row])
            .collect()
    }
}
