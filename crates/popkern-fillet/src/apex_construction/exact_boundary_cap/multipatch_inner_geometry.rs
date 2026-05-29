//! Geometry helpers for multi-patch inner pole candidates.
#![allow(dead_code)]

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use super::inner::{
    ExactBoundaryInnerPolePatchPrototype, ExactBoundaryInterfaceG1Report,
    ExactBoundaryStripInnerInterfaceCurve,
};
use super::strip::{ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryStripBoundaryDirection};

pub(super) const INTERFACE_SAMPLE_COUNT: usize = 17;
pub(super) const POLE_SAMPLE_COUNT: usize = 9;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct InterfaceResidual {
    pub(super) max_residual: f64,
    pub(super) sample_parameter: f64,
}

pub(super) fn strip_inner_boundary_curve(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
) -> NurbsCurve {
    let control_points = (0..candidate.surface.num_u)
        .map(|index| candidate.surface.control_points[index * candidate.surface.num_v + 1])
        .collect::<Vec<_>>();
    let weights = (0..candidate.surface.num_u)
        .map(|index| candidate.surface.weights[index * candidate.surface.num_v + 1])
        .collect::<Vec<_>>();
    NurbsCurve {
        degree: candidate.surface.degree_u,
        knots: candidate.surface.knots_u.clone(),
        control_points,
        weights,
    }
}

pub(super) fn curve_geometry_is_finite(curve: &NurbsCurve) -> bool {
    curve.knots.iter().all(|value| value.is_finite())
        && curve
            .weights
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
        && curve
            .control_points
            .iter()
            .all(|point| point_is_finite(*point))
}

pub(super) fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

pub(super) fn strip_inner_boundary_pcurve(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    curve: &NurbsCurve,
) -> Pcurve {
    let (t_start, t_end) = curve.parameter_range();
    match candidate.boundary_direction {
        ExactBoundaryStripBoundaryDirection::Forward => {
            Pcurve::linear((t_start, 1.0), (t_end, 1.0), t_start, t_end)
        }
        ExactBoundaryStripBoundaryDirection::Reversed => {
            Pcurve::linear((t_end, 1.0), (t_start, 1.0), t_start, t_end)
        }
    }
}

pub(super) fn interface_residual(
    candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    curve: &NurbsCurve,
    strip_pcurve: &Pcurve,
) -> InterfaceResidual {
    let (t_start, t_end) = curve.parameter_range();
    let mut report = InterfaceResidual {
        max_residual: 0.0,
        sample_parameter: t_start,
    };
    for index in 0..INTERFACE_SAMPLE_COUNT {
        let parameter =
            t_start + (t_end - t_start) * index as f64 / (INTERFACE_SAMPLE_COUNT - 1) as f64;
        let uv = strip_pcurve.point_at(parameter);
        let strip_point = candidate.surface.point_at(uv.x, uv.y);
        let curve_point = curve.point_at(uv.x);
        let residual = strip_point.distance_to(&curve_point);
        if residual > report.max_residual {
            report = InterfaceResidual {
                max_residual: residual,
                sample_parameter: parameter,
            };
        }
    }
    report
}

pub(super) fn inner_pole_surface(
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    pole: Point3,
) -> NurbsSurface {
    let mut control_points = Vec::with_capacity(interface.curve.control_points.len() * 3);
    let mut weights = Vec::with_capacity(interface.curve.weights.len() * 3);
    for (index, (point, weight)) in interface
        .curve
        .control_points
        .iter()
        .zip(&interface.curve.weights)
        .enumerate()
    {
        control_points.push(*point);
        control_points.push(interface_normal_row_point(
            strip_candidate,
            &interface.curve,
            index,
            *point,
            pole,
        ));
        control_points.push(pole);
        weights.push(*weight);
        weights.push(*weight);
        weights.push(1.0);
    }
    NurbsSurface {
        degree_u: interface.curve.degree,
        degree_v: 2,
        knots_u: interface.curve.knots.clone(),
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points,
        num_u: interface.curve.control_points.len(),
        num_v: 3,
        weights,
    }
}

fn interface_normal_row_point(
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    interface_curve: &NurbsCurve,
    index: usize,
    point: Point3,
    pole: Point3,
) -> Point3 {
    let parameter = control_row_parameter(interface_curve, index);
    let tangent = interface_curve
        .derivative_at(parameter, 1)
        .normalized()
        .unwrap_or_else(|| fallback_tangent(point, pole));
    let target_normal = strip_candidate
        .surface
        .normal_at(parameter, 1.0)
        .normalized()
        .unwrap_or_else(|| fallback_normal(tangent, point, pole));
    let direction = target_normal
        .cross(&tangent)
        .normalized()
        .unwrap_or_else(|| fallback_tangent(point, pole));
    let width = point.distance_to(&pole).max(1.0e-6) * 0.35;
    point + direction * width
}

fn control_row_parameter(curve: &NurbsCurve, index: usize) -> f64 {
    let (start, end) = curve.parameter_range();
    let count = curve.control_points.len();
    if count <= 1 {
        return start;
    }
    start + (end - start) * index as f64 / (count - 1) as f64
}

fn fallback_tangent(point: Point3, pole: Point3) -> Vec3 {
    (pole - point)
        .normalized()
        .unwrap_or_else(|| Vec3::new(1.0, 0.0, 0.0))
}

fn fallback_normal(tangent: Vec3, point: Point3, pole: Point3) -> Vec3 {
    tangent
        .cross(&(pole - point))
        .normalized()
        .unwrap_or_else(|| Vec3::new(0.0, 0.0, 1.0))
}

pub(super) fn inner_interface_pcurve(interface: &ExactBoundaryStripInnerInterfaceCurve) -> Pcurve {
    let (t_start, t_end) = interface.curve.parameter_range();
    if interface.strip_parameter_reversed {
        Pcurve::linear((t_end, 0.0), (t_start, 0.0), t_start, t_end)
    } else {
        Pcurve::linear((t_start, 0.0), (t_end, 0.0), t_start, t_end)
    }
}

pub(super) fn inner_interface_g0(
    surface: &NurbsSurface,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    pcurve: &Pcurve,
) -> InterfaceResidual {
    let (t_start, t_end) = interface.curve.parameter_range();
    let mut report = InterfaceResidual {
        max_residual: 0.0,
        sample_parameter: t_start,
    };
    for index in 0..INTERFACE_SAMPLE_COUNT {
        let parameter =
            t_start + (t_end - t_start) * index as f64 / (INTERFACE_SAMPLE_COUNT - 1) as f64;
        let uv = pcurve.point_at(parameter);
        let patch_point = surface.point_at(uv.x, uv.y);
        let interface_point = interface.curve.point_at(uv.x);
        let residual = patch_point.distance_to(&interface_point);
        if residual > report.max_residual {
            report = InterfaceResidual {
                max_residual: residual,
                sample_parameter: parameter,
            };
        }
    }
    report
}

pub(super) fn interface_g1_report(
    strip_candidate: &ExactBoundaryOuterStripBoundaryScaffold,
    interface: &ExactBoundaryStripInnerInterfaceCurve,
    prototype: &ExactBoundaryInnerPolePatchPrototype,
    tolerance: f64,
) -> ExactBoundaryInterfaceG1Report {
    let (t_start, t_end) = interface.curve.parameter_range();
    let mut report = ExactBoundaryInterfaceG1Report {
        max_interface_g1_deficit: 0.0,
        sample_parameter: t_start,
        sample_count: INTERFACE_SAMPLE_COUNT,
        tolerance,
    };
    for index in 0..INTERFACE_SAMPLE_COUNT {
        let parameter =
            t_start + (t_end - t_start) * index as f64 / (INTERFACE_SAMPLE_COUNT - 1) as f64;
        let strip_uv = interface.strip_pcurve.point_at(parameter);
        let inner_uv = prototype.interface_pcurve.point_at(parameter);
        let deficit = tangent_plane_deficit(
            &strip_candidate.surface,
            strip_uv.x,
            strip_uv.y,
            &prototype.surface,
            inner_uv.x,
            inner_uv.y,
        )
        .unwrap_or(f64::INFINITY);
        if deficit > report.max_interface_g1_deficit {
            report.max_interface_g1_deficit = deficit;
            report.sample_parameter = parameter;
        }
    }
    report
}

fn tangent_plane_deficit(
    strip_surface: &NurbsSurface,
    strip_u: f64,
    strip_v: f64,
    inner_surface: &NurbsSurface,
    inner_u: f64,
    inner_v: f64,
) -> Option<f64> {
    let strip_normal = safe_surface_normal(strip_surface, strip_u, strip_v)?;
    let inner_normal = safe_surface_normal(inner_surface, inner_u, inner_v)?;
    Some(1.0 - strip_normal.dot(&inner_normal).abs())
}

fn safe_surface_normal(surface: &NurbsSurface, u: f64, v: f64) -> Option<Vec3> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u = interior_parameter(u, u_min, u_max);
    let v = interior_parameter(v, v_min, v_max);
    surface.normal_at(u, v).normalized()
}

fn interior_parameter(parameter: f64, min: f64, max: f64) -> f64 {
    let span = max - min;
    if span.abs() <= 1.0e-12 {
        return parameter;
    }
    let inset = span.abs() * 1.0e-6;
    parameter.clamp(min + inset, max - inset)
}

pub(super) fn pole_g0(surface: &NurbsSurface, pole: Point3) -> InterfaceResidual {
    let ((t_start, t_end), _) = surface.parameter_range();
    let mut report = InterfaceResidual {
        max_residual: 0.0,
        sample_parameter: t_start,
    };
    for index in 0..POLE_SAMPLE_COUNT {
        let parameter = t_start + (t_end - t_start) * index as f64 / (POLE_SAMPLE_COUNT - 1) as f64;
        let patch_point = surface.point_at(parameter, 1.0);
        let residual = patch_point.distance_to(&pole);
        if residual > report.max_residual {
            report = InterfaceResidual {
                max_residual: residual,
                sample_parameter: parameter,
            };
        }
    }
    report
}
