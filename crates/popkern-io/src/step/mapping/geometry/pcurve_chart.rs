//! STEP → popkern pcurve chart maps for the analytic rotational
//! surfaces (split from `pcurves.rs`, 500-line cap).
//!
//! A STEP pcurve's coordinates live in its basis surface's chart: u is
//! an angle in the file's plane-angle unit measured from the
//! placement's ref direction; v is an angle (sphere latitude, torus)
//! or a length along the axis (cylinder, cone — from the placement,
//! where a cone's radius is R). popkern's charts measure u in radians
//! from `orthonormal_basis(axis)`, a sphere's v as colatitude and a
//! cone's v from the apex. The map between them is affine per axis
//! (`u' = u_scale·u + u_offset`, `v' = v_scale·v + v_offset`), then
//! the latitude reflection on the sphere.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_math::PortableMath;

/// STEP→popkern pcurve coordinate transform for one basis surface.
#[derive(Debug, Clone, Copy)]
pub(super) struct ChartTransform {
    pub u_scale: f64,
    pub u_offset: f64,
    pub v_scale: f64,
    pub v_offset: f64,
    pub latitude_flip: bool,
}

impl Default for ChartTransform {
    fn default() -> Self {
        ChartTransform {
            u_scale: 1.0,
            u_offset: 0.0,
            v_scale: 1.0,
            v_offset: 0.0,
            latitude_flip: false,
        }
    }
}

pub(super) fn apply_chart_transform(pcurve: Pcurve, transform: ChartTransform) -> Pcurve {
    let pcurve = map_pcurve_affine(pcurve, &transform);
    if !transform.latitude_flip {
        return pcurve;
    }
    flip_pcurve_latitude(pcurve)
}

/// `(u, v) → (u_scale·u + u_offset, v_scale·v + v_offset)`. Linear and
/// NURBS pcurves map exactly (an affine map of the control points);
/// with unit scales every form only translates. A scaled analytic
/// pcurve (an ellipse on a degree-unit chart) is mapped through its
/// NURBS form.
fn map_pcurve_affine(mut pcurve: Pcurve, t: &ChartTransform) -> Pcurve {
    let map = |u: f64, v: f64| (t.u_scale * u + t.u_offset, t.v_scale * v + t.v_offset);
    let unit_scale = t.u_scale == 1.0 && t.v_scale == 1.0;
    match &mut pcurve.geometry {
        PcurveGeometry::Linear { start, end } => {
            *start = map(start.0, start.1);
            *end = map(end.0, end.1);
        }
        PcurveGeometry::Nurbs(nurbs) => {
            for point in &mut nurbs.control_points {
                (point.x, point.y) = map(point.x, point.y);
            }
        }
        PcurveGeometry::Elliptical { center, .. } if unit_scale => {
            *center = map(center.0, center.1);
        }
        PcurveGeometry::Sinusoidal { u_offset, .. } if unit_scale && t.v_offset == 0.0 => {
            *u_offset += t.u_offset;
        }
        PcurveGeometry::CylinderArcViaRationalQuad {
            theta_a, theta_b, ..
        } if unit_scale && t.v_offset == 0.0 => {
            *theta_a += t.u_offset;
            *theta_b += t.u_offset;
        }
        _ => {
            let (mut nurbs, _) = pcurve.to_nurbs_within(1e-12);
            for point in &mut nurbs.control_points {
                (point.x, point.y) = map(point.x, point.y);
            }
            return Pcurve::from(nurbs);
        }
    }
    pcurve
}

/// v → π/2 − v (STEP sphere latitude to popkern colatitude). Import
/// produces Linear, NURBS and (from 2-D conics) Elliptical pcurves (see
/// `resolve_pcurve_on_surface`); all three are closed under this
/// reflection. An ellipse reflects exactly by negating its rotation and
/// sweep angles: `c + R(θ)(a cos α, b sin α)` maps to
/// `c' + R(−θ)(a cos(−α), b sin(−α))`.
fn flip_pcurve_latitude(mut pcurve: Pcurve) -> Pcurve {
    let half_pi = std::f64::consts::FRAC_PI_2;
    match &mut pcurve.geometry {
        PcurveGeometry::Linear { start, end } => {
            start.1 = half_pi - start.1;
            end.1 = half_pi - end.1;
        }
        PcurveGeometry::Nurbs(nurbs) => {
            for point in &mut nurbs.control_points {
                point.y = half_pi - point.y;
            }
        }
        PcurveGeometry::Elliptical {
            center,
            rotation,
            start_angle,
            end_angle,
            ..
        } => {
            center.1 = half_pi - center.1;
            *rotation = -*rotation;
            *start_angle = -*start_angle;
            *end_angle = -*end_angle;
        }
        // Not produced by import; leave unchanged rather than guessing
        // a convention.
        PcurveGeometry::Sinusoidal { .. } | PcurveGeometry::CylinderArcViaRationalQuad { .. } => {}
    }
    pcurve
}

pub(super) fn analytic_u_frame_offset(axis: Vec3, step_ref_direction: Vec3) -> Option<f64> {
    let axis = axis.normalized()?;
    let step_u = (step_ref_direction - axis * step_ref_direction.dot(&axis)).normalized()?;
    let (kernel_u, kernel_v) = orthonormal_basis(axis);
    Some(step_u.dot(&kernel_v).patan2(step_u.dot(&kernel_u)))
}
