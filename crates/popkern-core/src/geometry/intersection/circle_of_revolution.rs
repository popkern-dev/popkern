//! Circle-of-revolution emission shared by the analytic torus
//! intersection arms (PK-279).
//!
//! Every axis-aligned torus intersection class — coaxial cylinder,
//! cone, sphere or torus, and the axis-perpendicular plane — meets the
//! torus in exact circles centred on the shared axis. This module owns
//! the shared emission machinery:
//!
//! - the exact rational-quadratic NURBS circle;
//! - the co-parameterized chart pcurve fit (the `plane_sphere`
//!   pattern: chart the 3D curve at its own parameters, unwrap by
//!   continuity, interpolate degree-1 with knots AT those parameters,
//!   and densify until mid-knot verification passes) generalized to
//!   any analytic chart with a closed-form point inverse;
//! - the exact affine plane pcurve projection;
//! - the periodic window containment test used to honour the caller's
//!   parameter windows for circle stations.

use crate::geometry::curve::CurveEval;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::{FRAC_1_SQRT_2, TAU};

/// Exact rational-quadratic full circle (parameter 0..1) centred at
/// `c` in the plane spanned by the orthonormal pair (`a1`, `a2`).
pub(crate) fn nurbs_circle(c: Point3, a1: Vec3, a2: Vec3, r: f64) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            c + a1 * r,
            c + a1 * r + a2 * r,
            c + a2 * r,
            c - a1 * r + a2 * r,
            c - a1 * r,
            c - a1 * r - a2 * r,
            c - a2 * r,
            c + a1 * r - a2 * r,
            c + a1 * r,
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Build the pcurve of `curve_3d` on `surface`, co-parameterized so
/// that `surface(uv(t)) ≈ curve_3d(t)` for the same `t`.
///
/// Planes get the exact affine control-point projection (degree,
/// knots, and weights preserved — projection commutes with the
/// rational combination). Every other analytic chart is fitted with
/// the adaptive degree-1 interpolation of [`fit_chart_uv_curve`].
/// Returns `None` when the chart inverse fails at the first sample —
/// callers treat that as "decline the analytic arm" and fall back to
/// marching.
pub(crate) fn pcurve_on(
    surface: &SurfaceGeometry,
    curve_3d: &NurbsCurve,
    u_window: (f64, f64),
    v_window: (f64, f64),
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    match surface {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => Some(plane_uv_of_curve(curve_3d, *origin, *u_axis, *v_axis)),
        _ => {
            let mut pcurve = fit_chart_uv_curve(surface, curve_3d, tolerance)?;
            let (period_u, period_v) = surface.parameter_periodicity();
            align_pcurve_to_windows(&mut pcurve, u_window, v_window, period_u, period_v);
            Some(pcurve)
        }
    }
}

/// Shift a fitted chart pcurve by whole periods per axis so it lands
/// in the CALLER's chart branch. The marching path wraps every sample
/// into the caller's parameter window (`wrap_parameter` in
/// `marching.rs`), and downstream UV clipping matches pcurves against
/// face trim domains in that same branch — a canonical-branch pcurve
/// a period away silently misses the face and the boolean split stage
/// fails endpoint matching (measured: rim-fillet strip windows
/// authored on a negative-v branch vs the chart-canonical `[0, TAU)`
/// fit). The shift is a single constant per axis, preserving the
/// fit's continuity.
pub(crate) fn align_pcurve_to_windows(
    pcurve: &mut NurbsCurve,
    u_window: (f64, f64),
    v_window: (f64, f64),
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    let shift = |values: &mut dyn Iterator<Item = &mut f64>, window: (f64, f64), period: f64| {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        let mut collected: Vec<&mut f64> = Vec::new();
        for v in values {
            lo = lo.min(*v);
            hi = hi.max(*v);
            collected.push(v);
        }
        if !lo.is_finite() || !hi.is_finite() {
            return;
        }
        let window_center = 0.5 * (window.0 + window.1);
        let chain_center = 0.5 * (lo + hi);
        let k = ((window_center - chain_center) / period).round();
        if k != 0.0 {
            for v in collected {
                *v += k * period;
            }
        }
    };
    if let Some(p) = period_u {
        shift(
            &mut pcurve.control_points.iter_mut().map(|cp| &mut cp.x),
            u_window,
            p,
        );
    }
    if let Some(p) = period_v {
        shift(
            &mut pcurve.control_points.iter_mut().map(|cp| &mut cp.y),
            v_window,
            p,
        );
    }
}

/// Project a 3D NURBS curve into a plane's UV chart by transforming
/// each control point (exact: the chart map is affine).
pub(crate) fn plane_uv_of_curve(
    curve: &NurbsCurve,
    origin: Point3,
    u_axis: Vec3,
    v_axis: Vec3,
) -> NurbsCurve {
    let u_len_sq = u_axis.dot(&u_axis);
    let v_len_sq = v_axis.dot(&v_axis);
    let cps = curve
        .control_points
        .iter()
        .map(|p| {
            let d = *p - origin;
            let u = if u_len_sq > 1e-30 {
                d.dot(&u_axis) / u_len_sq
            } else {
                d.dot(&u_axis)
            };
            let v = if v_len_sq > 1e-30 {
                d.dot(&v_axis) / v_len_sq
            } else {
                d.dot(&v_axis)
            };
            Point3::new(u, v, 0.0)
        })
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: cps,
        weights: curve.weights.clone(),
    }
}

/// Fit the chart image of `curve_3d` on `surface` as a
/// co-parameterized degree-1 NURBS — the shared
/// [`crate::geometry::chart_pcurve::fit_chart_pcurve`] over the
/// curve's own parameter range.
pub(crate) fn fit_chart_uv_curve(
    surface: &SurfaceGeometry,
    curve_3d: &NurbsCurve,
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    crate::geometry::chart_pcurve::fit_chart_pcurve(
        surface,
        curve_3d,
        curve_3d.parameter_range(),
        tolerance,
    )
}

/// Assemble a full-circle intersection curve with co-parameterized
/// pcurves on both operand surfaces (in operand order), each aligned
/// to its caller parameter windows. Returns `None` when either pcurve
/// cannot be charted — callers decline the analytic arm and keep the
/// marching fallback.
#[allow(clippy::too_many_arguments)]
pub(crate) fn circle_curve_on(
    surface1: &SurfaceGeometry,
    windows1: ((f64, f64), (f64, f64)),
    surface2: &SurfaceGeometry,
    windows2: ((f64, f64), (f64, f64)),
    center: Point3,
    a1: Vec3,
    a2: Vec3,
    radius: f64,
    tolerance: &Tolerance,
) -> Option<CurveOnSurface> {
    let curve_3d = nurbs_circle(center, a1, a2, radius);
    let curve_uv1 = pcurve_on(surface1, &curve_3d, windows1.0, windows1.1, tolerance)?;
    let curve_uv2 = pcurve_on(surface2, &curve_3d, windows2.0, windows2.1, tolerance)?;
    Some(CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed: true,
    })
}

/// True when `x` (or some period translate of it) lies inside
/// `window`. For non-periodic parameters pass `period = None` — the
/// test degenerates to plain interval containment. The window is
/// normalized so a reversed pair still bounds correctly, and windows
/// spanning a full period accept everything.
pub(crate) fn station_in_window(x: f64, window: (f64, f64), period: Option<f64>) -> bool {
    let lo = window.0.min(window.1);
    let hi = window.0.max(window.1);
    match period {
        None => x >= lo && x <= hi,
        Some(p) => {
            if hi - lo >= p - 1e-12 {
                return true;
            }
            // Shift x by whole periods into [lo, lo + p) and test.
            let shifted = x - p * ((x - lo) / p).floor();
            shifted <= hi
        }
    }
}

/// Angular station of a profile point in the torus tube chart —
/// `atan2` of (axial offset, radial distance minus major radius),
/// matching the kernel torus chart convention (`project_point_to_uv`).
/// `rho` is the point's UNSIGNED distance from the axis.
pub(crate) fn torus_tube_station(rho: f64, h: f64, major_radius: f64) -> f64 {
    let mut v = h.atan2(rho - major_radius);
    if v < 0.0 {
        v += TAU;
    }
    v
}
