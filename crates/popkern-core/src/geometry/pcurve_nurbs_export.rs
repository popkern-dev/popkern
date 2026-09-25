//! NURBS export helpers for analytic pcurve variants.

use std::f64::consts::{PI, TAU};

use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::Pcurve;

pub(super) fn sampled_to_nurbs_approx(pcurve: &Pcurve) -> NurbsCurve {
    // `interpolate_curve` requires ascending parameters; a pcurve
    // trimmed against a descending domain (t_start > t_end) is
    // fitted through the domain-swapped clone instead. Analytic
    // variants evaluate by FRACTION of (t_start → t_end), so the
    // swapped clone traverses the same UV path in the same order —
    // only its knot domain ascends, which is what interpolation
    // needs. (PK-215: the descending case previously fed the solve
    // a descending knot vector, failed, and the silent fallback
    // shipped a chord for a quarter-arc plane trim — 2.9e-2 sagitta
    // in STEP.)
    if pcurve.t_start > pcurve.t_end {
        let ascending = Pcurve {
            geometry: pcurve.geometry.clone(),
            t_start: pcurve.t_end,
            t_end: pcurve.t_start,
        };
        return sampled_to_nurbs_with_points(&ascending, sampled_points(&ascending, 17));
    }
    sampled_to_nurbs_with_points(pcurve, sampled_points(pcurve, 17))
}

/// Degree of the co-parameterized interpolant for analytic pcurves.
/// Uniform-node interpolation error falls as `h^(degree+1)`: a degree-5
/// fit of a radius-20 full circle reaches 1e-9 at 64 spans, where the
/// cubic 16-span fit was 6e-3 (PK-506 measurement).
const ANALYTIC_EXPORT_DEGREE: usize = 5;
/// Starting span count; halving the span improves a degree-5 fit 64x.
const ANALYTIC_EXPORT_INITIAL_SPANS: usize = 8;
/// Span cap: the dense interpolation solve is cubic in the node count,
/// and 256 spans settle any sweep of a smooth analytic trace far below
/// the tolerance floor.
const ANALYTIC_EXPORT_MAX_SPANS: usize = 256;
/// Error probes per span when measuring a fit against its source.
const ERROR_PROBES_PER_SPAN: usize = 4;

/// A NURBS form of an analytic pcurve that keeps its parameterization
/// (so it stays co-parameterized with the edge curve) and is measured
/// against the analytic source: degree-5 interpolation at uniform
/// nodes, with the span count doubled until the largest deviation at
/// interior probes is within `uv_tolerance` or the span cap is reached.
/// Returns the fit and its measured deviation (the best fit found when
/// the tolerance is out of reach).
pub(super) fn analytic_to_nurbs_within(pcurve: &Pcurve, uv_tolerance: f64) -> (NurbsCurve, f64) {
    let ascending = if pcurve.t_start > pcurve.t_end {
        Pcurve {
            geometry: pcurve.geometry.clone(),
            t_start: pcurve.t_end,
            t_end: pcurve.t_start,
        }
    } else {
        pcurve.clone()
    };
    let mut best: Option<(NurbsCurve, f64)> = None;
    let mut spans = ANALYTIC_EXPORT_INITIAL_SPANS;
    loop {
        let points = sampled_points(&ascending, spans + 1);
        let params = uniform_params(&ascending, spans + 1);
        let Ok(fit) =
            crate::geometry::fitting::interpolate_curve(&points, &params, ANALYTIC_EXPORT_DEGREE)
        else {
            break;
        };
        let deviation = max_deviation(&fit, &ascending, spans);
        if best.as_ref().is_none_or(|(_, d)| deviation < *d) {
            best = Some((fit, deviation));
        }
        if deviation <= uv_tolerance || spans >= ANALYTIC_EXPORT_MAX_SPANS {
            break;
        }
        spans *= 2;
    }
    best.unwrap_or_else(|| {
        let fallback = sampled_to_nurbs_approx(&ascending);
        let deviation = max_deviation(&fallback, &ascending, 16);
        (fallback, deviation)
    })
}

fn uniform_params(pcurve: &Pcurve, samples: usize) -> Vec<f64> {
    (0..samples)
        .map(|i| {
            let s = i as f64 / (samples as f64 - 1.0);
            pcurve.t_start + s * (pcurve.t_end - pcurve.t_start)
        })
        .collect()
}

/// Largest uv distance between `fit` and `source` at the same
/// parameter, probed inside each of `spans` uniform spans.
fn max_deviation(fit: &NurbsCurve, source: &Pcurve, spans: usize) -> f64 {
    use crate::geometry::curve::CurveEval;
    let probes = spans * ERROR_PROBES_PER_SPAN;
    let mut worst = 0.0_f64;
    for k in 0..probes {
        let s = (k as f64 + 0.5) / probes as f64;
        let t = source.t_start + s * (source.t_end - source.t_start);
        worst = worst.max(fit.point_at(t).distance_to(&source.point_at(t)));
    }
    worst
}

pub(super) fn cylinder_arc_via_rational_quad_to_nurbs(pcurve: &Pcurve) -> NurbsCurve {
    let start = pcurve.point_at(pcurve.t_start);
    let mut end = pcurve.point_at(pcurve.t_end);
    unwrap_u_relative_to(start.x, &mut end.x);
    NurbsCurve {
        degree: 1,
        knots: vec![pcurve.t_start, pcurve.t_start, pcurve.t_end, pcurve.t_end],
        control_points: vec![start, end],
        weights: vec![1.0, 1.0],
    }
}

fn sampled_to_nurbs_with_points(pcurve: &Pcurve, control_points: Vec<Point3>) -> NurbsCurve {
    let params = (0..control_points.len())
        .map(|i| {
            let s = i as f64 / (control_points.len() as f64 - 1.0);
            pcurve.t_start + s * (pcurve.t_end - pcurve.t_start)
        })
        .collect::<Vec<_>>();
    match crate::geometry::fitting::interpolate_curve(&control_points, &params, 3) {
        Ok(c) => c,
        // Last-resort fallback: a degree-1 polyline THROUGH ALL
        // samples, not a two-point chord — a chord silently degrades
        // curved UV trims (a quarter-arc plane trim ships with its
        // full sagitta as pcurve error). The polyline stays within
        // the inter-sample sagitta of the true trace. Reached only
        // when interpolation fails on ascending parameters, which
        // requires degenerate/duplicated samples.
        Err(_) => polyline_nurbs(pcurve, control_points, params),
    }
}

fn polyline_nurbs(pcurve: &Pcurve, control_points: Vec<Point3>, params: Vec<f64>) -> NurbsCurve {
    let n = control_points.len();
    let strictly_ascending = n >= 2 && params.windows(2).all(|w| w[1] > w[0]);
    if !strictly_ascending {
        // Degenerate parameter domain: nothing better than the
        // endpoint segment exists.
        return NurbsCurve {
            degree: 1,
            knots: vec![pcurve.t_start, pcurve.t_start, pcurve.t_end, pcurve.t_end],
            control_points: vec![
                control_points[0],
                control_points.last().copied().unwrap_or(control_points[0]),
            ],
            weights: vec![1.0, 1.0],
        };
    }
    let mut knots = Vec::with_capacity(n + 2);
    knots.push(params[0]);
    knots.extend(params.iter().copied());
    knots.push(params[n - 1]);
    NurbsCurve {
        degree: 1,
        knots,
        weights: vec![1.0; n],
        control_points,
    }
}

fn sampled_points(pcurve: &Pcurve, samples: usize) -> Vec<Point3> {
    let samples = samples.max(2);
    (0..samples)
        .map(|i| {
            let s = i as f64 / (samples as f64 - 1.0);
            let t = pcurve.t_start + s * (pcurve.t_end - pcurve.t_start);
            pcurve.point_at(t)
        })
        .collect()
}

fn unwrap_u_relative_to(reference_u: f64, u: &mut f64) {
    while *u - reference_u > PI {
        *u -= TAU;
    }
    while *u - reference_u < -PI {
        *u += TAU;
    }
}

#[cfg(test)]
#[path = "pcurve_nurbs_export_tests.rs"]
mod tests;
