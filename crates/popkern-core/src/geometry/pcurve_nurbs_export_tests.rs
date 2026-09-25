//! PK-550: analytic pcurves convert to co-parameterized NURBS within a
//! requested uv tolerance, checked densely against the analytic source.

use std::f64::consts::{FRAC_PI_2, PI, TAU};

use super::*;
use crate::geometry::curve::CurveEval;

/// Largest distance between `fit` and `source` at equal parameters,
/// over a probe grid far denser than the fit's own error probes.
fn dense_deviation(fit: &NurbsCurve, source: &Pcurve) -> f64 {
    let (lo, hi) = if source.t_start <= source.t_end {
        (source.t_start, source.t_end)
    } else {
        (source.t_end, source.t_start)
    };
    let ascending = Pcurve {
        geometry: source.geometry.clone(),
        t_start: lo,
        t_end: hi,
    };
    (0..=4001)
        .map(|k| {
            let t = lo + (hi - lo) * k as f64 / 4001.0;
            fit.point_at(t).distance_to(&ascending.point_at(t))
        })
        .fold(0.0, f64::max)
}

fn cases() -> Vec<(&'static str, Pcurve)> {
    vec![
        (
            "quarter arc r3",
            Pcurve::elliptical(
                (20.0, -27.0),
                3.0,
                3.0,
                0.0,
                PI,
                3.0 * FRAC_PI_2,
                PI,
                1.5 * PI,
            ),
        ),
        (
            "half arc r7",
            Pcurve::elliptical((0.0, 0.0), 7.0, 7.0, 0.3, 0.0, PI, 0.0, 1.0),
        ),
        (
            "full circle r20",
            Pcurve::elliptical((-20.0, 5.0), 20.0, 20.0, 0.0, 0.0, TAU, 0.0, TAU),
        ),
        (
            "descending ellipse",
            Pcurve::elliptical((1.0, 2.0), 9.0, 4.0, 0.7, 0.2, 4.0, 3.0, -1.0),
        ),
        (
            "sinusoid",
            Pcurve::sinusoidal(0.4, 12.0, 5.0, -3.0, 0.0, TAU),
        ),
    ]
}

#[test]
fn analytic_pcurves_convert_within_the_requested_tolerance() {
    for (name, pcurve) in cases() {
        for tolerance in [1e-6, 1e-9] {
            let (fit, reported) = pcurve.to_nurbs_within(tolerance);
            let measured = dense_deviation(&fit, &pcurve);
            assert!(
                reported <= tolerance,
                "{name} @ {tolerance:e}: reported {reported:e}"
            );
            // The dense check may exceed the probes' estimate only by a
            // small factor (the probes sit mid-span, where the error of
            // a uniform interpolant peaks).
            assert!(
                measured <= 2.0 * tolerance,
                "{name} @ {tolerance:e}: measured {measured:e}"
            );
        }
    }
}

#[test]
fn the_coarse_cubic_misses_the_bound_the_adaptive_fit_meets() {
    // Negative control: the former export (`to_nurbs`, 17-sample cubic)
    // is far outside 1e-9 on the arcs that opened re-import seam gaps.
    for (name, pcurve) in cases().into_iter().take(3) {
        let coarse = dense_deviation(&pcurve.to_nurbs(), &pcurve);
        assert!(
            coarse > 1e-7,
            "{name}: coarse cubic unexpectedly tight ({coarse:e})"
        );
    }
}

#[test]
fn exact_forms_convert_without_deviation() {
    let line = Pcurve::linear((0.0, 1.0), (2.0, 3.0), 0.0, 1.0);
    let (fit, deviation) = line.to_nurbs_within(1e-12);
    assert_eq!(deviation, 0.0);
    assert!(dense_deviation(&fit, &line) < 1e-14);
}
