//! A co-parameterized C1 cubic chart pcurve, for curves the degree-1
//! chart fit cannot hold within the linear tolerance.
//!
//! The degree-1 fits (`chart_pcurve::fit_chart_pcurve`, the sphere arms'
//! `fit_sphere_uv_curve`) double their samples until every mid-knot
//! deviation passes the tolerance, capped at 4097 points. A chord's sag
//! falls only with the square of its step, and the cap is reached on
//! large parts: two balls of radius 1000 crossing in a circle 600 across
//! kept a pcurve whose image strayed 7.4e-4 from the exact circle, past
//! the model's 1e-4 imprecision, and every boolean declined (generated
//! corpus SS-X at scale 1000). Here each sample's chart velocity (a
//! central difference of the chart along the curve) makes a cubic
//! Hermite span between samples, whose error falls with the fourth power
//! of the step; knots stay at the sample parameters, so the pcurve keeps
//! the curve's own parameterization.

use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;

/// The first sample count tried, doubled (less one) up to the cap.
const FIRST_SAMPLES: usize = 17;
/// The most samples a fit takes.
const CAP_SAMPLES: usize = 4097;

/// A chart value on the period branch nearest a previous one.
type ChartAt<'a> = &'a dyn Fn(f64, Option<(f64, f64)>) -> Option<(f64, f64)>;

/// A curve to chart: its point at a parameter, the chart value of a
/// parameter on the period branch nearest a previous value, and the
/// surface point at a chart value.
pub(crate) struct ChartedCurve<'a> {
    /// The 3-D curve's point at `t`.
    pub(crate) curve: &'a dyn Fn(f64) -> Point3,
    /// The chart value at `t`, unwrapped next to the given value.
    pub(crate) chart: ChartAt<'a>,
    /// The surface point at a chart value.
    pub(crate) surface: &'a dyn Fn(f64, f64) -> Point3,
}

/// The C1 cubic chart pcurve of `c` over `range` whose image passes
/// within `gate` of the curve at every span's quarter points, or the
/// densest tried when the cap is reached first, kept only where its
/// worst gap is under `beat` (the caller's capped degree-1 fit's: beside
/// a singular row the chart velocity runs away and the cubic is the
/// worse fit, a ball's seam meridian to its pole); `None` otherwise, or
/// where a chart value or velocity is not finite.
pub(crate) fn fit_hermite_chart_pcurve(
    c: &ChartedCurve<'_>,
    (t0, t1): (f64, f64),
    gate: f64,
    beat: f64,
) -> Option<NurbsCurve> {
    let delta = (t1 - t0).abs() * 1.0e-7;
    let mut n = FIRST_SAMPLES;
    loop {
        let mut params = Vec::with_capacity(n);
        let mut uvs: Vec<(f64, f64)> = Vec::with_capacity(n);
        let mut ders: Vec<(f64, f64)> = Vec::with_capacity(n);
        for i in 0..n {
            let t = t0 + (t1 - t0) * i as f64 / (n - 1) as f64;
            let uv = (c.chart)(t, uvs.last().copied())?;
            let (a, b) = ((t - delta).max(t0.min(t1)), (t + delta).min(t0.max(t1)));
            let (ua, ub) = ((c.chart)(a, Some(uv))?, (c.chart)(b, Some(uv))?);
            let d = ((ub.0 - ua.0) / (b - a), (ub.1 - ua.1) / (b - a));
            if ![uv.0, uv.1, d.0, d.1].iter().all(|x| x.is_finite()) {
                return None;
            }
            params.push(t);
            uvs.push(uv);
            ders.push(d);
        }
        let fit = hermite(&params, &uvs, &ders);
        let worst = params
            .windows(2)
            .flat_map(|w| [0.25, 0.5, 0.75].map(|f| w[0] + (w[1] - w[0]) * f))
            .map(|t| {
                use crate::geometry::curve::CurveEval as _;
                let q = fit.point_at(t);
                (c.surface)(q.x, q.y).distance_to(&(c.curve)(t))
            })
            .fold(0.0_f64, f64::max);
        if worst <= gate || n >= CAP_SAMPLES {
            return (worst < beat).then_some(fit);
        }
        n = (n - 1) * 2 + 1;
    }
}

/// The C1 piecewise cubic through `uvs` at `params` with velocities
/// `ders`, in Bézier form: each interior knot triple.
fn hermite(params: &[f64], uvs: &[(f64, f64)], ders: &[(f64, f64)]) -> NurbsCurve {
    let n = params.len();
    let p = |(u, v): (f64, f64)| Point3::new(u, v, 0.0);
    let mut control_points = Vec::with_capacity(3 * (n - 1) + 1);
    let mut knots = vec![params[0]; 4];
    for i in 0..n - 1 {
        let h = (params[i + 1] - params[i]) / 3.0;
        control_points.push(p(uvs[i]));
        control_points.push(p((uvs[i].0 + h * ders[i].0, uvs[i].1 + h * ders[i].1)));
        control_points.push(p((
            uvs[i + 1].0 - h * ders[i + 1].0,
            uvs[i + 1].1 - h * ders[i + 1].1,
        )));
        let k = if i + 1 == n - 1 { 4 } else { 3 };
        knots.extend(std::iter::repeat_n(params[i + 1], k));
    }
    control_points.push(p(uvs[n - 1]));
    let weights = vec![1.0; control_points.len()];
    NurbsCurve {
        degree: 3,
        knots,
        control_points,
        weights,
    }
}
