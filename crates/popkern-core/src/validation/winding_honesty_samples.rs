//! Polygon samples of loop pcurves for the winding audit (split from
//! `winding_honesty.rs` under the 500-line file budget).

/// Uniform sample count per half-edge pcurve when polygonizing the
/// loop.
const SAMPLES_PER_HALF_EDGE: usize = 8;

/// Ascending fractions of `(t0, t1)` to polygonize `pcurve` at:
/// `SAMPLES_PER_HALF_EDGE` uniform ones, and for a NURBS pcurve every
/// knot inside the range, so each polynomial piece contributes its
/// ends. A piece lies in its control points' hull, so the polygon
/// follows any turn the pieces resolve; eight samples of the whole
/// curve need not (PK-1055: a hole passing 0.05 from a sphere's poles,
/// whose pcurve sweeps half a turn of u at each, polygonized as an
/// octagon of area +3.8 for its true -13.4, and its face declined as
/// wound like its outer loop).
pub(super) fn sample_fractions(
    pcurve: &crate::geometry::pcurve::Pcurve,
    (t0, t1): (f64, f64),
) -> Vec<f64> {
    let mut fractions: Vec<f64> = (0..SAMPLES_PER_HALF_EDGE)
        .map(|k| k as f64 / (SAMPLES_PER_HALF_EDGE - 1) as f64)
        .collect();
    if let Some(n) = pcurve.as_nurbs()
        && t1 != t0
    {
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        fractions.extend(
            n.knots
                .iter()
                .filter(|&&k| k > lo && k < hi)
                .map(|&k| (k - t0) / (t1 - t0)),
        );
        fractions.sort_by(f64::total_cmp);
        fractions.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
    }
    fractions
}

/// Append a UV sample, shifting by whole periods so consecutive
/// samples stay chart-continuous across periodic seams.
pub(super) fn push_unwrapped(
    polygon: &mut Vec<(f64, f64)>,
    uv: (f64, f64),
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    let Some(&(prev_u, prev_v)) = polygon.last() else {
        polygon.push(uv);
        return;
    };
    let unwrap = |value: f64, prev: f64, period: Option<f64>| -> f64 {
        let Some(p) = period else {
            return value;
        };
        if p <= 0.0 {
            return value;
        }
        let shift = ((prev - value) / p).round();
        value + shift * p
    };
    polygon.push((
        unwrap(uv.0, prev_u, period_u),
        unwrap(uv.1, prev_v, period_v),
    ));
}
