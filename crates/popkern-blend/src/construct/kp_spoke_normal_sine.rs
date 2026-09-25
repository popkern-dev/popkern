//! Spoke G1 residual of the K-P bi-nonic ring: the sine of the angle
//! between the two patches' normals along a shared spoke (PK-539).
//!
//! With `u`, `v` the two patches' cross tangents at a spoke sample
//! and `w` the spoke tangent, the patch normals are `w × u` and
//! `w × v`, and by `(w×u)×(w×v) = w (w·(u×v))`
//!
//! `sine = |w| |T| / (|w×u| |w×v|)`, `T = (u × v) · w`.
//!
//! The ring used the normalized triple product
//! `|T| / (|u| |v| |w|) = sine · sin∠(w,u) · sin∠(w,v)` as both its
//! optimizer objective and its promotion gate. That product reads
//! small whenever a cross tangent leans toward the spoke, so the
//! optimizer reduced it by skewing the chart instead of closing the
//! crease: measured on cube_extreme (0.02/0.10/0.20), a spoke sample
//! at 1.7° between cross tangent and spoke read 6.9e-4 while the
//! normals stood 3.9° apart (sine 6.7e-2), and the served ring claimed
//! G1 at a worst sine of 1.3e-1. The sine is invariant to that skew.
//!
//! Every ring *gate* now measures this sine. The optimizers keep the
//! triple product as their search objective: optimizing the sine
//! itself was measured to reach true G1 residuals (4.6e-3 sum of
//! squares) only with candidates whose neighbouring patches nearly
//! cross, and each such candidate cost ~1 minute in the ring-crossing
//! SSI gate, at a cost G2 production requests would pay too. The
//! honest objective with crossing-aware constraints is a filed
//! follow-up; until then a creased candidate is declined at the gate,
//! never served as G1.
//!
//! The sine alone is orientation-free: a fold across the spoke
//! (normals 180° apart) would read 0. The ring measure reads a sample
//! whose relative orientation disagrees with the rest of its spoke as
//! 1 (a fold); a chart mirrored along the whole spoke is legitimate.
//!
//! A chart whose cross tangent is parallel to the spoke has no normal
//! there. The denominators carry a relative floor
//! [`CHART_ANGLE_FLOOR`] so the measure stays finite and reads large
//! (never small) as the chart degenerates.

use popkern_core::geometry::Vec3;

/// Relative floor on `|w × u| / (|w| |u|)` (the chart angle's sine):
/// below it the chart counts as degenerate.
pub(crate) const CHART_ANGLE_FLOOR: f64 = 1e-6;

/// `|a × b|` with the relative floor `CHART_ANGLE_FLOOR · |a| |b|`
/// folded in.
fn floored_cross_norm(a: Vec3, b: Vec3) -> f64 {
    let floor = CHART_ANGLE_FLOOR * a.length() * b.length();
    (a.cross(&b).length_squared() + floor * floor).sqrt()
}

/// The spoke normal sine for cross tangents `u`, `v` and spoke
/// tangent `w` (0 when every vector vanishes). Orientation-free: see
/// [`spoke_orientation`] for folds.
pub(crate) fn spoke_normal_sine(u: Vec3, v: Vec3, w: Vec3) -> f64 {
    let t = u.cross(&v).dot(&w);
    let (nu, nv) = (floored_cross_norm(w, u), floored_cross_norm(w, v));
    let denom = nu * nv;
    if denom > 1e-300 {
        w.length() * t.abs() / denom
    } else {
        0.0
    }
}

/// The spoke normal sine and its gradients with respect to `u`, `v`,
/// `w` (orientation-free, as [`spoke_normal_sine`]). With
/// `s = sign(T) |w| / (|w×u| |w×v|)` and `f = s T`:
/// - `∂f/∂u = s (v×w) − f ((w×u)×w) / |w×u|²`
/// - `∂f/∂v = s (w×u) − f ((w×v)×w) / |w×v|²`
/// - `∂f/∂w = f w/|w|² + s (u×v) − f (u×(w×u))/|w×u|² − f (v×(w×v))/|w×v|²`
///
/// (the chart-angle floor's own derivative is dropped: it only acts
/// on degenerate charts). `None` when every vector vanishes.
pub(crate) fn spoke_normal_sine_gradient(u: Vec3, v: Vec3, w: Vec3) -> Option<(Vec3, Vec3, Vec3)> {
    let t = u.cross(&v).dot(&w);
    let (nu, nv) = (floored_cross_norm(w, u), floored_cross_norm(w, v));
    let w_sq = w.length_squared();
    if nu * nv <= 1e-300 || w_sq <= 1e-300 {
        return None;
    }
    let (nu_sq, nv_sq) = (nu * nu, nv * nv);
    let sign = if t >= 0.0 { 1.0 } else { -1.0 };
    let s = sign * w.length() / (nu * nv);
    let f = s * t;
    let (wu, wv) = (w.cross(&u), w.cross(&v));
    let du = v.cross(&w) * s - wu.cross(&w) * (f / nu_sq);
    let dv = wu * s - wv.cross(&w) * (f / nv_sq);
    let dw =
        w * (f / w_sq) + u.cross(&v) * s - u.cross(&wu) * (f / nu_sq) - v.cross(&wv) * (f / nv_sq);
    Some((du, dv, dw))
}

/// Relative orientation of the two sides' normals at a spoke sample:
/// `true` when `(w×u)·(w×v) ≥ 0`. A ring's charts may legitimately
/// be mirrored against each other (the face orientation carries
/// that), so an absolute sign means nothing; a sample whose
/// orientation disagrees with the rest of its spoke is a fold.
pub(crate) fn spoke_orientation(u: Vec3, v: Vec3, w: Vec3) -> bool {
    w.cross(&u).dot(&w.cross(&v)) >= 0.0
}

/// Worst spoke G1 over one spoke's samples `(u, v, w)`: the normal
/// sine where the sample agrees with the spoke's majority orientation,
/// 1 where it disagrees (a fold across the spoke, which the
/// orientation-free sine would read as G1).
pub(crate) fn worst_spoke_measure(samples: &[(Vec3, Vec3, Vec3)]) -> f64 {
    let agreeing = samples
        .iter()
        .filter(|(u, v, w)| spoke_orientation(*u, *v, *w))
        .count();
    let majority = 2 * agreeing >= samples.len();
    samples
        .iter()
        .map(|&(u, v, w)| {
            if spoke_orientation(u, v, w) == majority {
                spoke_normal_sine(u, v, w)
            } else {
                1.0
            }
        })
        .fold(0.0, f64::max)
}

#[cfg(test)]
#[path = "kp_spoke_normal_sine_tests.rs"]
mod tests;
