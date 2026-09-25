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
//! Every ring *gate* measures this sine. It is orientation-free: a fold
//! across the spoke (normals 180° apart) would read 0. So a gate reads
//! a sample whose relative orientation disagrees with the spoke's sense
//! as 1. The sense is the majority on the promotion gate, and the
//! patches' reference normals ([`spoke_senses`]) in the repair and the
//! optimizers, since a candidate can have almost a whole spoke
//! reversed.
//!
//! The honest objective is the oriented chord [`spoke_normal_chord`],
//! `2 sin(θ/2)` (PK-544). It is skew-invariant like the sine, and
//! monotone over the whole half-turn, so a minimizer is pushed out of a
//! fold, never into one. Around it (PK-547):
//! - a chart-lean barrier ([`lean_penalty`]) supplies what skew
//!   invariance leaves out;
//! - [`whole_spoke_parameters`] carries every honest measure into the
//!   spoke end zones.
//!
//! The normalized triple product remains only as a search heuristic,
//! whose candidates are gated and repaired honestly.
//!
//! A chart whose cross tangent is parallel to the spoke has no normal
//! there. The denominators carry a relative floor
//! [`CHART_ANGLE_FLOOR`] so the measure stays finite and reads large
//! (never small) as the chart degenerates.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

/// The commercial ring G1 gate on the spoke normal sine (≈ 2.9°): the
/// K-P promotion's ring gate and the bar the ring optimizers aim for.
pub(crate) const RING_SPOKE_G1_GATE: f64 = 5.0e-2;

/// The promotion window's 50 spoke parameters on `[0.05, 0.95]`: the
/// phase-2/3 search heuristics' sampling.
pub(crate) fn window_spoke_parameters() -> Vec<f64> {
    (0..50).map(|k| 0.05 + 0.90 * (k as f64) / 49.0).collect()
}

/// The whole-spoke parameters every honest ring measure samples (PK-547):
/// the promotion window plus both end zones, clustered toward each end
/// down to 1e-4. The window alone hid a 5%-long crease at the rim end
/// of the served extreme-ladder rings: spoke sine 4.4e-2 at t = 0.95,
/// 0.25 at 0.99 and 0.325 at 1, with a healthy chart. The exact corner
/// points are left out, since there the rational-Gregory evaluators
/// take conventions rather than limits.
pub(crate) fn whole_spoke_parameters() -> Vec<f64> {
    let mut out = window_spoke_parameters();
    for &d in &[1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 3e-2, 4e-2] {
        out.push(d);
        out.push(1.0 - d);
    }
    out
}

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

/// The oriented spoke normal chord `|n̂_a − σ n̂_b| = 2 sin(θ/2)`
/// (PK-544). Here `n̂_a`, `n̂_b` are the unit normals `w × u` and
/// `w × v`, `σ = ±1` is the spoke's sense (`true` = `+1`), and `θ` is
/// the angle between `n̂_a` and `σ n̂_b`.
///
/// Unlike the sine, it is monotone over the whole `θ ∈ [0, π]`: 0 only
/// at true G1, 2 at a fold across the spoke. An optimizer minimizing it
/// is pushed out of a fold, never into one. It shares the sine's
/// invariance to chart skew, since only the normal directions enter.
/// The chart-angle floor biases it by at most `√2 · CHART_ANGLE_FLOOR`
/// at G1. Reads 0 when every vector vanishes, like
/// [`spoke_normal_sine`].
pub(crate) fn spoke_normal_chord(u: Vec3, v: Vec3, w: Vec3, sense: bool) -> f64 {
    let (na, nb) = (floored_cross_norm(w, u), floored_cross_norm(w, v));
    if na * nb <= 1e-300 {
        return 0.0;
    }
    let sign = if sense { 1.0 } else { -1.0 };
    let c = sign * w.cross(&u).dot(&w.cross(&v)) / (na * nb);
    (2.0 * (1.0 - c)).max(0.0).sqrt()
}

/// Gradient of [`spoke_normal_chord`] with respect to `u`, `v`, `w`.
/// With `A = w × u`, `B = w × v` and `c = σ A·B / (|A||B|)`:
/// - `∂c/∂A = σ B/(|A||B|) − c A/|A|²`, and `∂c/∂B` likewise;
/// - `∂c/∂u = ∂c/∂A × w`, `∂c/∂v = ∂c/∂B × w`, and
///   `∂c/∂w = u × ∂c/∂A + v × ∂c/∂B`;
/// - `f = √(2 − 2c)`, so `∂f = −∂c / f`.
///
/// `None` at the chord's non-differentiable zero (`f` below 1e-12) or
/// on vanishing vectors. The chart-angle floor's derivative is
/// dropped, as for the sine.
pub(crate) fn spoke_normal_chord_gradient(
    u: Vec3,
    v: Vec3,
    w: Vec3,
    sense: bool,
) -> Option<(Vec3, Vec3, Vec3)> {
    let (na, nb) = (floored_cross_norm(w, u), floored_cross_norm(w, v));
    if na * nb <= 1e-300 {
        return None;
    }
    let sign = if sense { 1.0 } else { -1.0 };
    let (a, b) = (w.cross(&u), w.cross(&v));
    let c = sign * a.dot(&b) / (na * nb);
    let f = (2.0 * (1.0 - c)).max(0.0).sqrt();
    if f < 1e-12 {
        return None;
    }
    let ga = b * (sign / (na * nb)) - a * (c / (na * na));
    let gb = a * (sign / (na * nb)) - b * (c / (nb * nb));
    let scale = -1.0 / f;
    Some((
        ga.cross(&w) * scale,
        gb.cross(&w) * scale,
        (u.cross(&ga) + v.cross(&gb)) * scale,
    ))
}

/// Chart-angle sine `|w × u| / (|w| |u|)` below which the lean
/// barrier engages (≈ 5.7°, PK-547). The chord is invariant to chart
/// skew, so an optimizer minimizing it alone is free to lean a cross
/// tangent onto the spoke. The triple-product heuristic actively
/// rewards that lean: on pyramid3_uniform its ring reached a chart
/// sine of 2.6e-3, and the served chart's Jacobian dipped 70×
/// mid-spoke, leaving half a seam uncertifiable.
pub(crate) const LEAN_FLOOR: f64 = 0.1;

/// Hinge penalty on a cross tangent `u` leaning onto the spoke `w`
/// past [`LEAN_FLOOR`]: 0 above it, rising linearly to 1 at a
/// degenerate chart.
pub(crate) fn lean_penalty(u: Vec3, w: Vec3) -> f64 {
    let scale = u.length() * w.length();
    if scale.is_nan() || scale <= 0.0 {
        return 1.0;
    }
    ((LEAN_FLOOR - u.cross(&w).length() / scale) / LEAN_FLOOR).max(0.0)
}

/// Gradient of [`lean_penalty`] with respect to `u` and `w`. With
/// `c = u × w` and `s = |c| / (|u||w|)`:
/// - `∂s/∂u = (w × c)/(|c||u||w|) − s u/|u|²`;
/// - `∂s/∂w = (c × u)/(|c||u||w|) − s w/|w|²`;
/// - the penalty's gradient is `−∂s / LEAN_FLOOR`.
///
/// `(0, 0)` when inactive, `None` on a degenerate configuration.
pub(crate) fn lean_penalty_gradient(u: Vec3, w: Vec3) -> Option<(Vec3, Vec3)> {
    let (nu, nw) = (u.length(), w.length());
    let c = u.cross(&w);
    let nc = c.length();
    let product = nu * nw;
    if product.is_nan() || product <= 0.0 || nc.is_nan() || nc <= 0.0 {
        return None;
    }
    let s = nc / (nu * nw);
    if s >= LEAN_FLOOR {
        return Some((Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)));
    }
    let scale = -1.0 / LEAN_FLOOR;
    let ds_du = w.cross(&c) * (1.0 / (nc * nu * nw)) - u * (s / (nu * nu));
    let ds_dw = c.cross(&u) * (1.0 / (nc * nu * nw)) - w * (s / (nw * nw));
    Some((ds_du * scale, ds_dw * scale))
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

/// Expected relative chart orientation along each spoke `k`: patch
/// `k`'s `u = 1` side against patch `k + 1`'s `v = 0` side. The
/// spoke's `(w × u) · (w × v)` equals `n_a · n_b`, the dot product of
/// the two chart normals there (PK-544). An orientation-consistent
/// ring therefore carries the sign of the two patches' reference
/// (mid-face) normals' dot product along the whole spoke.
///
/// A per-spoke majority vote is no substitute: the exact-injection
/// candidate on cube_extreme 1:20 has 49 of 50 samples reversed on
/// one spoke.
pub(crate) fn spoke_senses(references: &[(Vec3, f64)]) -> Vec<bool> {
    let n = references.len();
    (0..n)
        .map(|k| references[k].0.dot(&references[(k + 1) % n].0) >= 0.0)
        .collect()
}

/// Spoke samples whose relative orientation opposes their spoke's
/// sense: folds across the spoke, where the normals stand more than
/// 90° apart, over [`whole_spoke_parameters`].
/// Inside a fold-free region the spoke normal sine is monotone in the
/// normal angle, so a polish constrained here optimizes it honestly.
pub(crate) fn spoke_fold_count(ring: &[RectangularGregoryBiNonic], senses: &[bool]) -> usize {
    let parameters = whole_spoke_parameters();
    let n = ring.len();
    let mut folds = 0usize;
    for (k, &sense) in senses.iter().enumerate().take(n) {
        let (a, b) = (&ring[k], &ring[(k + 1) % n]);
        for &t in &parameters {
            let (u, v, w) = (
                a.partial_u(1.0, t),
                b.partial_v(1.0 - t, 0.0),
                a.partial_v(1.0, t),
            );
            if spoke_orientation(u, v, w) != sense {
                folds += 1;
            }
        }
    }
    folds
}

/// [`spoke_senses`] of a ring from its own mid-face chart normals.
pub(crate) fn ring_spoke_senses(ring: &[RectangularGregoryBiNonic]) -> Vec<bool> {
    let references: Vec<(Vec3, f64)> = ring
        .iter()
        .map(|p| (p.partial_u(0.5, 0.5).cross(&p.partial_v(0.5, 0.5)), 1.0))
        .collect();
    spoke_senses(&references)
}

/// Worst spoke G1 of a ring against fixed spoke senses: the normal
/// sine per sample, 1 where the sample folds against its spoke's sense.
/// This is [`worst_spoke_measure`] with the orientation taken from the
/// patches' reference normals rather than a per-spoke majority.
pub(crate) fn worst_oriented_ring_sine(ring: &[RectangularGregoryBiNonic], senses: &[bool]) -> f64 {
    let parameters = whole_spoke_parameters();
    let n = ring.len();
    let mut worst = 0.0_f64;
    for (k, &sense) in senses.iter().enumerate().take(n) {
        let (a, b) = (&ring[k], &ring[(k + 1) % n]);
        for &t in &parameters {
            let (u, v, w) = (
                a.partial_u(1.0, t),
                b.partial_v(1.0 - t, 0.0),
                a.partial_v(1.0, t),
            );
            let measure = if spoke_orientation(u, v, w) == sense {
                spoke_normal_sine(u, v, w)
            } else {
                1.0
            };
            worst = worst.max(measure);
        }
    }
    worst
}

/// Sum of squared oriented spoke normal chords over a ring's spoke
/// samples: the honest ring residual the K-P optimizers compare
/// candidates by (PK-544).
pub(crate) fn ring_chord_residual(ring: &[RectangularGregoryBiNonic], senses: &[bool]) -> f64 {
    let parameters = whole_spoke_parameters();
    let n = ring.len();
    let mut sum = 0.0;
    for (k, &sense) in senses.iter().enumerate().take(n) {
        let (a, b) = (&ring[k], &ring[(k + 1) % n]);
        for &t in &parameters {
            let chord = spoke_normal_chord(
                a.partial_u(1.0, t),
                b.partial_v(1.0 - t, 0.0),
                a.partial_v(1.0, t),
                sense,
            );
            sum += chord * chord;
        }
    }
    sum
}

#[cfg(test)]
#[path = "kp_spoke_normal_sine_tests.rs"]
mod tests;
