//! The ring test of `contact_classify`: the gap between two tangent
//! surfaces sampled on a circle about the contact, in its tangent plane.

use super::{ContactSide, SurfaceContactKind, regular};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, PortableMath, Vec3};

/// Samples on the ring.
const RING_SAMPLES: usize = 96;

/// The gap's pattern on one ring.
pub(super) struct RingPattern {
    /// Per sample: +1, −1, or 0 (within the tolerance).
    signs: Vec<i8>,
    pub(super) sign_changes: usize,
}

impl RingPattern {
    /// What the pattern says, for comparing two radii.
    pub(super) fn summary(&self) -> (usize, usize, Option<i8>) {
        let zero_runs = runs(&self.signs, |s| s == 0);
        let nonzero: Vec<i8> = self.signs.iter().copied().filter(|&s| s != 0).collect();
        let one_sign = nonzero
            .first()
            .copied()
            .filter(|&f| nonzero.iter().all(|&s| s == f));
        (self.sign_changes, zero_runs.len(), one_sign)
    }

    pub(super) fn kind(&self, e1: Vec3, e2: Vec3) -> SurfaceContactKind {
        let side_of = |s: i8| {
            if s > 0 {
                ContactSide::AgainstNormal
            } else {
                ContactSide::AlongNormal
            }
        };
        let (changes, zero_runs, one_sign) = self.summary();
        let n = self.signs.len();
        let angle = |k: f64| std::f64::consts::TAU * k / n as f64;
        let dir = |a: f64| e1 * a.pcos() + e2 * a.psin();
        if self.signs.iter().all(|&s| s == 0) {
            return SurfaceContactKind::Coincident;
        }
        if changes == 0 && zero_runs == 0 {
            return SurfaceContactKind::TangentPoint {
                side: side_of(one_sign.unwrap_or(1)),
            };
        }
        if changes == 0 && zero_runs == 2 {
            // A tangent curve through the point: the two zero arcs are
            // where it crosses the ring.
            let mids = runs(&self.signs, |s| s == 0);
            let a = angle(mids[0]);
            return SurfaceContactKind::TangentCurve {
                direction: dir(a),
                side: one_sign.map(side_of),
            };
        }
        if changes > 0 && changes % 2 == 0 {
            // Each branch crosses the ring twice: the changes' angles,
            // the first half of them, are the branch directions.
            let at = change_angles(&self.signs);
            let branches = at
                .iter()
                .take(changes / 2)
                .map(|&k| dir(angle(k)))
                .collect();
            return SurfaceContactKind::CrossingPoint { branches };
        }
        SurfaceContactKind::Undecided {
            reason: format!(
                "the ring's gap changes sign {changes} times over {zero_runs} zero arcs"
            ),
        }
    }
}

/// The midpoints (as fractional sample indices) of the cyclic runs of
/// samples matching `pred`.
fn runs(signs: &[i8], pred: impl Fn(i8) -> bool) -> Vec<f64> {
    let n = signs.len();
    if n == 0 || signs.iter().all(|&s| pred(s)) {
        return Vec::new();
    }
    // Start just after a sample that does not match, so no run wraps.
    let start = (0..n).find(|&i| !pred(signs[i])).unwrap_or(0);
    let mut out = Vec::new();
    let mut k = 1;
    while k <= n {
        let i = (start + k) % n;
        if pred(signs[i]) {
            let first = start + k;
            while k <= n && pred(signs[(start + k) % n]) {
                k += 1;
            }
            let last = start + k - 1;
            out.push(((first + last) as f64 * 0.5) % n as f64);
        } else {
            k += 1;
        }
    }
    out
}

/// Where the nonzero sign flips around the ring, as fractional sample
/// indices (the middle of the zero run between, if any).
fn change_angles(signs: &[i8]) -> Vec<f64> {
    let n = signs.len();
    let nz: Vec<usize> = (0..n).filter(|&i| signs[i] != 0).collect();
    let mut out = Vec::new();
    for (k, &i) in nz.iter().enumerate() {
        let j = nz[(k + 1) % nz.len()];
        if signs[i] != signs[j] {
            let gap = (j + n - i) % n;
            out.push((i as f64 + 0.5 * gap as f64) % n as f64);
        }
    }
    out.sort_by(f64::total_cmp);
    out
}

/// The gap `h₁ − h₂` on the ring of radius `rho` about `p` in the
/// tangent plane, as signs against the tolerance; `None` if a height
/// does not converge.
#[allow(clippy::too_many_arguments)]
pub(super) fn ring_pattern(
    s1: &SurfaceGeometry,
    uv1: (f64, f64),
    s2: &SurfaceGeometry,
    uv2: (f64, f64),
    (n, e1, e2): (Vec3, Vec3, Vec3),
    p: Point3,
    rho: f64,
    linear: f64,
) -> Option<RingPattern> {
    let mut signs = Vec::with_capacity(RING_SAMPLES);
    let (mut w1, mut w2) = (off_pole(s1, uv1), off_pole(s2, uv2));
    for k in 0..RING_SAMPLES {
        let a = std::f64::consts::TAU * k as f64 / RING_SAMPLES as f64;
        let q = p + (e1 * a.pcos() + e2 * a.psin()) * rho;
        // Each height from the surface's own closest point to the
        // sample (robust where the chart is not: beside a pole), else
        // from the last sample's.
        let start1 = s1
            .closest_point(&q)
            .filter(|&c| regular(s1, c))
            .unwrap_or(w1);
        let start2 = s2
            .closest_point(&q)
            .filter(|&c| regular(s2, c))
            .unwrap_or(w2);
        let (h1, u1) = height(s1, start1, q, (n, e1, e2))?;
        let (h2, u2) = height(s2, start2, q, (n, e1, e2))?;
        (w1, w2) = (u1, u2);
        let g = h1 - h2;
        signs.push(if g.abs() <= linear {
            0
        } else if g > 0.0 {
            1
        } else {
            -1
        });
    }
    let nz: Vec<i8> = signs.iter().copied().filter(|&s| s != 0).collect();
    let sign_changes = (0..nz.len())
        .filter(|&i| nz[i] != nz[(i + 1) % nz.len()])
        .count();
    Some(RingPattern {
        signs,
        sign_changes,
    })
}

/// A parameter beside `uv` where the chart is regular (the pole itself
/// has no Jacobian for Newton).
fn off_pole(s: &SurfaceGeometry, uv: (f64, f64)) -> (f64, f64) {
    if regular(s, uv) {
        return uv;
    }
    let ((_, _), (v0, v1)) = s.parameter_range();
    let step = if v0.is_finite() && v1.is_finite() {
        1e-3 * (v1 - v0)
    } else {
        1e-3
    };
    let v = if uv.1 - step >= v0 {
        uv.1 - step
    } else {
        uv.1 + step
    };
    (uv.0, v)
}

/// The height along `n` at which the line through `q` along `n` meets
/// `s`, by Newton from `start`, and the parameter there.
fn height(
    s: &SurfaceGeometry,
    start: (f64, f64),
    q: Point3,
    (n, e1, e2): (Vec3, Vec3, Vec3),
) -> Option<(f64, (f64, f64))> {
    let (mut u, mut v) = start;
    let scale = (q - s.point_at(u, v)).length().max(1e-300);
    for _ in 0..40 {
        let d = s.point_at(u, v) - q;
        let (f1, f2) = (d.dot(&e1), d.dot(&e2));
        if f1.abs().max(f2.abs()) <= 1e-13 * scale.max(1.0) {
            return Some((d.dot(&n), (u, v)));
        }
        let su = s.partial_u(u, v);
        let sv = s.partial_v(u, v);
        let (a, b, c, dd) = (su.dot(&e1), sv.dot(&e1), su.dot(&e2), sv.dot(&e2));
        let det = a * dd - b * c;
        if det.abs() <= 1e-300 {
            return None;
        }
        u -= (dd * f1 - b * f2) / det;
        v -= (-c * f1 + a * f2) / det;
        if !(u.is_finite() && v.is_finite()) {
            return None;
        }
    }
    let d = s.point_at(u, v) - q;
    (d.dot(&e1).abs().max(d.dot(&e2).abs()) <= 1e-9 * scale.max(1.0)).then(|| (d.dot(&n), (u, v)))
}
