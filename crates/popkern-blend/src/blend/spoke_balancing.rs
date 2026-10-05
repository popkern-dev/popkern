//! Balanced spoke tangents at the ring center (PK-240).
//!
//! # The chart-vs-geometry mismatch
//!
//! The N-quad decomposition fans the 2N-gon boundary ring from its
//! central point `C` with spokes aimed straight at alternating hex
//! vertices. The KP guided promotion rungs (binonic, biquintic,
//! c-twist) parameterize that fan through the REGULAR characteristic
//! map — their τ/µ tables are valence-parameterized with
//! `c = cos(2π/N)`, i.e. every sector is assumed to occupy a 2π/N
//! wedge around `C`. At extreme radius asymmetry the hex vertices
//! cluster angularly (measured cube_extreme_20, radii 1:20: sector
//! spans 45.7° / 140.7° / 173.6° at `C`, tiny arc side 0.014 vs big
//! sides ~0.24 — while `C` stays roughly equidistant, 0.14–0.19,
//! from every vertex). With geometry at 173° against a 120° chart,
//! the guide-composition rows distort hardest along the flat
//! sector's two spoke-to-center bands — exactly where the fold
//! slivers were localized — and every optimizer arm's product
//! self-intersects (candidate=guide), so the binonic promotion
//! declines and the cap ships creased plain Gregory (g1 1.9e-1 at
//! 1:20 vs the binonic class 1e-3).
//!
//! # The fix
//!
//! The spoke's direction AT `C` is a free design parameter — only
//! its endpoints and the ring boundary are constrained. When the
//! natural sector split is degenerate, give the N spoke tangents at
//! `C` an equal 2π/N spacing in the tangent plane (global rotation =
//! circular mean of the natural deviations, cyclic order preserved)
//! and let each spoke CURVE out to its hex vertex as a cubic. Every
//! sector's corner angle at `C` becomes 2π/N — the geometry then
//! matches the regular chart the KP tables assume. Below the
//! engagement threshold the decomposition is byte-identical to the
//! straight-spoke path.

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_math::PortableMath;

/// Natural sector span (radians) above which balancing engages.
/// The guided rung's fold class was measured from 173.6° down to
/// ~140°-span sectors at 1:10 serving fold-free; 150° sits between
/// the served and folding classes with margin on both sides.
const MAX_NATURAL_SECTOR_SPAN: f64 = 150.0 * std::f64::consts::PI / 180.0;

/// Chosen balanced-fan geometry for one ring decomposition.
pub(super) struct BalancedSpokePlan {
    /// Unit spoke tangent directions at `C` (3-D, in the tangent
    /// plane of `n_c`), one per spoke `k` — spoke `k` targets hex
    /// vertex `V_{2(k+1) mod 2N}`.
    pub dirs: Vec<Vec3>,
    /// Worst natural sector span (radians) that triggered balancing.
    pub natural_max_span: f64,
    /// Worst sector span after balancing — `2π/N` at strength 1,
    /// interpolated toward the natural span at lower strengths.
    pub balanced_max_span: f64,
}

/// Measure the natural sector split at `C` and produce a balanced
/// plan when it is degenerate. Returns `None` when the split is
/// already regular enough (the caller keeps straight spokes), when
/// the normal is unusable, or when any spoke direction has no
/// tangent-plane component.
///
/// `strength` in `(0, 1]` interpolates each spoke tangent between
/// its natural direction (0) and the fully regular 2π/N spacing (1).
/// The fold-triggered retry ladder walks it downward: full
/// regularization can over-distort the smallest sector (its outer
/// arc is tiny while its angular span grows to 2π/N — measured
/// pre-solve folds of the initial fan at 1:20/1:30), while a partial
/// balance can already clear the guided rung's fold threshold.
pub(super) fn balanced_spoke_plan(
    hex_vertices: &[Point3],
    c: Point3,
    n_c: Vec3,
    n: usize,
    strength: f64,
) -> Option<BalancedSpokePlan> {
    if !(strength > 0.0 && strength <= 1.0) {
        return None;
    }
    if n < 3 || hex_vertices.len() != 2 * n {
        return None;
    }
    let nn = n_c.normalized()?;
    let t1 = super::surface_math::pick_perpendicular(&nn);
    let t2 = nn.cross(&t1);

    let mut theta = Vec::with_capacity(n);
    for k in 0..n {
        let v = hex_vertices[(2 * k + 2) % (2 * n)];
        let d = super::gregory_parallelogram::sub3(v, c);
        let (x, y) = (d.dot(&t1), d.dot(&t2));
        if (x * x + y * y).sqrt() < 1e-12 {
            return None;
        }
        theta.push(y.patan2(x));
    }

    let tau = std::f64::consts::TAU;
    let natural_max_span = (0..n)
        .map(|k| {
            let mut span = theta[k] - theta[(k + n - 1) % n];
            while span <= 0.0 {
                span += tau;
            }
            while span > tau {
                span -= tau;
            }
            span
        })
        .fold(0.0_f64, f64::max);
    if natural_max_span <= MAX_NATURAL_SECTOR_SPAN {
        return None;
    }

    // Balanced angles ψ_k = φ + k·2π/N with the global rotation φ
    // chosen as the circular mean of the per-spoke deviations —
    // the least-squares rotation on the circle.
    let delta = tau / n as f64;
    let (mut sin_sum, mut cos_sum) = (0.0_f64, 0.0_f64);
    for (k, &th) in theta.iter().enumerate() {
        let dev = th - k as f64 * delta;
        sin_sum += dev.psin();
        cos_sum += dev.pcos();
    }
    let phi = sin_sum.patan2(cos_sum);

    // Interpolate each tangent from its natural angle toward the
    // regular target by `strength`, through the short way around.
    let wrap_to_pi = |mut a: f64| {
        while a > std::f64::consts::PI {
            a -= tau;
        }
        while a < -std::f64::consts::PI {
            a += tau;
        }
        a
    };
    let psi: Vec<f64> = theta
        .iter()
        .enumerate()
        .map(|(k, &th)| {
            let target = phi + k as f64 * delta;
            th + strength * wrap_to_pi(target - th)
        })
        .collect();
    let balanced_max_span = (0..n)
        .map(|k| {
            let mut span = psi[k] - psi[(k + n - 1) % n];
            while span <= 0.0 {
                span += tau;
            }
            while span > tau {
                span -= tau;
            }
            span
        })
        .fold(0.0_f64, f64::max);
    let mut dirs = Vec::with_capacity(n);
    for &angle in &psi {
        dirs.push(t1 * angle.pcos() + t2 * angle.psin());
    }
    Some(BalancedSpokePlan {
        dirs,
        natural_max_span,
        balanced_max_span,
    })
}

/// Cubic Bezier spoke from hex vertex `v` to the center `c`,
/// arriving at `c` travelling along `-dir_at_c` (so the reversed
/// orientation departs `C` along `dir_at_c`). The vertex-end leg
/// follows the straight chord; both legs are one third of the chord
/// length — the straight spoke's own control layout, bent only at
/// the `C` end.
pub(super) fn curved_spoke(v: Point3, c: Point3, dir_at_c: Vec3) -> Option<NurbsCurve> {
    let chord = super::gregory_parallelogram::sub3(c, v);
    let len = chord.length();
    if len < 1e-12 {
        return None;
    }
    let leg = len / 3.0;
    let chord_dir = chord * (1.0 / len);
    let p0 = v;
    let p1 = v + chord_dir * leg;
    let p2 = c + dir_at_c * leg;
    let p3 = c;
    Some(NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![p0, p1, p2, p3],
        weights: vec![1.0; 4],
    })
}

/// Reverse a clamped Bezier spoke (same geometry, opposite
/// orientation). Valid for the uniform clamped knot vector
/// [`curved_spoke`] constructs — the knots are symmetric, so only
/// the control points and weights reverse.
pub(super) fn reversed(curve: &NurbsCurve) -> NurbsCurve {
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points,
        weights,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn regular_hexagon() -> Vec<Point3> {
        (0..6)
            .map(|i| {
                let a = std::f64::consts::TAU * i as f64 / 6.0;
                Point3::new(a.cos(), a.sin(), 0.0)
            })
            .collect()
    }

    #[test]
    fn regular_ring_declines_balancing() {
        let hex = regular_hexagon();
        let plan = balanced_spoke_plan(
            &hex,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            3,
            1.0,
        );
        assert!(
            plan.is_none(),
            "a 120-degree split must keep straight spokes"
        );
    }

    #[test]
    fn clustered_ring_balances_to_regular_spacing() {
        // Angularly clustered even vertices: targets at 10, 30, 200
        // degrees — natural spans 20/170/170.
        let deg = |d: f64| d.to_radians();
        let mut hex = Vec::new();
        for a in [
            deg(350.0),
            deg(0.0),
            deg(10.0),
            deg(20.0),
            deg(30.0),
            deg(115.0),
        ] {
            hex.push(Point3::new(a.cos(), a.sin(), 0.0));
        }
        // Spoke targets are indices 2, 4, 0: angles 10, 30, 350.
        let plan = balanced_spoke_plan(
            &hex,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            3,
            1.0,
        )
        .expect("clustered split must balance");
        assert!(plan.natural_max_span > deg(150.0));
        // Balanced dirs must be exactly 120 degrees apart.
        for k in 0..3 {
            let a = plan.dirs[k];
            let b = plan.dirs[(k + 1) % 3];
            let cos = a.dot(&b);
            assert!(
                (cos - deg(120.0).cos()).abs() < 1e-12,
                "spoke pair {k} spacing cos {cos}"
            );
        }
    }

    #[test]
    fn curved_spoke_endpoints_and_tangent() {
        use popkern_core::geometry::CurveEval;
        let v = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(0.0, 1.0, 0.0);
        let spoke = curved_spoke(v, c, dir).expect("spoke");
        assert!(spoke.point_at(0.0).distance_to(&v) < 1e-14);
        assert!(spoke.point_at(1.0).distance_to(&c) < 1e-14);
        // Arrives at C travelling along -dir.
        let near_end = spoke.point_at(1.0 - 1e-6);
        let arrive = super::super::gregory_parallelogram::sub3(c, near_end);
        let arrive_n = arrive.normalized().expect("nonzero");
        assert!(
            arrive_n.dot(&(dir * -1.0)) > 0.999,
            "arrival direction must be -dir_at_c"
        );
        // Reversal preserves geometry.
        let rev = reversed(&spoke);
        assert!(rev.point_at(0.0).distance_to(&c) < 1e-14);
        assert!(rev.point_at(1.0).distance_to(&v) < 1e-14);
        assert!(
            rev.point_at(0.25).distance_to(&spoke.point_at(0.75)) < 1e-12,
            "reversed spoke must trace the same geometry"
        );
    }
}
