//! Co-parameterized cylinder pcurves for plane-cylinder conics.
//!
//! The 3D conic is a rational quadratic whose speed varies along each
//! quarter span, so its azimuth is NOT linear (nor rational) in the
//! curve parameter. A pcurve authored as `u = 2πt` or from projected
//! control points therefore disagrees with the 3D curve at the same
//! `t` by up to 0.012 rad mid-quarter — and every consumer that
//! clips the record at one parameter window cuts the 3D curve and its
//! cylinder pcurve at DIFFERENT points (measured: a cap circle
//! clipped to a second cylinder's disk ended its pcurve 0.7° short of
//! the rim vertex the 3D end reached, and the face split's crossing
//! inventory against the rim decided the union's topology by luck).
//!
//! The axial chart coordinate `v` is an affine function of the 3D
//! point, so it stays exact through the conic's own control net. The
//! azimuth is interpolated as a C1 cubic Hermite with the curve's
//! exact azimuth rate at every sample, doubling the sample count
//! until every interval's interior check points land back on the 3D
//! curve within a fraction of the linear tolerance.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Fraction of the linear tolerance the azimuth fit may deviate from
/// the 3D curve, in model space.
const FIT_FRACTION: f64 = 1.0e-3;

/// Interior check points per interval, as fractions of the interval.
const CHECKS: [f64; 3] = [0.25, 0.5, 0.75];

/// Sample-count cap (intervals); the fit converges at fourth order,
/// so a smooth conic never approaches it.
const MAX_INTERVALS: usize = 4096;

/// The cylinder chart pcurve of `curve` (a closed or open conic lying
/// on the cylinder `origin`/`axis`/`radius`), co-parameterized with
/// it: `cylinder(pcurve(t)) == curve(t)` within `FIT_FRACTION` of the
/// linear tolerance for every `t`. `None` when a sample sits on the
/// axis (no azimuth) or the fit cannot reach the gate under the cap;
/// the caller then keeps its control-net pcurve.
pub(super) fn cylinder_conic_pcurve(
    curve: &NurbsCurve,
    origin: Point3,
    axis: Vec3,
    radius: f64,
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    let axis_len_sq = axis.dot(&axis);
    if !(axis_len_sq > 0.0 && axis_len_sq.is_finite()) {
        return None;
    }
    let (p, q) = crate::geometry::analytic_surfaces::orthonormal_basis(axis);
    let (t0, t1) = curve.parameter_range();
    if !(t1 - t0).is_finite() || t1 - t0 <= tolerance.numerical_floor() {
        return None;
    }
    let gate = FIT_FRACTION * tolerance.linear();
    let chart = Chart {
        origin,
        axis,
        axis_len_sq,
        p,
        q,
        radius,
    };
    let mut n = 16usize;
    while n <= MAX_INTERVALS {
        let samples = chart.samples(curve, t0, t1, n)?;
        if chart.worst_deviation(curve, &samples) <= gate {
            return Some(hermite_nurbs(&samples));
        }
        n *= 2;
    }
    None
}

/// One chart sample: parameter, unwrapped azimuth and its rate, axial
/// height.
#[derive(Clone, Copy)]
struct Sample {
    t: f64,
    u: f64,
    du: f64,
    v: f64,
    dv: f64,
}

struct Chart {
    origin: Point3,
    axis: Vec3,
    axis_len_sq: f64,
    p: Vec3,
    q: Vec3,
    radius: f64,
}

impl Chart {
    /// Azimuth (principal), its rate, height and its rate at `t`.
    fn sample(&self, curve: &NurbsCurve, t: f64) -> Option<(f64, f64, f64, f64)> {
        use popkern_math::PortableMath;
        let d = curve.point_at(t) - self.origin;
        let dt = curve.derivative_at(t, 1);
        // Chart conventions of `SurfaceGeometry::Cylinder::point_at`:
        // basis from the raw axis, height in units of the raw axis.
        let v = d.dot(&self.axis) / self.axis_len_sq;
        let radial = d - self.axis * v;
        let r2 = radial.dot(&radial);
        let off_axis = r2 > 1.0e-24 * self.radius.max(1.0).powi(2);
        if !off_axis {
            return None;
        }
        let u = radial.dot(&self.q).patan2(radial.dot(&self.p));
        let du = self.axis.dot(&radial.cross(&dt)) / r2;
        Some((u, du, v, dt.dot(&self.axis) / self.axis_len_sq))
    }

    fn samples(&self, curve: &NurbsCurve, t0: f64, t1: f64, n: usize) -> Option<Vec<Sample>> {
        let tau = std::f64::consts::TAU;
        let mut out: Vec<Sample> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = if i == n {
                t1
            } else {
                t0 + (t1 - t0) * i as f64 / n as f64
            };
            let (mut u, du, v, dv) = self.sample(curve, t)?;
            match out.last() {
                Some(prev) => u -= tau * ((u - prev.u) / tau).round(),
                None => u = u.rem_euclid(tau),
            }
            if !(u.is_finite() && du.is_finite() && v.is_finite() && dv.is_finite()) {
                return None;
            }
            out.push(Sample { t, u, du, v, dv });
        }
        Some(out)
    }

    fn point(&self, u: f64, v: f64) -> Point3 {
        use popkern_math::PortableMath;
        self.origin
            + self.p * (self.radius * u.pcos())
            + self.q * (self.radius * u.psin())
            + self.axis * v
    }

    fn worst_deviation(&self, curve: &NurbsCurve, samples: &[Sample]) -> f64 {
        let mut worst = 0.0f64;
        for w in samples.windows(2) {
            for f in CHECKS {
                let (u, v) = hermite(&w[0], &w[1], f);
                let t = w[0].t + (w[1].t - w[0].t) * f;
                worst = worst.max(self.point(u, v).distance_to(&curve.point_at(t)));
            }
        }
        worst
    }
}

/// Cubic Hermite interpolation of `(u, v)` at fraction `f` of the
/// interval `[a, b]`.
fn hermite(a: &Sample, b: &Sample, f: f64) -> (f64, f64) {
    let h = b.t - a.t;
    let (f2, f3) = (f * f, f * f * f);
    let h00 = 2.0 * f3 - 3.0 * f2 + 1.0;
    let h10 = f3 - 2.0 * f2 + f;
    let h01 = -2.0 * f3 + 3.0 * f2;
    let h11 = f3 - f2;
    (
        h00 * a.u + h10 * h * a.du + h01 * b.u + h11 * h * b.du,
        h00 * a.v + h10 * h * a.dv + h01 * b.v + h11 * h * b.dv,
    )
}

/// The Hermite spline as a non-rational cubic NURBS: one Bézier
/// segment per interval, interior knots of full multiplicity 3 (the
/// spline is C1 by construction; the knots only delimit segments).
fn hermite_nurbs(samples: &[Sample]) -> NurbsCurve {
    let mut knots = vec![samples[0].t; 4];
    let mut control_points = vec![Point3::new(samples[0].u, samples[0].v, 0.0)];
    for (i, w) in samples.windows(2).enumerate() {
        let (a, b) = (&w[0], &w[1]);
        let h = (b.t - a.t) / 3.0;
        control_points.push(Point3::new(a.u + h * a.du, a.v + h * a.dv, 0.0));
        control_points.push(Point3::new(b.u - h * b.du, b.v - h * b.dv, 0.0));
        control_points.push(Point3::new(b.u, b.v, 0.0));
        let last = i + 2 == samples.len();
        knots.extend(std::iter::repeat_n(b.t, if last { 4 } else { 3 }));
    }
    let weights = vec![1.0; control_points.len()];
    NurbsCurve {
        degree: 3,
        knots,
        control_points,
        weights,
    }
}

#[cfg(test)]
#[path = "plane_cylinder_pcurve_tests.rs"]
mod tests;
