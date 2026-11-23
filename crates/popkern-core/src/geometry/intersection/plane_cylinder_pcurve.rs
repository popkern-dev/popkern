//! Co-parameterized cylinder and cone pcurves for plane sections.
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
//!
//! A cone's chart is the cylinder's with the radius growing along the
//! axis (`radius + taper · v`): its plane sections take the same fit
//! (`revolved_conic_pcurve`). The shared degree-1 chart fit could not
//! reach the tolerance on a near-parabolic ellipse at scale 1000 within
//! its sample cap and returned a pcurve 3.4 off the section.

use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Fraction of the linear tolerance the azimuth fit may deviate from
/// the 3D curve, in model space.
const FIT_FRACTION: f64 = 1.0e-3;

/// Interior check points per interval, as fractions of the interval.
const CHECKS: [f64; 3] = [0.25, 0.5, 0.75];

/// Interval-count cap; the fit converges at fourth order and refines
/// only where it strays, so a smooth conic never approaches it.
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
    revolved_conic_pcurve(curve, origin, axis, (radius, 0.0), None, tolerance)
}

/// [`cylinder_conic_pcurve`] on the chart of radius `radius + taper · v`
/// about `axis` from `origin` (a cone: its apex, radius 0 and taper the
/// tangent of its half angle). With `heights`, the chart's v window, the
/// fit holds the gate where the curve passes within it (padded by a
/// tenth) and a millionth of its extent beyond: a near-parabolic
/// ellipse cut from a frustum at scale 1000 runs 3e5 out, where no face
/// reaches and the gate is past the arithmetic.
pub(super) fn revolved_conic_pcurve(
    curve: &NurbsCurve,
    origin: Point3,
    axis: Vec3,
    (radius, taper): (f64, f64),
    heights: Option<(f64, f64)>,
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
    let gate_at = |a: &Sample, b: &Sample| match heights {
        Some((lo, hi)) => {
            let (lo, hi) = (lo.min(hi), lo.max(hi));
            let pad = 0.1 * (hi - lo);
            let (below, above) = (a.v.min(b.v), a.v.max(b.v));
            if above < lo - pad || below > hi + pad {
                gate.max(1.0e-6 * (hi - lo))
            } else {
                gate
            }
        }
        None => gate,
    };
    let chart = Chart {
        origin,
        axis,
        axis_len_sq,
        p,
        q,
        radius,
        taper,
    };
    // Intervals whose check points stray are halved, the rest kept: a
    // near-parabolic section's far end turns within a sliver of its
    // parameter (a box face's ellipse on a frustum at scale 1000), and
    // uniform doubling spent its whole cap on the flat stretches.
    let mut samples = chart.samples(curve, t0, t1, 16)?;
    while samples.len() <= MAX_INTERVALS + 1 {
        let mut refined = Vec::with_capacity(2 * samples.len());
        refined.push(samples[0]);
        let mut split = false;
        for w in samples.windows(2) {
            if chart.interval_deviation(curve, &w[0], &w[1]) > gate_at(&w[0], &w[1]) {
                refined.push(chart.sample_beside(curve, 0.5 * (w[0].t + w[1].t), &w[0])?);
                split = true;
            }
            refined.push(w[1]);
        }
        if !split {
            return Some(hermite_nurbs(&samples));
        }
        samples = refined;
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
    taper: f64,
}

impl Chart {
    /// Azimuth (principal), height, and their rates at `t` along the
    /// derivative `dt`.
    fn chart_at(&self, curve: &NurbsCurve, t: f64, dt: Vec3) -> Option<(f64, f64, f64, f64)> {
        use popkern_math::PortableMath;
        let d = curve.point_at(t) - self.origin;
        // Chart conventions of `SurfaceGeometry::Cylinder::point_at`:
        // basis from the raw axis, height in units of the raw axis.
        let v = d.dot(&self.axis) / self.axis_len_sq;
        let radial = d - self.axis * v;
        let r2 = radial.dot(&radial);
        let scale = (self.radius + self.taper * v).abs().max(1.0);
        let off_axis = r2 > 1.0e-24 * scale * scale;
        if !off_axis {
            return None;
        }
        let u = radial.dot(&self.q).patan2(radial.dot(&self.p));
        let du = self.axis.dot(&radial.cross(&dt)) / r2;
        Some((u, du, v, dt.dot(&self.axis) / self.axis_len_sq))
    }

    /// The sample at `t`, its azimuth on the branch nearest `near` (or
    /// in `[0, 2π)`).
    fn sample_near(&self, curve: &NurbsCurve, t: f64, near: Option<f64>) -> Option<Sample> {
        let tau = std::f64::consts::TAU;
        let (u, du, v, dv) = self.chart_at(curve, t, curve.derivative_at(t, 1))?;
        let u = match near {
            Some(n) => u - tau * ((u - n) / tau).round(),
            None => u.rem_euclid(tau),
        };
        [u, du, v, dv]
            .iter()
            .all(|x| x.is_finite())
            .then_some(Sample { t, u, du, v, dv })
    }

    fn point(&self, u: f64, v: f64) -> Point3 {
        use popkern_math::PortableMath;
        let r = self.radius + self.taper * v;
        self.origin + self.p * (r * u.pcos()) + self.q * (r * u.psin()) + self.axis * v
    }

    /// Samples at `n` even steps.
    fn samples(&self, curve: &NurbsCurve, t0: f64, t1: f64, n: usize) -> Option<Vec<Sample>> {
        let mut out: Vec<Sample> = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let t = if i == n {
                t1
            } else {
                t0 + (t1 - t0) * i as f64 / n as f64
            };
            let s = self.sample_near(curve, t, out.last().map(|p| p.u))?;
            out.push(s);
        }
        Some(out)
    }

    /// The sample at `t`, its azimuth on the branch nearest `near`'s.
    fn sample_beside(&self, curve: &NurbsCurve, t: f64, near: &Sample) -> Option<Sample> {
        self.sample_near(curve, t, Some(near.u))
    }

    /// How far the Hermite interval from `a` to `b` strays from the
    /// curve at its check points.
    fn interval_deviation(&self, curve: &NurbsCurve, a: &Sample, b: &Sample) -> f64 {
        CHECKS
            .iter()
            .map(|&f| {
                let (u, v) = hermite(a, b, f);
                let t = a.t + (b.t - a.t) * f;
                self.point(u, v).distance_to(&curve.point_at(t))
            })
            .fold(0.0, f64::max)
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
