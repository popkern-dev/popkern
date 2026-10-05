//! Precomputed frame, feasibility root-finding, branch sampling, and
//! serving verification for the plane×torus spiric solver — split
//! from [`super::plane_torus_spiric`] to honour the 500-line file
//! cap. See that module's docs for the mathematics.

use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;
use std::f64::consts::TAU;

use super::plane_torus_spiric::SCAN_SAMPLES;

/// Precomputed frame and coefficients for one plane×torus pair.
pub(super) struct SpiricFrame {
    plane_origin: Point3,
    n_hat: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    center: Point3,
    axis_n: Vec3,
    p1: Vec3,
    p2: Vec3,
    big_r: f64,
    r: f64,
    n1: f64,
    n2: f64,
    na: f64,
    d0: f64,
}

impl SpiricFrame {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        plane_origin: Point3,
        n_hat: Vec3,
        plane_u: Vec3,
        plane_v: Vec3,
        center: Point3,
        axis_n: Vec3,
        big_r: f64,
        r: f64,
    ) -> Option<Self> {
        let (p1, p2) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_n);
        let n1 = n_hat.dot(&p1);
        let n2 = n_hat.dot(&p2);
        let na = n_hat.dot(&axis_n);
        let d0 = (plane_origin - center).dot(&n_hat);
        if !(n1.is_finite() && n2.is_finite() && na.is_finite() && d0.is_finite()) {
            return None;
        }
        // The radial normal component must be meaningful — the
        // axis-perpendicular case is owned by the circle arm.
        if (n1 * n1 + n2 * n2).sqrt() < 1.0e-9 {
            return None;
        }
        Some(Self {
            plane_origin,
            n_hat,
            plane_u,
            plane_v,
            center,
            axis_n,
            p1,
            p2,
            big_r,
            r,
            n1,
            n2,
            na,
            d0,
        })
    }

    fn m_at(&self, u: f64) -> f64 {
        self.n1 * u.pcos() + self.n2 * u.psin()
    }

    pub(super) fn f_at(&self, u: f64) -> f64 {
        let m = self.m_at(u);
        let a = self.d0 - self.big_r * m;
        let k_sq = self.r * self.r * (m * m + self.na * self.na);
        a * a - k_sq
    }

    fn f_deriv_at(&self, u: f64) -> f64 {
        let m = self.m_at(u);
        let m_prime = -self.n1 * u.psin() + self.n2 * u.pcos();
        let a = self.d0 - self.big_r * m;
        -2.0 * m_prime * (self.big_r * a + self.r * self.r * m)
    }

    /// All roots of `F` on `[0, TAU)`, by dense sign-change bracketing
    /// plus Newton–bisection refinement.
    pub(super) fn feasibility_roots(&self) -> Vec<f64> {
        let mut roots = Vec::new();
        let dt = TAU / SCAN_SAMPLES as f64;
        let mut prev_u = 0.0;
        let mut prev_f = self.f_at(0.0);
        for i in 1..=SCAN_SAMPLES {
            let u = i as f64 * dt;
            let f = self.f_at(u);
            if prev_f == 0.0 {
                roots.push(prev_u);
            } else if f != 0.0
                && prev_f.is_sign_positive() != f.is_sign_positive()
                && let Some(root) = self.refine_root(prev_u, u)
            {
                roots.push(root);
            }
            prev_u = u;
            prev_f = f;
        }
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        roots.dedup_by(|a, b| (*a - *b).abs() < 1.0e-10);
        roots
    }

    fn refine_root(&self, mut lo: f64, mut hi: f64) -> Option<f64> {
        let mut f_lo = self.f_at(lo);
        let mut u = 0.5 * (lo + hi);
        for _ in 0..80 {
            let f = self.f_at(u);
            if f == 0.0 || (hi - lo) < 1.0e-13 {
                return Some(u);
            }
            if f.is_sign_positive() == f_lo.is_sign_positive() {
                lo = u;
                f_lo = f;
            } else {
                hi = u;
            }
            // Newton step when it stays inside the bracket, else bisect.
            let df = self.f_deriv_at(u);
            let newton = if df != 0.0 { u - f / df } else { f64::NAN };
            u = if newton.is_finite() && newton > lo && newton < hi {
                newton
            } else {
                0.5 * (lo + hi)
            };
        }
        Some(u)
    }

    /// Maximal `F < 0` arcs between consecutive roots (cyclic; the
    /// wrap arc runs past TAU).
    pub(super) fn active_arcs(&self, roots: &[f64]) -> Vec<(f64, f64)> {
        let mut arcs = Vec::new();
        for (i, &ua) in roots.iter().enumerate() {
            let ub = if i + 1 < roots.len() {
                roots[i + 1]
            } else {
                roots[0] + TAU
            };
            if ub - ua < 1.0e-9 {
                continue;
            }
            if self.f_at(0.5 * (ua + ub)) < 0.0 {
                arcs.push((ua, ub));
            }
        }
        arcs
    }

    /// Tube angle for branch `sign` at `u`, continuity-unwrapped
    /// against `prev`.
    fn v_branch(&self, u: f64, sign: f64, prev: Option<f64>) -> f64 {
        let m = self.m_at(u);
        let a = self.d0 - self.big_r * m;
        let k = self.r * (m * m + self.na * self.na).sqrt();
        let ratio = if k > 0.0 {
            (a / k).clamp(-1.0, 1.0)
        } else {
            1.0
        };
        let phi = self.na.patan2(m);
        let mut v = phi + sign * ratio.pacos();
        if let Some(p) = prev {
            v -= TAU * ((v - p) / TAU).round();
        }
        v
    }

    fn sample_at(&self, u: f64, v: f64, swapped: bool) -> MarchSample {
        let rho = self.big_r + self.r * v.pcos();
        let point = self.center
            + self.p1 * (rho * u.pcos())
            + self.p2 * (rho * u.psin())
            + self.axis_n * (self.r * v.psin());
        let d = point - self.plane_origin;
        let ul = self.plane_u.dot(&self.plane_u);
        let vl = self.plane_v.dot(&self.plane_v);
        let pu = if ul > 1e-30 {
            d.dot(&self.plane_u) / ul
        } else {
            d.dot(&self.plane_u)
        };
        let pv = if vl > 1e-30 {
            d.dot(&self.plane_v) / vl
        } else {
            d.dot(&self.plane_v)
        };
        if swapped {
            MarchSample {
                u1: u,
                v1: v,
                u2: pu,
                v2: pv,
                point,
            }
        } else {
            MarchSample {
                u1: pu,
                v1: pv,
                u2: u,
                v2: v,
                point,
            }
        }
    }

    /// One full-period loop on tube branch `sign` (the root-free case:
    /// the plane cuts every tube circle). The chain ends exactly one
    /// axis period after it starts, closing the loop.
    pub(super) fn full_period_branch(
        &self,
        sign: f64,
        n: usize,
        swapped: bool,
    ) -> Vec<MarchSample> {
        let mut samples = Vec::with_capacity(n);
        let mut prev_v: Option<f64> = None;
        for i in 0..n {
            let u = TAU * i as f64 / (n - 1) as f64;
            let v = self.v_branch(u, sign, prev_v);
            prev_v = Some(v);
            samples.push(self.sample_at(u, v, swapped));
        }
        samples
    }

    /// One closed loop over the active arc `[ua, ub]`: forward along
    /// the `v+` branch, back along `v−`, the branches meeting at the
    /// arc ends where `|A| = K`.
    ///
    /// The tube angle has a square-root turning point at each arc end
    /// (`v − φ ~ ±√(u − u_root)`), so uniform-in-`u` samples starve
    /// exactly the high-curvature turns and the fit plateaus around
    /// 1e-4 (measured on the 30°-oblique fixture). Cosine clustering
    /// (`u = ua + (ub−ua)·(1−cos πs)/2`) puts sample density ∝ 1/√ of
    /// the distance to each end and restores fit convergence.
    pub(super) fn arc_loop(&self, ua: f64, ub: f64, n: usize, swapped: bool) -> Vec<MarchSample> {
        let u_at = |i: usize| {
            let s = i as f64 / (n - 1) as f64;
            ua + (ub - ua) * 0.5 * (1.0 - (std::f64::consts::PI * s).pcos())
        };
        let mut samples = Vec::with_capacity(2 * n);
        let mut prev_v: Option<f64> = None;
        for i in 0..n {
            let u = u_at(i);
            let v = self.v_branch(u, 1.0, prev_v);
            prev_v = Some(v);
            samples.push(self.sample_at(u, v, swapped));
        }
        for i in (1..n - 1).rev() {
            let u = u_at(i);
            let v = self.v_branch(u, -1.0, prev_v);
            prev_v = Some(v);
            samples.push(self.sample_at(u, v, swapped));
        }
        // Close the loop on the starting sample (period-consistent in
        // the tube angle).
        let v_close = self.v_branch(ua, -1.0, prev_v);
        samples.push(self.sample_at(ua, v_close, swapped));
        samples
    }

    /// Verify a fitted intersection curve SERVES: the 3D curve lies on
    /// both surfaces within `gate`, and BOTH pcurves are
    /// co-parameterized with it (`S(uv(t)) ≈ curve(t)`). The pcurve
    /// checks matter as much as the 3D ones — downstream UV clipping
    /// evaluates trim endpoints through the pcurves, and an unverified
    /// pcurve fit left a boolean split endpoint 1.8e-5 off the rim
    /// circle (measured: rim-fillet strip × cuboid wall; the split
    /// stage's edge insertion rejected it as PointOffCurve). 257
    /// check parameters are coprime with every sample ladder size so
    /// fit error cannot alias into the checks.
    pub(super) fn fitted_curve_serves(
        &self,
        curve: &crate::geometry::curve_on_surface::CurveOnSurface,
        swapped: bool,
        gate: f64,
    ) -> bool {
        let checks = 257;
        for i in 0..=checks {
            let t = i as f64 / checks as f64;
            let p = curve.curve_3d.point_at(t);
            let plane_err = (p - self.plane_origin).dot(&self.n_hat).abs();
            let d = p - self.center;
            let h = d.dot(&self.axis_n);
            let rho = (d - self.axis_n * h).length();
            let torus_err =
                (((rho - self.big_r) * (rho - self.big_r) + h * h).sqrt() - self.r).abs();
            if plane_err > gate || torus_err > gate {
                return false;
            }
            let (uv_torus, uv_plane) = if swapped {
                (curve.curve_uv1.point_at(t), curve.curve_uv2.point_at(t))
            } else {
                (curve.curve_uv2.point_at(t), curve.curve_uv1.point_at(t))
            };
            let rho_t = self.big_r + self.r * uv_torus.y.pcos();
            let torus_point = self.center
                + self.p1 * (rho_t * uv_torus.x.pcos())
                + self.p2 * (rho_t * uv_torus.x.psin())
                + self.axis_n * (self.r * uv_torus.y.psin());
            let plane_point =
                self.plane_origin + self.plane_u * uv_plane.x + self.plane_v * uv_plane.y;
            if torus_point.distance_to(&p) > gate || plane_point.distance_to(&p) > gate {
                return false;
            }
        }
        true
    }
}
