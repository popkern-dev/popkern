//! Reparameterization-aware G1/G2/G3 join analysis between two curves
//! (PK-507). A curvature comb reports jumps in |kappa|; a join test
//! must also see DIRECTION — two equal-radius arcs bending opposite
//! ways have equal |kappa| and are not G2 — and must not mistake a
//! nonlinear reparameterization for a break.
//!
//! With `x1, x2, x3` the derivatives of the first curve at its end and
//! `y1, y2, y3` of the second at its start, a Gk join is a
//! reparameterization `t(s)` with `t' = beta > 0`, `t'' = gamma`,
//! `t''' = delta` matching them: `y1 = beta x1`, `y2 = beta^2 x2 +
//! gamma x1`, `y3 = beta^3 x3 + 3 beta gamma x2 + delta x1`. Projecting
//! on `x1` (`q = x1 . x1`) gives `beta = (y1 . x1) / q`, `gamma = ((y2 -
//! beta^2 x2) . x1) / q`, `delta = ((y3 - beta^3 x3 - 3 beta gamma x2) .
//! x1) / q`; the residuals of the identities are the classifier's G2
//! and G3 evidence. Alongside it the geometric measures the identities
//! imply, each its own number: the gap, the tangent angle, the jump of
//! the curvature VECTOR `kappa N = (q a - h t) / q^2` and of its
//! arc-length derivative `[q^2 j - 3 q h a - q h' t + 4 h^2 t] / (q^3
//! sqrt q)` (`t, a, j` = first to third derivatives, `h = t . a`,
//! `h' = a . a + t . j`). Inflections (`kappa = 0`) need no special
//! case: the vectors are zero there, not undefined.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::{Point3, PortableMath, Vec3};
use serde::Serialize;

/// A curve's point and first three derivatives at one parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JoinJet {
    /// The point.
    pub point: Point3,
    /// First, second and third derivatives.
    pub d: [Vec3; 3],
}

impl JoinJet {
    /// The jet of `curve` at `t`: exact for NURBS
    /// ([`popkern_core::geometry::NurbsCurve::derivatives`]), the
    /// curve's own derivatives otherwise.
    pub fn of(curve: &CurveGeometry, t: f64) -> JoinJet {
        if let CurveGeometry::Nurbs(n) = curve {
            let d = n.derivatives(t, 3);
            return JoinJet {
                point: Point3::new(d[0].x, d[0].y, d[0].z),
                d: [d[1], d[2], d[3]],
            };
        }
        JoinJet {
            point: curve.point_at(t),
            d: [
                curve.derivative_at(t, 1),
                curve.derivative_at(t, 2),
                curve.derivative_at(t, 3),
            ],
        }
    }

    /// The curvature vector `kappa N` (zero at an inflection).
    pub fn curvature_vector(&self) -> Option<Vec3> {
        let [t, a, _] = self.d;
        let q = t.dot(&t);
        (q > 0.0).then(|| (a * q - t * t.dot(&a)) * (1.0 / (q * q)))
    }

    /// The curvature vector's arc-length derivative.
    pub fn curvature_derivative(&self) -> Option<Vec3> {
        let [t, a, j] = self.d;
        let q = t.dot(&t);
        if q <= 0.0 {
            return None;
        }
        let (h, hp) = (t.dot(&a), a.dot(&a) + t.dot(&j));
        Some(
            (j * (q * q) - a * (3.0 * q * h) - t * (q * hp) + t * (4.0 * h * h))
                * (1.0 / (q * q * q * q.sqrt())),
        )
    }
}

/// Thresholds a join is judged by.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JoinTolerances {
    /// Largest gap (length) for G0.
    pub gap: f64,
    /// Largest tangent angle (radians) for G1.
    pub angle: f64,
    /// Largest curvature-vector jump (1/length) for G2.
    pub curvature: f64,
    /// Largest curvature-derivative jump (1/length^2) for G3.
    pub curvature_derivative: f64,
}

/// How two curves meet.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct JoinReport {
    /// Distance between the end and the start.
    pub gap: f64,
    /// Angle between the tangents (radians; pi for a cusp).
    pub tangent_angle: f64,
    /// The reparameterization rates `beta, gamma, delta` (see the
    /// module docs; `beta <= 0` means the tangents oppose).
    pub beta: f64,
    /// See `beta`.
    pub gamma: f64,
    /// See `beta`.
    pub delta: f64,
    /// |y2 - beta^2 x2 - gamma x1| / |y2|-scale: the G2 identity's
    /// relative residual.
    pub g2_residual: f64,
    /// The G3 identity's relative residual.
    pub g3_residual: f64,
    /// |kappa_a N_a - kappa_b N_b|.
    pub curvature_jump: f64,
    /// |d(kappa N)/ds| jump.
    pub curvature_derivative_jump: f64,
    /// The continuity reached, 0..=3 (`-1`: not even G0).
    pub continuity: i32,
}

/// Judge the join of a curve ending in `end` to one starting in
/// `start`.
pub fn analyze_join(end: &JoinJet, start: &JoinJet, tol: &JoinTolerances) -> JoinReport {
    let [x1, x2, x3] = end.d;
    let [y1, y2, y3] = start.d;
    let gap = end.point.distance_to(&start.point);
    let q = x1.dot(&x1);
    let tangent_angle = {
        let c = x1.cross(&y1).length();
        c.patan2(x1.dot(&y1))
    };
    let (beta, gamma, delta) = if q > 0.0 {
        let beta = y1.dot(&x1) / q;
        let gamma = (y2 - x2 * (beta * beta)).dot(&x1) / q;
        let delta = (y3 - x3 * beta.powi(3) - x2 * (3.0 * beta * gamma)).dot(&x1) / q;
        (beta, gamma, delta)
    } else {
        (0.0, 0.0, 0.0)
    };
    let scale = |v: Vec3, w: Vec3| v.length().max(w.length()).max(1e-300);
    let r2 = y2 - x2 * (beta * beta) - x1 * gamma;
    let r3 = y3 - x3 * beta.powi(3) - x2 * (3.0 * beta * gamma) - x1 * delta;
    let g2_residual = r2.length() / scale(y2, x2 * (beta * beta));
    let g3_residual = r3.length() / scale(y3, x3 * beta.powi(3));
    let jump = |a: Option<Vec3>, b: Option<Vec3>| match (a, b) {
        (Some(a), Some(b)) => (a - b).length(),
        _ => f64::INFINITY,
    };
    let curvature_jump = jump(end.curvature_vector(), start.curvature_vector());
    let curvature_derivative_jump = jump(end.curvature_derivative(), start.curvature_derivative());
    let continuity = if gap > tol.gap {
        -1
    } else if beta <= 0.0 || tangent_angle > tol.angle {
        0
    } else if curvature_jump > tol.curvature {
        1
    } else if curvature_derivative_jump > tol.curvature_derivative {
        2
    } else {
        3
    };
    JoinReport {
        gap,
        tangent_angle,
        beta,
        gamma,
        delta,
        g2_residual,
        g3_residual,
        curvature_jump,
        curvature_derivative_jump,
        continuity,
    }
}

#[cfg(test)]
#[path = "curve_join_tests.rs"]
mod tests;
