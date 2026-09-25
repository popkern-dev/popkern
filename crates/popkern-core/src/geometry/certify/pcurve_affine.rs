//! Degree-1 NURBS pcurves that are one affine map up to roundoff
//! (PK-542).
//!
//! The kernel authors circle rails on cylinders as fitted polylines:
//! the fit samples the exact angle-linear chart line, so the control
//! points are collinear and knot-proportional and the pcurve IS a
//! linear map (an explicit `Linear` arm was removed under PK-292). As
//! a NURBS it has one knot span per sample (64 on the corpus ring
//! rails), and every certifier interval straddles a knot: the branch
//! union then runs every seam quantity per span.
//!
//! Enclosing all spans at once is sound. Let `L` be the affine map
//! through the end control points and, per span `k`, `L_k = L + e_k`
//! its polynomial. With `δ` the certified worst node deviation
//! `|P_k − L(t_k)|` and `δ′` the worst span-slope deviation
//! `|slope_k − L′|`, every `e_k` extended over the domain (and a
//! domain width beyond it) satisfies `|e_k| ≤ δ + 2 δ′ (t₁ − t₀)` and
//! `|e_k′| ≤ δ′`, and all higher derivatives vanish. The single jet
//! `L ± that`, `L′ ± δ′`, `0` therefore contains every span
//! polynomial's jet everywhere a certifier interval reaches, so each
//! span's Taylor forms — valid, as the span is a polynomial — are
//! contained in the forms computed from it.

use super::interval::Iv;
use super::interval_vec::IvVec2;
use super::pcurve_jet::PCURVE_JET_TERMS;
use crate::geometry::nurbs_curve::NurbsCurve;

/// Largest node deviation, relative to the control-point magnitude,
/// for which a degree-1 pcurve is merged into one affine branch: a
/// roundoff-level departure from collinearity, not a genuine polyline.
const AFFINE_MERGE_RELATIVE: f64 = 1e-12;

/// A degree-1 pcurve enclosed as one affine map.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct AffinePolyline {
    t: (f64, f64),
    start: (f64, f64),
    end: (f64, f64),
    /// Per-axis value slack.
    value_dev: (f64, f64),
    /// Per-axis slope slack.
    slope_dev: (f64, f64),
}

fn pt(x: f64) -> Iv {
    Iv::point(x)
}

impl AffinePolyline {
    /// The merged enclosure of a non-rational degree-1 NURBS pcurve,
    /// or `None` when it is a genuine polyline (or not degree 1).
    pub(super) fn of(curve: &NurbsCurve) -> Option<Self> {
        let n = curve.control_points.len();
        if curve.degree != 1 || n < 3 || curve.knots.len() != n + 2 {
            return None;
        }
        if curve.weights.iter().any(|&w| w != curve.weights[0]) {
            return None;
        }
        let nodes: Vec<f64> = (0..n).map(|k| curve.knots[k + 1]).collect();
        let (t0, t1) = (nodes[0], nodes[n - 1]);
        if t1 <= t0 {
            return None;
        }
        let cp = |k: usize| (curve.control_points[k].x, curve.control_points[k].y);
        let (p0, p1) = (cp(0), cp(n - 1));
        let span = pt(t1) - pt(t0);
        let slope = ((pt(p1.0) - pt(p0.0)) / span, (pt(p1.1) - pt(p0.1)) / span);
        let line = |t: f64, axis: usize| {
            let (a, d) = if axis == 0 {
                (p0.0, slope.0)
            } else {
                (p0.1, slope.1)
            };
            pt(a) + d * (pt(t) - pt(t0))
        };
        let (mut dev, mut slope_dev) = ([0.0_f64; 2], [0.0_f64; 2]);
        for k in 0..n {
            let (x, y) = cp(k);
            dev[0] = dev[0].max((pt(x) - line(nodes[k], 0)).mag());
            dev[1] = dev[1].max((pt(y) - line(nodes[k], 1)).mag());
            if k + 1 < n {
                let width = pt(nodes[k + 1]) - pt(nodes[k]);
                let (x1, y1) = cp(k + 1);
                slope_dev[0] = slope_dev[0].max(((pt(x1) - pt(x)) / width - slope.0).mag());
                slope_dev[1] = slope_dev[1].max(((pt(y1) - pt(y)) / width - slope.1).mag());
            }
        }
        let magnitude = curve
            .control_points
            .iter()
            .map(|p| p.x.abs().max(p.y.abs()))
            .fold(1.0_f64, f64::max);
        if dev[0].max(dev[1]) > AFFINE_MERGE_RELATIVE * magnitude || !slope_dev[0].is_finite() {
            return None;
        }
        let reach = (pt(2.0) * span).mag();
        let slack = |axis: usize| (pt(dev[axis]) + pt(slope_dev[axis]) * pt(reach)).mag();
        Some(Self {
            t: (t0, t1),
            start: p0,
            end: p1,
            value_dev: (slack(0), slack(1)),
            slope_dev: (slope_dev[0], slope_dev[1]),
        })
    }

    /// Jet over the parameter interval `t`.
    pub(super) fn jet(&self, t: Iv) -> [IvVec2; PCURVE_JET_TERMS] {
        let widen = |x: Iv, d: f64| x + Iv::new(-d, d);
        let span = pt(self.t.1) - pt(self.t.0);
        let slope = (
            (pt(self.end.0) - pt(self.start.0)) / span,
            (pt(self.end.1) - pt(self.start.1)) / span,
        );
        let offset = t - pt(self.t.0);
        let mut out = [IvVec2::new(pt(0.0), pt(0.0)); PCURVE_JET_TERMS];
        out[0] = IvVec2::new(
            widen(pt(self.start.0) + slope.0 * offset, self.value_dev.0),
            widen(pt(self.start.1) + slope.1 * offset, self.value_dev.1),
        );
        out[1] = IvVec2::new(
            widen(slope.0, self.slope_dev.0),
            widen(slope.1, self.slope_dev.1),
        );
        out
    }
}
