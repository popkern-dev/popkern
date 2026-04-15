//! §14.8 / §15 curvature comb analysis.
//!
//! A **curvature comb** is a diagnostic visualisation for a parametric
//! curve: at each of `n` sample points along the curve, draw a "hair"
//! perpendicular to the curve whose length is proportional to the
//! local curvature `κ`. The set of hair tips traces out a second
//! curve — the *comb envelope* — whose smoothness reflects the
//! original curve's Cⁿ continuity.
//!
//! Curvature combs are the primary diagnostic tool for Class-A
//! surfacing (§15): a G2-continuous blend produces a smooth comb,
//! while a G1 blend produces a comb with a **jump discontinuity** at
//! the blend boundary (the sudden change from surface curvature to
//! rolling-ball curvature `1/r`).
//!
//! This module computes combs from:
//!
//! - a single [`NurbsCurve`] sampled uniformly in its parameter
//!   range ([`curve_curvature_comb`])
//! - a cross-section of a surface along a parameter isoline — useful
//!   for diagnosing continuity across a blend surface by taking a
//!   slice perpendicular to the blend spine
//!   ([`surface_isoline_curvature_comb`])
//!
//! Callers typically feed the resulting [`CurvatureComb`] into a
//! visualisation layer. For programmatic analysis the [`CurvatureComb::jumps`]
//! helper returns the derivative of the curvature at sample-to-sample
//! spacing, which spikes at a G1 break.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

/// One sample of a curvature comb.
///
/// The "hair" at a sample point is drawn from `point` along
/// `hair_direction` with length `hair_length`. A zero-curvature point
/// has zero hair length (a straight segment).
#[derive(Debug, Clone, Copy)]
pub struct CurvatureSample {
    /// Arc-length parameter `t` at which this sample was taken.
    /// Reported in the curve's native parameter space (not arc length).
    pub t: f64,
    /// 3D point on the curve.
    pub point: Point3,
    /// Unit vector perpendicular to the curve, pointing toward the
    /// centre of osculation (the "into the curve" direction). On a
    /// straight segment this falls back to any perpendicular unit
    /// vector, but `curvature` will be zero so the hair won't draw.
    pub hair_direction: Vec3,
    /// Curvature `κ` at this sample. A straight segment has `κ = 0`;
    /// a circle of radius `r` has `κ = 1/r`.
    pub curvature: f64,
}

impl CurvatureSample {
    /// Length of the hair drawn from [`CurvatureSample::point`] along
    /// [`CurvatureSample::hair_direction`], scaled by an externally
    /// chosen factor for visualisation.
    ///
    /// Arguments:
    /// - `scale` — visual scale factor (e.g. `1.0` for raw curvature,
    ///   or a higher value for small-curvature curves to make the
    ///   comb visible).
    pub fn hair_length(&self, scale: f64) -> f64 {
        self.curvature * scale
    }
}

/// A sampled curvature comb.
///
/// The comb is just a collection of [`CurvatureSample`]s in parameter
/// order. Aggregate accessors on this struct characterise the comb's
/// shape — max/min curvature, monotonicity, and the jump-discontinuity
/// score that §15 uses to flag G1 breaks.
#[derive(Debug, Clone)]
pub struct CurvatureComb {
    /// Samples in parameter order.
    pub samples: Vec<CurvatureSample>,
}

impl CurvatureComb {
    /// Minimum curvature across all samples (inclusive of zero for
    /// straight segments).
    pub fn min_curvature(&self) -> Option<f64> {
        self.samples
            .iter()
            .map(|s| s.curvature)
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Maximum curvature across all samples.
    pub fn max_curvature(&self) -> Option<f64> {
        self.samples
            .iter()
            .map(|s| s.curvature)
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Ratio `max(κ) / min(κ)`. §15 requires this to stay below `10`
    /// for corner patches to count as Class-A. Returns `None` if the
    /// comb has no samples or the minimum curvature is zero (straight
    /// segment).
    pub fn curvature_ratio(&self) -> Option<f64> {
        let min = self.min_curvature()?;
        let max = self.max_curvature()?;
        if min <= 1e-12 {
            return None;
        }
        Some(max / min)
    }

    /// Sample-to-sample curvature deltas (first forward differences).
    /// A long flat run has near-zero deltas; a G1 jump produces one
    /// large spike.
    ///
    /// Returns a vector with `samples.len() - 1` entries; empty if
    /// there are fewer than 2 samples.
    pub fn jumps(&self) -> Vec<f64> {
        if self.samples.len() < 2 {
            return Vec::new();
        }
        self.samples
            .windows(2)
            .map(|w| (w[1].curvature - w[0].curvature).abs())
            .collect()
    }

    /// Maximum single-step curvature jump. Combined with a threshold
    /// (e.g. `tol_angular` for the matching-angle requirement, or a
    /// fraction of `max_curvature` for a relative test), this flags
    /// G1 breaks.
    pub fn max_jump(&self) -> f64 {
        self.jumps().into_iter().fold(0.0_f64, f64::max)
    }

    /// Number of sign changes in the curvature signed by the curve's
    /// principal-normal orientation. §15 allows at most 2 for
    /// inflecting blends and 0 for convex blends. This implementation
    /// counts *unsigned* sign changes of the forward-difference of
    /// curvature — good enough for detecting oscillation in the
    /// Class-A monotonicity test.
    ///
    /// Returns 0 for combs with fewer than 3 samples.
    pub fn sign_changes(&self) -> usize {
        if self.samples.len() < 3 {
            return 0;
        }
        let diffs: Vec<f64> = self
            .samples
            .windows(2)
            .map(|w| w[1].curvature - w[0].curvature)
            .collect();
        diffs.windows(2).filter(|w| w[0] * w[1] < 0.0).count()
    }
}

/// Build a curvature comb by sampling a NURBS curve uniformly in its
/// parameter range.
///
/// `n_samples` must be `≥ 2`; smaller values are silently clamped to
/// `2` (endpoints only).
///
/// **Numerical method.** Curvature is estimated from a 3-point
/// osculating circle using `point_at` only — we never call
/// `derivative_at` because `NurbsCurve`'s default implementation uses
/// finite differences with `h = 1e-7`, and the cascaded FD required
/// for a second derivative loses ~8 digits of precision (a circle of
/// radius 2 collapses to κ ≈ 5e-5 instead of 0.5). The 3-point method
/// is well-conditioned for any smooth curve and matches the analytic
/// curvature to within the sampling density.
///
/// The `hair_direction` for each sample is the unit vector in the
/// osculating plane perpendicular to the curve tangent, pointing
/// toward the centre of the circumscribed circle. On a straight
/// segment this falls back to an arbitrary perpendicular; since
/// `curvature = 0` the hair won't draw regardless.
pub fn curve_curvature_comb(curve: &NurbsCurve, n_samples: usize) -> CurvatureComb {
    let n = n_samples.max(2);
    let (t_min, t_max) = curve.parameter_range();
    let step = (t_max - t_min) / ((n - 1) as f64);

    // Neighbour offset for the 3-point stencil. Half of one sample
    // spacing is well-separated from `t` without leaving the parameter
    // domain at the endpoints.
    let h = step * 0.5;

    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let t = t_min + step * (i as f64);
        let point = curve.point_at(t);

        let t_lo = (t - h).max(t_min);
        let t_hi = (t + h).min(t_max);
        let p_lo = curve.point_at(t_lo);
        let p_hi = curve.point_at(t_hi);

        let (kappa, hair_direction) = osculating_circle_at(&p_lo, &point, &p_hi);
        samples.push(CurvatureSample {
            t,
            point,
            hair_direction,
            curvature: kappa,
        });
    }

    CurvatureComb { samples }
}

/// 3-point osculating-circle curvature + hair direction.
///
/// Given three nearby points `a`, `b`, `c` on a smooth curve, returns
/// `(κ, n)` where `κ = 1/R` is the curvature at `b` and `n` is the
/// unit vector in the plane of the triangle pointing from `b` toward
/// the centre of the circumscribed circle. On a near-degenerate (near-
/// collinear) configuration returns `(0, arbitrary perpendicular)`.
///
/// Reference: κ = 4·Area(△abc) / (|ab|·|bc|·|ca|).
fn osculating_circle_at(a: &Point3, b: &Point3, c: &Point3) -> (f64, Vec3) {
    let ab = *b - *a;
    let bc = *c - *b;
    let ca = *a - *c;
    let len_ab = ab.length();
    let len_bc = bc.length();
    let len_ca = ca.length();

    // Degenerate: any side has zero length (samples collapse at a
    // boundary).
    if len_ab < 1e-15 || len_bc < 1e-15 || len_ca < 1e-15 {
        return (0.0, Vec3::new(0.0, 0.0, 1.0));
    }

    // Area of the triangle via the cross product magnitude.
    let cross_ab_bc = ab.cross(&bc);
    let two_area = cross_ab_bc.length();
    let denom = len_ab * len_bc * len_ca;
    if denom < 1e-30 {
        return (0.0, Vec3::new(0.0, 0.0, 1.0));
    }
    let kappa = 2.0 * two_area / denom;

    // Hair direction: in the plane of the triangle, perpendicular to
    // the tangent at `b` (estimated as the midpoint of `ab` and `bc`
    // directions, i.e. `(c - a).normalized()`), pointing toward the
    // circumscribed-circle centre.
    let plane_normal = match cross_ab_bc.normalized() {
        Some(n) => n,
        None => return (0.0, Vec3::new(0.0, 0.0, 1.0)),
    };
    let tangent_mid = (*c - *a)
        .normalized()
        .unwrap_or_else(|| ab.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0)));
    let hair = plane_normal.cross(&tangent_mid);
    let hair_unit = hair
        .normalized()
        .unwrap_or_else(|| any_perpendicular_unit(&tangent_mid));

    // Orient hair toward the arc centre. The arc centre lies on the
    // concave side of the triangle — i.e. opposite to `ab × bc`
    // rotated into the plane. We test by checking the sign of
    // `hair_unit · (midpoint(a,c) - b)`: if negative, flip.
    let mid_ac = Point3::new(0.5 * (a.x + c.x), 0.5 * (a.y + c.y), 0.5 * (a.z + c.z));
    let to_mid = mid_ac - *b;
    let oriented = if hair_unit.dot(&to_mid) < 0.0 {
        hair_unit * -1.0
    } else {
        hair_unit
    };

    (kappa, oriented)
}

/// Build a curvature comb along one parameter isoline of a surface.
///
/// For a tensor-product surface `S(u, v)`, this samples `S(u_fixed, v)`
/// for `v ∈ v_range` (or, when `along = CombAlong::U`, `S(u, v_fixed)`).
/// Curvature is measured on the 3D curve traced by the isoline, which
/// is the relevant quantity for §15 "highlight lines cross blend
/// boundaries without kink": a smooth iso-line means smooth highlights
/// along that parameter direction.
///
/// The `hair_direction` for each sample is the principal normal of
/// the iso-curve in 3D space.
pub fn surface_isoline_curvature_comb(
    surface: &dyn SurfaceEval,
    along: CombAlong,
    fixed: f64,
    sweep_range: (f64, f64),
    n_samples: usize,
) -> CurvatureComb {
    let n = n_samples.max(2);
    let step = (sweep_range.1 - sweep_range.0) / ((n - 1) as f64);
    let h = step * 0.25;

    let eval = |s: f64| -> Point3 {
        match along {
            CombAlong::V => surface.point_at(fixed, s),
            CombAlong::U => surface.point_at(s, fixed),
        }
    };

    let mut samples = Vec::with_capacity(n);
    for i in 0..n {
        let s = sweep_range.0 + step * (i as f64);
        let point = eval(s);

        // Central-difference derivatives clamped to the sweep range.
        let (s_bwd, s_fwd) = ((s - h).max(sweep_range.0), (s + h).min(sweep_range.1));
        let dt = (s_fwd - s_bwd).max(1e-12);
        let p_bwd = eval(s_bwd);
        let p_fwd = eval(s_fwd);
        let tangent = Vec3::new(
            (p_fwd.x - p_bwd.x) / dt,
            (p_fwd.y - p_bwd.y) / dt,
            (p_fwd.z - p_bwd.z) / dt,
        );

        // Second derivative via three-point stencil, clamped.
        let second = Vec3::new(
            (p_fwd.x - 2.0 * point.x + p_bwd.x) / (dt * dt * 0.25),
            (p_fwd.y - 2.0 * point.y + p_bwd.y) / (dt * dt * 0.25),
            (p_fwd.z - 2.0 * point.z + p_bwd.z) / (dt * dt * 0.25),
        );

        let t_len = tangent.length();
        let kappa = if t_len > 1e-12 {
            let cross = tangent.cross(&second);
            cross.length() / (t_len * t_len * t_len)
        } else {
            0.0
        };

        let hair_direction = if t_len > 1e-12 {
            let projected = second - tangent * (second.dot(&tangent) / (t_len * t_len));
            projected
                .normalized()
                .unwrap_or_else(|| any_perpendicular_unit(&tangent))
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };

        samples.push(CurvatureSample {
            t: s,
            point,
            hair_direction,
            curvature: kappa,
        });
    }

    CurvatureComb { samples }
}

/// Which parameter direction a surface iso-line sweeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CombAlong {
    /// Sweep along `v` with `u = fixed`.
    V,
    /// Sweep along `u` with `v = fixed`.
    U,
}

fn any_perpendicular_unit(t: &Vec3) -> Vec3 {
    // Pick an arbitrary unit vector not parallel to `t` and cross.
    let helper = if t.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    t.cross(&helper)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

#[cfg(test)]
#[path = "curvature_comb_tests.rs"]
mod tests;
