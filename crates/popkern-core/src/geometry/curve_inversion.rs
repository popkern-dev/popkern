//! Curve point-inversion: given a 3-D point, find the curve parameter
//! `t` such that `curve.point_at(t) ≈ point`.
//!
//! Per ADR `2026-05-02-phase-6c-a2-constructive-apex-builder.md` §14.2
//! (Newton point-inversion robustness) and §14.5 (NURBS-vs-Arc detection).
//!
//! # Tiers
//!
//! 1. **Analytic fast path** — `Line`, `Arc`, `Circle`. Closed-form
//!    inversion; always exact within numerical precision or returns
//!    [`PointInversionError::PointOffCurve`].
//! 2. **NURBS recognised as a line or circular arc** — semantic line
//!    detection ([`CurveGeometry::as_line_segment_in_range`]) or sample-based
//!    arc detection
//!    ([`CurveGeometry::as_circular_arc`]). If the NURBS sampled
//!    points lie on a single line/circle within tolerance, fall back into
//!    Tier 1's line/arc inversion.
//! 3. **Tier 2 Newton on NURBS** — stubbed for Session 4. Returns
//!    [`PointInversionError::UnsupportedCurve`]. Implement when a real
//!    NURBS-not-an-arc case arises.
//! 4. **Tier 3 unsupported curves** — `Ellipse`, `Parabola`, `Hyperbola`
//!    return [`PointInversionError::UnsupportedCurve`]. Out of Session
//!    4 scope (not used for cylinder V-end arcs).

use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

pub use crate::geometry::curve_nearest_param::nearest_param_on_curve;

#[path = "curve_inversion_circular_arc.rs"]
mod curve_inversion_circular_arc;

#[cfg(test)]
use curve_inversion_circular_arc::circle_through_three_points;

/// Errors returned by [`invert_curve_at_point`].
#[derive(Debug, Clone, PartialEq)]
pub enum PointInversionError {
    /// The point is not within `tolerance.linear()` of any point on the curve.
    /// `distance` is the residual at the best inverted parameter.
    PointOffCurve {
        /// Distance from `point` to the curve at the best-found parameter.
        distance: f64,
    },
    /// Newton iteration diverged or hit the max-iteration cap.
    /// Currently unused (Tier 2 is stubbed); reserved for future expansion.
    NewtonDiverged {
        /// Iteration count reached before bailing.
        iters: u32,
        /// Final residual before bailing.
        last_residual: f64,
    },
    /// Multiple curve parameters give points within tolerance of `point`
    /// — caller must disambiguate (e.g. by selecting the closest).
    /// Currently unused (each Tier 1 inversion has a unique result on its
    /// expected domain); reserved for future expansion.
    AmbiguousRoots {
        /// All candidate parameters within tolerance.
        roots: Vec<f64>,
    },
    /// Inverted parameter lies outside the curve-on-edge's
    /// `[t_start, t_end]` range. The geometry-level inversion succeeded,
    /// but the resulting `t` is outside the edge's contribution to the
    /// curve.
    ParameterOutOfRange {
        /// The inverted parameter (outside range).
        t: f64,
        /// The expected `(t_start, t_end)` range.
        range: (f64, f64),
    },
    /// Curve type is not supported by the current implementation.
    /// `kind` names the curve variant.
    UnsupportedCurve {
        /// The curve kind that is unsupported.
        kind: &'static str,
    },
}

/// Parameters describing a circular arc embedded in 3-D.
///
/// Returned by [`CurveGeometry::as_circular_arc`] when a curve (analytic
/// or NURBS-stored) is recognised as a circular arc.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularArcParams {
    /// Centre of the circle the arc lies on.
    pub center: Point3,
    /// Radius of the circle.
    pub radius: f64,
    /// Unit normal of the arc's plane (right-hand rule with sweep direction).
    pub normal: Vec3,
    /// Angle of the start endpoint, measured around `normal` in the plane.
    pub start_angle: f64,
    /// Signed sweep from `start_angle` to the end endpoint
    /// (positive = CCW around `normal`).
    pub sweep_angle: f64,
    /// 3-D position of the start endpoint (redundant with center+radius+angle;
    /// kept for direct equality checks).
    pub start: Point3,
    /// 3-D position of the end endpoint.
    pub end: Point3,
}

/// Find parameter `t` such that `curve.point_at(t) ≈ point`.
///
/// Returns `Ok(t)` on success with `t` in `[curve.t_start, curve.t_end]`.
/// On any failure mode, returns the typed [`PointInversionError`] explaining
/// why; the caller chooses how to react (e.g. try a different sub-edge).
///
/// # Tier dispatch
///
/// - `Line` / `Arc` / `Circle` — analytic, exact.
/// - `Nurbs` — attempts line-segment detection, then circular-arc detection via
///   [`CurveGeometry::as_circular_arc`]; if recognised, dispatches to the
///   Tier 1 analytic path. Otherwise returns `UnsupportedCurve`
///   (Tier 2 Newton-on-NURBS is stubbed for Session 4).
/// - `Ellipse` / `Parabola` / `Hyperbola` — return `UnsupportedCurve`
///   (Tier 3, deferred).
pub fn invert_curve_at_point(
    curve: &CurveOnEdge,
    point: Point3,
    tolerance: &Tolerance,
) -> Result<f64, PointInversionError> {
    let edge_range = (curve.t_start, curve.t_end);
    match &curve.geometry {
        CurveGeometry::Line { origin, direction } => {
            invert_line(*origin, *direction, point, tolerance, edge_range)
        }
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => invert_arc_general(
            *center,
            *normal,
            *radius,
            0.0,
            std::f64::consts::TAU,
            point,
            tolerance,
            edge_range,
            true,
        ),
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => invert_arc_general(
            *center,
            *normal,
            *radius,
            *start_angle,
            *end_angle,
            point,
            tolerance,
            edge_range,
            false,
        ),
        CurveGeometry::Nurbs(_) => {
            if let Some(line) = curve
                .geometry
                .as_line_segment_in_range(edge_range, tolerance.linear())
            {
                let direction = (line.end - line.start) * (1.0 / (line.t_end - line.t_start));
                return invert_line(line.start, direction, point, tolerance, edge_range);
            }
            // Circular-arc fast path — VERIFIED against the actual
            // curve: the reconstruction measures its angles in its
            // OWN circle_basis frame, which can sit a phase apart
            // from the NURBS knot frame (measured: an exact-arc rim
            // NURBS with knots in absolute angle inverted π/2 off,
            // the "interior" parameter evaluated 21 mm away, and a
            // junction pre-split silently declined). A misphased
            // answer falls through to the generic inversion.
            {
                use crate::geometry::curve::CurveEval as _;
                if let Some(arc) = curve.geometry.as_circular_arc(tolerance)
                    && let Ok(t) = invert_arc_params(&arc, point, tolerance, edge_range)
                    && curve.geometry.point_at(t).distance_to(&point) <= tolerance.linear() * 10.0
                {
                    return Ok(t);
                }
            }
            invert_nurbs_generic(curve, point, tolerance)
        }
        CurveGeometry::Ellipse { .. } => Err(PointInversionError::UnsupportedCurve {
            kind: "Ellipse (Tier 1 deferred)",
        }),
        CurveGeometry::Parabola { .. } => Err(PointInversionError::UnsupportedCurve {
            kind: "Parabola (Tier 3 deferred)",
        }),
        CurveGeometry::Hyperbola { .. } => Err(PointInversionError::UnsupportedCurve {
            kind: "Hyperbola (Tier 3 deferred)",
        }),
    }
}

/// Invert a point onto an arbitrary NURBS curve over the edge window:
/// 64-sample coarse scan for the nearest bracket, then ternary
/// refinement of the squared distance. Unimodal within one sample
/// bracket for the tolerances and curvatures this kernel handles;
/// the caller-facing contract stays "on-curve or PointOffCurve".
/// (Callers that must resolve MULTI-SHEET curves — a many-turn sweep
/// seam — use [`nearest_param_on_curve`], whose coarse density adapts
/// to the knot-span count; this site keeps its measured behavior.)
fn invert_nurbs_generic(
    curve: &CurveOnEdge,
    point: Point3,
    tolerance: &Tolerance,
) -> Result<f64, PointInversionError> {
    use crate::geometry::curve::CurveEval as _;
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    if !(hi - lo).is_finite() || hi - lo <= 0.0 {
        return Err(PointInversionError::PointOffCurve {
            distance: f64::INFINITY,
        });
    }
    let d_at = |t: f64| curve.geometry.point_at(t).distance_to(&point);
    const N: usize = 64;
    let mut best = (lo, f64::INFINITY);
    for i in 0..=N {
        let t = lo + (hi - lo) * i as f64 / N as f64;
        let d = d_at(t);
        if d < best.1 {
            best = (t, d);
        }
    }
    let step = (hi - lo) / N as f64;
    let (mut a, mut b) = ((best.0 - step).max(lo), (best.0 + step).min(hi));
    for _ in 0..64 {
        let m1 = a + (b - a) / 3.0;
        let m2 = b - (b - a) / 3.0;
        if d_at(m1) <= d_at(m2) {
            b = m2;
        } else {
            a = m1;
        }
    }
    let t = 0.5 * (a + b);
    let d = d_at(t);
    if d <= tolerance.linear() * 10.0 {
        Ok(t)
    } else {
        Err(PointInversionError::PointOffCurve { distance: d })
    }
}

// ---------------------------------------------------------------------------
// Line
// ---------------------------------------------------------------------------

fn invert_line(
    origin: Point3,
    direction: Vec3,
    point: Point3,
    tolerance: &Tolerance,
    edge_range: (f64, f64),
) -> Result<f64, PointInversionError> {
    let dir_sq = direction.dot(&direction);
    if dir_sq < 1e-30 {
        return Err(PointInversionError::UnsupportedCurve {
            kind: "Line with zero-length direction",
        });
    }
    let diff = point - origin;
    let t = diff.dot(&direction) / dir_sq;
    // Residual: distance from point to the line at parameter t.
    let on_line = origin + direction * t;
    let residual = point.distance_to(&on_line);
    if residual > tolerance.linear() {
        return Err(PointInversionError::PointOffCurve { distance: residual });
    }
    in_range_or_out(t, edge_range)
}

// ---------------------------------------------------------------------------
// Arc / Circle (general angular inversion)
// ---------------------------------------------------------------------------

/// Invert a point onto a circular arc (or full circle).
///
/// `is_full_circle` selects the angle-wrapping policy:
/// - `true` (Circle): result is wrapped into `[start_angle, start_angle + 2π)`.
/// - `false` (Arc): result must lie in `[start_angle, end_angle]` modulo 2π;
///   otherwise returns `ParameterOutOfRange`.
#[allow(clippy::too_many_arguments)]
fn invert_arc_general(
    center: Point3,
    normal: Vec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    point: Point3,
    tolerance: &Tolerance,
    edge_range: (f64, f64),
    is_full_circle: bool,
) -> Result<f64, PointInversionError> {
    if radius <= 0.0 {
        return Err(PointInversionError::UnsupportedCurve {
            kind: "Arc/Circle with non-positive radius",
        });
    }
    let n = normal
        .normalized()
        .ok_or(PointInversionError::UnsupportedCurve {
            kind: "Arc/Circle with zero normal",
        })?;
    let (u, v) = circle_basis_for(n);

    // Project point onto the arc plane: subtract its component along n.
    let delta = point - center;
    let along_normal = delta.dot(&n);
    let in_plane = delta - n * along_normal;

    // Compute angle in (u, v) basis.
    let cos_t = in_plane.dot(&u);
    let sin_t = in_plane.dot(&v);
    let raw_angle = sin_t.atan2(cos_t);
    let in_plane_radius = (cos_t * cos_t + sin_t * sin_t).sqrt();

    // Verify on circle: in-plane radius matches AND off-plane component is 0.
    let radial_residual = (in_plane_radius - radius).abs();
    let total_residual = (radial_residual * radial_residual + along_normal * along_normal).sqrt();
    if total_residual > tolerance.linear() {
        return Err(PointInversionError::PointOffCurve {
            distance: total_residual,
        });
    }

    // Wrap raw_angle into the arc's angular domain.
    let t = wrap_angle_into(raw_angle, start_angle, end_angle, is_full_circle)?;

    in_range_or_out(t, edge_range)
}

fn invert_arc_params(
    arc: &CircularArcParams,
    point: Point3,
    tolerance: &Tolerance,
    edge_range: (f64, f64),
) -> Result<f64, PointInversionError> {
    let end_angle = arc.start_angle + arc.sweep_angle;
    invert_arc_general(
        arc.center,
        arc.normal,
        arc.radius,
        arc.start_angle.min(end_angle),
        arc.start_angle.max(end_angle),
        point,
        tolerance,
        edge_range,
        false,
    )
}

/// Wrap `raw_angle` (in `[-π, π]`) into the arc's angular domain.
///
/// For full circles, returns `raw_angle` shifted into
/// `[start_angle, start_angle + 2π)`.
///
/// For arcs, accepts any angle that lies within
/// `[start_angle, end_angle]` modulo 2π (with a small tolerance for
/// endpoints); otherwise returns `ParameterOutOfRange`.
fn wrap_angle_into(
    raw_angle: f64,
    start_angle: f64,
    end_angle: f64,
    is_full_circle: bool,
) -> Result<f64, PointInversionError> {
    use std::f64::consts::TAU;
    let (lo, hi) = if start_angle <= end_angle {
        (start_angle, end_angle)
    } else {
        (end_angle, start_angle)
    };

    if is_full_circle {
        // Map raw_angle (in [-π, π]) into [start_angle, start_angle + 2π).
        let mut t = raw_angle;
        while t < start_angle {
            t += TAU;
        }
        while t >= start_angle + TAU {
            t -= TAU;
        }
        return Ok(t);
    }

    // Arc: try raw_angle, raw_angle + 2π, raw_angle - 2π; pick the one
    // closest to the arc's domain.
    let candidates = [raw_angle, raw_angle + TAU, raw_angle - TAU];
    let endpoint_eps = 1e-12;
    let mut in_range: Option<f64> = None;
    let mut best_outside: f64 = f64::INFINITY;
    let mut best_outside_t: f64 = raw_angle;
    for &t in &candidates {
        if t >= lo - endpoint_eps && t <= hi + endpoint_eps {
            in_range = Some(t.clamp(lo, hi));
            break;
        }
        let dist_outside = (t - lo).abs().min((t - hi).abs());
        if dist_outside < best_outside {
            best_outside = dist_outside;
            best_outside_t = t;
        }
    }
    in_range.ok_or(PointInversionError::ParameterOutOfRange {
        t: best_outside_t,
        range: (lo, hi),
    })
}

/// Basis (u, v) for a circle/arc lying in a plane with the given unit normal.
///
/// `u × v` points along `normal` (right-hand rule). Both are unit vectors.
fn circle_basis_for(normal: Vec3) -> (Vec3, Vec3) {
    // Pick an axis not nearly parallel to normal.
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = normal
        .cross(&axis)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

fn in_range_or_out(t: f64, range: (f64, f64)) -> Result<f64, PointInversionError> {
    let (lo, hi) = if range.0 <= range.1 {
        range
    } else {
        (range.1, range.0)
    };
    let edge_eps = 1e-9 * (hi - lo).abs().max(1.0);
    if t < lo - edge_eps || t > hi + edge_eps {
        Err(PointInversionError::ParameterOutOfRange { t, range })
    } else {
        Ok(t.clamp(lo, hi))
    }
}
#[cfg(test)]
#[path = "curve_inversion_tests.rs"]
mod tests;
