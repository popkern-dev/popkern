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

use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

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
                invert_line(line.start, direction, point, tolerance, edge_range)
            } else if let Some(arc) = curve.geometry.as_circular_arc(tolerance) {
                invert_arc_params(&arc, point, tolerance, edge_range)
            } else {
                Err(PointInversionError::UnsupportedCurve {
                    kind: "Nurbs (Tier 2 Newton not yet implemented)",
                })
            }
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

// ---------------------------------------------------------------------------
// CurveGeometry::as_circular_arc — sample-based detector
// ---------------------------------------------------------------------------

impl CurveGeometry {
    /// Detect whether this curve is geometrically a circular arc.
    ///
    /// Returns `Some(params)` if 5 evenly-spaced sample points all lie on
    /// the circle fitted through 3 of them (within `tolerance.linear()`);
    /// otherwise `None`.
    ///
    /// Works for `Arc` / `Circle` (returns the analytic representation
    /// directly) and for `Nurbs` curves that geometrically trace a
    /// circular arc (returns the fitted representation).  Returns `None`
    /// for `Line` / `Ellipse` / `Parabola` / `Hyperbola`.
    pub fn as_circular_arc(&self, tolerance: &Tolerance) -> Option<CircularArcParams> {
        match self {
            CurveGeometry::Line { .. } => None,
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => {
                let n = normal.normalized()?;
                let (u, _v) = circle_basis_for(n);
                let start = *center + u * *radius;
                Some(CircularArcParams {
                    center: *center,
                    radius: *radius,
                    normal: n,
                    start_angle: 0.0,
                    sweep_angle: std::f64::consts::TAU,
                    start,
                    end: start,
                })
            }
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            } => {
                let n = normal.normalized()?;
                let (u, v) = circle_basis_for(n);
                let start_dir = u * start_angle.cos() + v * start_angle.sin();
                let end_dir = u * end_angle.cos() + v * end_angle.sin();
                Some(CircularArcParams {
                    center: *center,
                    radius: *radius,
                    normal: n,
                    start_angle: *start_angle,
                    sweep_angle: end_angle - start_angle,
                    start: *center + start_dir * *radius,
                    end: *center + end_dir * *radius,
                })
            }
            CurveGeometry::Ellipse { .. }
            | CurveGeometry::Parabola { .. }
            | CurveGeometry::Hyperbola { .. } => None,
            CurveGeometry::Nurbs(_) => detect_circular_arc_from_samples(self, tolerance),
        }
    }
}

/// Sample-based circular-arc detector.  See `as_circular_arc`.
fn detect_circular_arc_from_samples(
    curve: &CurveGeometry,
    tolerance: &Tolerance,
) -> Option<CircularArcParams> {
    let (t_lo, t_hi) = curve.parameter_range();
    if !t_lo.is_finite() || !t_hi.is_finite() || t_hi <= t_lo {
        return None;
    }
    let frac = |x: f64| t_lo + (t_hi - t_lo) * x;
    let samples = [
        curve.point_at(frac(0.0)),
        curve.point_at(frac(0.25)),
        curve.point_at(frac(0.5)),
        curve.point_at(frac(0.75)),
        curve.point_at(frac(1.0)),
    ];

    let (center, radius, normal) = circle_through_three_points(samples[0], samples[2], samples[4])?;

    // Verify all 5 samples lie on the fitted circle within tolerance.
    let lin_tol = tolerance.linear();
    for (i, p) in samples.iter().enumerate() {
        let radial = (*p - center).length() - radius;
        let along_normal = (*p - center).dot(&normal).abs();
        let total = (radial * radial + along_normal * along_normal).sqrt();
        if total > lin_tol {
            // Sample i is off the candidate circle. Not an arc within tolerance.
            let _ = i;
            return None;
        }
    }

    // Compute start/end angles in the (u, v) basis.
    let (u, v) = circle_basis_for(normal);
    let start_vec = samples[0] - center;
    let end_vec = samples[4] - center;
    let start_angle = start_vec.dot(&v).atan2(start_vec.dot(&u));
    let end_angle = end_vec.dot(&v).atan2(end_vec.dot(&u));

    // Determine sweep direction by checking which sense puts the midpoint
    // of the parameter range close to the midpoint of the arc.
    let mid_vec = samples[2] - center;
    let mid_angle = mid_vec.dot(&v).atan2(mid_vec.dot(&u));

    let sweep_ccw = signed_arc_sweep(start_angle, end_angle, mid_angle);

    Some(CircularArcParams {
        center,
        radius,
        normal,
        start_angle,
        sweep_angle: sweep_ccw,
        start: samples[0],
        end: samples[4],
    })
}

/// Compute the signed sweep from `start` to `end` (around the chosen
/// `normal`, CCW positive) such that `mid_angle` lies between them.
///
/// Returns the sweep in radians; negate the caller's normal to flip
/// direction if needed.
fn signed_arc_sweep(start: f64, end: f64, mid: f64) -> f64 {
    use std::f64::consts::TAU;
    // Shift everything so start is 0; CCW sweep is in [0, 2π).
    let end_shifted = (end - start).rem_euclid(TAU);
    let mid_shifted = (mid - start).rem_euclid(TAU);
    if mid_shifted <= end_shifted {
        // Mid lies in CCW arc from start to end; sweep is +end_shifted.
        end_shifted
    } else {
        // Mid lies in the other arc; sweep is the CW direction.
        end_shifted - TAU
    }
}

/// Fit a circle through three points in 3-D.
///
/// Returns `(center, radius, unit_normal)` where `unit_normal` follows
/// the right-hand rule for traversal `p0 → p1 → p2`.  Returns `None`
/// if the points are collinear.
fn circle_through_three_points(p0: Point3, p1: Point3, p2: Point3) -> Option<(Point3, f64, Vec3)> {
    let v1 = p1 - p0;
    let v2 = p2 - p0;
    let n = v1.cross(&v2);
    let n_sq = n.dot(&n);
    if n_sq < 1e-30 {
        return None; // collinear or coincident
    }
    let v1_sq = v1.dot(&v1);
    let v2_sq = v2.dot(&v2);
    // Circumcenter formula in 3-D:
    //   centre = p0 + ( (|v2|² · (n × v1) + |v1|² · (v2 × n)) / (2 · |n|²) )
    let to_center = (n.cross(&v1) * v2_sq + v2.cross(&n) * v1_sq) * (1.0 / (2.0 * n_sq));
    let center = p0 + to_center;
    let radius = (p0 - center).length();
    let normal = n.normalized()?;
    Some((center, radius, normal))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::nurbs_curve::NurbsCurve;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    fn line_curve(origin: Point3, direction: Vec3, t_start: f64, t_end: f64) -> CurveOnEdge {
        CurveOnEdge {
            geometry: CurveGeometry::Line { origin, direction },
            t_start,
            t_end,
        }
    }

    fn arc_curve(
        center: Point3,
        normal: Vec3,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> CurveOnEdge {
        CurveOnEdge {
            geometry: CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle,
                end_angle,
            },
            t_start: start_angle,
            t_end: end_angle,
        }
    }

    fn circle_curve(center: Point3, normal: Vec3, radius: f64) -> CurveOnEdge {
        CurveOnEdge {
            geometry: CurveGeometry::Circle {
                center,
                normal,
                radius,
            },
            t_start: 0.0,
            t_end: TAU,
        }
    }

    // ---------- Line tier 1 ----------

    #[test]
    fn line_inversion_at_origin() {
        let line = line_curve(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            -10.0,
            10.0,
        );
        let t = invert_curve_at_point(&line, Point3::new(2.5, 0.0, 0.0), &tol()).unwrap();
        assert!((t - 2.5).abs() < 1e-12);
    }

    #[test]
    fn line_inversion_off_curve() {
        let line = line_curve(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            -10.0,
            10.0,
        );
        let result = invert_curve_at_point(&line, Point3::new(2.5, 0.5, 0.0), &tol());
        assert!(matches!(
            result,
            Err(PointInversionError::PointOffCurve { .. })
        ));
    }

    #[test]
    fn line_inversion_out_of_range() {
        let line = line_curve(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            1.0,
        );
        let result = invert_curve_at_point(&line, Point3::new(5.0, 0.0, 0.0), &tol());
        assert!(matches!(
            result,
            Err(PointInversionError::ParameterOutOfRange { .. })
        ));
    }

    // ---------- Arc tier 1 ----------

    /// Reconstruct the expected point at a given parameter via `point_at`,
    /// then invert and check round-trip.
    fn check_round_trip(curve: &CurveOnEdge, t_expected: f64) {
        let p = curve.geometry.point_at(t_expected);
        let t = invert_curve_at_point(curve, p, &tol())
            .unwrap_or_else(|e| panic!("invert failed at t={t_expected}: {e:?}"));
        assert!(
            (t - t_expected).abs() < 1e-9,
            "round-trip t={t}, expected={t_expected}, point={p:?}"
        );
    }

    #[test]
    fn arc_inversion_at_endpoints() {
        let arc = arc_curve(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            0.0,
            FRAC_PI_2,
        );
        check_round_trip(&arc, 0.0);
        check_round_trip(&arc, FRAC_PI_2);
    }

    #[test]
    fn arc_inversion_at_midpoint() {
        let arc = arc_curve(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            0.0,
            FRAC_PI_2,
        );
        check_round_trip(&arc, FRAC_PI_4);
    }

    #[test]
    fn arc_inversion_off_curve() {
        let arc = arc_curve(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            0.0,
            FRAC_PI_2,
        );
        // Off-plane: pick a point on the arc, then nudge it along normal.
        let p_on_arc = arc.geometry.point_at(FRAC_PI_4);
        let p_off_plane = Point3::new(p_on_arc.x, p_on_arc.y, 0.5);
        let result = invert_curve_at_point(&arc, p_off_plane, &tol());
        assert!(
            matches!(result, Err(PointInversionError::PointOffCurve { .. })),
            "off-plane should be off-curve, got {result:?}"
        );
        // In-plane but wrong radius (point at half the radius from centre).
        let half_radius = Point3::new(p_on_arc.x * 0.5, p_on_arc.y * 0.5, 0.0);
        let result = invert_curve_at_point(&arc, half_radius, &tol());
        assert!(
            matches!(result, Err(PointInversionError::PointOffCurve { .. })),
            "wrong radius should be off-curve, got {result:?}"
        );
    }

    #[test]
    fn arc_inversion_out_of_arc_range() {
        let arc = arc_curve(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            0.0,
            FRAC_PI_2,
        );
        // Pick a point on the same circle but at parameter t = π
        // (outside the [0, π/2] arc range).
        let circle = CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let p_on_circle_outside_arc = circle.point_at(PI);
        let result = invert_curve_at_point(&arc, p_on_circle_outside_arc, &tol());
        assert!(
            matches!(result, Err(PointInversionError::ParameterOutOfRange { .. })),
            "expected ParameterOutOfRange, got {result:?}"
        );
    }

    // ---------- Circle tier 1 ----------

    #[test]
    fn circle_inversion_round_trip_at_quarter() {
        let circle = circle_curve(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0), 1.0);
        // Round-trip at multiple parameters.
        for t_expected in [0.0_f64, FRAC_PI_2, PI, FRAC_PI_2 * 3.0] {
            let p = circle.geometry.point_at(t_expected);
            let t = invert_curve_at_point(&circle, p, &tol()).unwrap();
            // For a full circle the result is wrapped into [0, 2π).
            // 0 and 2π are the same point; allow either.
            let diff = (t - t_expected).abs().min((t - t_expected - TAU).abs());
            assert!(
                diff < 1e-9,
                "circle round-trip: got {t}, expected {t_expected}"
            );
        }
    }

    // ---------- as_circular_arc detector ----------

    #[test]
    fn as_circular_arc_recognises_arc_directly() {
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let params = arc.as_circular_arc(&tol()).unwrap();
        assert!((params.radius - 2.0).abs() < 1e-12);
        assert!((params.start_angle - 0.0).abs() < 1e-12);
        assert!((params.sweep_angle - FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn as_circular_arc_returns_none_for_line() {
        let line = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        assert!(line.as_circular_arc(&tol()).is_none());
    }

    #[test]
    fn as_circular_arc_recognises_nurbs_quarter_arc() {
        // A rational quadratic NURBS quarter circle in xy-plane, radius 1.
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let nurbs = CurveGeometry::Nurbs(arc.to_nurbs());
        let params = nurbs
            .as_circular_arc(&tol())
            .expect("NURBS quarter arc should be detected as circular");
        assert!(
            (params.radius - 1.0).abs() < 1e-6,
            "fitted radius should be ≈1, got {}",
            params.radius
        );
        // The fitted normal should be parallel to (0, 0, 1) (sign may flip
        // depending on sweep direction).
        let z_align = params.normal.dot(&Vec3::new(0.0, 0.0, 1.0)).abs();
        assert!(
            z_align > 0.999,
            "normal should be ±z, got {:?}",
            params.normal
        );
    }

    #[test]
    fn as_circular_arc_returns_none_for_non_arc_nurbs() {
        // Cubic NURBS through scattered (non-circular) points.
        let nurbs = NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.5, 0.0),
                Point3::new(1.0, 0.5, 0.0),
                Point3::new(1.5, 1.0, 0.0),
                Point3::new(2.0, 0.3, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0, 1.0],
        };
        let geom = CurveGeometry::Nurbs(nurbs);
        assert!(geom.as_circular_arc(&tol()).is_none());
    }

    #[test]
    fn invert_nurbs_arc_falls_back_to_arc_path() {
        let arc = CurveGeometry::Arc {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            start_angle: 0.0,
            end_angle: FRAC_PI_2,
        };
        let nurbs_curve = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(arc.to_nurbs()),
            t_start: 0.0,
            t_end: FRAC_PI_2,
        };
        // Pick the point at the analytic arc's midpoint and confirm the
        // NURBS-detected-as-arc tier 1b path produces the same parameter
        // (up to small fitting / detector tolerance).
        let mid_point = arc.point_at(FRAC_PI_4);
        let result = invert_curve_at_point(&nurbs_curve, mid_point, &tol());
        let t = result.expect("NURBS-stored arc should invert via Tier 1b");
        assert!(
            (t - FRAC_PI_4).abs() < 1e-3,
            "NURBS-detected arc midpoint inversion: got t={t}, expected ≈π/4"
        );
    }

    #[test]
    fn invert_nurbs_line_uses_analytic_line_path() {
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
                weights: vec![1.0, 1.0],
            }),
            t_start: 0.0,
            t_end: 1.0,
        };
        let t = invert_curve_at_point(&curve, Point3::new(2.0, 0.0, 0.0), &tol()).unwrap();
        assert!((t - 0.5).abs() < 1e-12);
    }

    #[test]
    fn invert_non_arc_nurbs_returns_unsupported() {
        let nurbs = NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.5, 0.0),
                Point3::new(1.0, 0.5, 0.0),
                Point3::new(1.5, 1.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0, 1.0],
        };
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(nurbs),
            t_start: 0.0,
            t_end: 1.0,
        };
        let result = invert_curve_at_point(&curve, Point3::new(0.5, 0.5, 0.0), &tol());
        assert!(matches!(
            result,
            Err(PointInversionError::UnsupportedCurve { .. })
        ));
    }

    #[test]
    fn ellipse_returns_unsupported() {
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Ellipse {
                center: Point3::new(0.0, 0.0, 0.0),
                u_basis: Vec3::new(2.0, 0.0, 0.0),
                v_basis: Vec3::new(0.0, 1.0, 0.0),
                start_angle: 0.0,
                end_angle: FRAC_PI_2,
            },
            t_start: 0.0,
            t_end: FRAC_PI_2,
        };
        let result = invert_curve_at_point(&curve, Point3::new(2.0, 0.0, 0.0), &tol());
        assert!(matches!(
            result,
            Err(PointInversionError::UnsupportedCurve { .. })
        ));
    }

    // ---------- circle_through_three_points ----------

    #[test]
    fn fit_circle_through_three_points_xy_plane() {
        let p0 = Point3::new(1.0, 0.0, 0.0);
        let p1 = Point3::new(0.0, 1.0, 0.0);
        let p2 = Point3::new(-1.0, 0.0, 0.0);
        let (center, radius, normal) = circle_through_three_points(p0, p1, p2).unwrap();
        assert!(center.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-12);
        assert!((radius - 1.0).abs() < 1e-12);
        assert!((normal.dot(&Vec3::new(0.0, 0.0, 1.0)) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn fit_circle_returns_none_for_collinear() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(2.0, 0.0, 0.0);
        assert!(circle_through_three_points(p0, p1, p2).is_none());
    }
}
