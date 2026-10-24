//! Exact pcurves for circles that are chart iso-lines of a surface of
//! revolution: the latitude circles of cylinders, cones, spheres and
//! tori (`v` fixed, `u` the azimuth) and the meridian circles of
//! spheres and tori (`u` fixed, `v` the polar or tube angle).
//!
//! An angle-parameterised `Arc` advances the chart angle by exactly
//! its own parameter, so these pcurves are `Linear` in UV and exact
//! to machine precision — no sampled fit, no fit-tolerance drift.
//! That drift is not academic: a chart fit at the 1e-3 tolerance
//! seeded the mesher's point projection a hair off a sphere's seam,
//! the projection answered on the far period sheet, and the face
//! over-filled (PK-418). Every arm recognises its class by geometry
//! (centre on the axis, normal along the axis, radius matching the
//! chart) at a relative gate and returns `None` otherwise; the caller
//! verifies whatever pcurve it receives before rewriting topology.

use popkern_math::PortableMath;
use std::f64::consts::TAU;

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

/// Relative recognition gate: a scale-free `1e-9` of the arc radius.
const GATE: f64 = 1.0e-9;

/// The exact `Linear` pcurve of `curve` on `surface` when the curve is
/// an `Arc` that is a chart iso-circle of a cylinder, cone, sphere or
/// torus; `None` for every other pairing (the caller fits a chart
/// pcurve instead). `range` is the arc's parameter window.
pub fn iso_circle_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    range: (f64, f64),
) -> Option<Pcurve> {
    let CurveGeometry::Arc {
        center,
        normal,
        radius,
        ..
    } = curve
    else {
        return None;
    };
    let gate = GATE * radius.max(1.0);
    let arc = ArcClass {
        center: *center,
        normal: *normal,
        radius: *radius,
        gate,
    };
    match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: r,
        } => {
            let axial = axis.normalized()?;
            let station = (*center - *origin).dot(&axial);
            arc.latitude(curve, range, *origin, axial, station, *r)
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let axial = axis.normalized()?;
            let station = (*center - *apex).dot(&axial);
            let chart_radius = (station * half_angle.ptan()).abs();
            arc.latitude(curve, range, *apex, axial, station, chart_radius)
        }
        SurfaceGeometry::Sphere {
            center: sc,
            radius: sr,
            axis,
        } => {
            let axial = *axis;
            let d = *center - *sc;
            let polar = (d.dot(&axial) / sr).clamp(-1.0, 1.0).pacos();
            let chart_radius = sr * polar.psin();
            arc.latitude(curve, range, *sc, axial, polar, chart_radius)
                .or_else(|| arc.sphere_meridian(curve, range, surface, *sc, *sr, axial))
        }
        SurfaceGeometry::Torus {
            center: tc,
            axis,
            major_radius,
            minor_radius,
        } => {
            let axial = axis.normalized()?;
            let d = *center - *tc;
            let axial_offset = d.dot(&axial);
            // A latitude circle's radius is its distance from the
            // axis; its tube angle follows from the offset along the
            // axis and the excess over the major radius.
            let v = axial_offset.patan2(*radius - *major_radius);
            let on_tube = (*major_radius + minor_radius * v.pcos() - *radius).abs() <= gate
                && (minor_radius * v.psin() - axial_offset).abs() <= gate;
            let latitude = if on_tube {
                arc.latitude(curve, range, *tc, axial, v, *radius)
            } else {
                None
            };
            latitude.or_else(|| {
                arc.torus_meridian(
                    curve,
                    range,
                    surface,
                    *tc,
                    axial,
                    *major_radius,
                    *minor_radius,
                )
            })
        }
        _ => None,
    }
}

/// The arc under recognition.
struct ArcClass {
    center: Point3,
    normal: Vec3,
    radius: f64,
    gate: f64,
}

impl ArcClass {
    /// `u` runs, `v = station`: the arc's centre lies on the axis
    /// through `origin`, its normal is the axis, and its radius is
    /// the chart's radius at that station.
    #[allow(clippy::too_many_arguments)]
    fn latitude(
        &self,
        curve: &CurveGeometry,
        range: (f64, f64),
        origin: Point3,
        axial: Vec3,
        station: f64,
        chart_radius: f64,
    ) -> Option<Pcurve> {
        let aligned = self.normal.dot(&axial).abs() > 1.0 - GATE;
        let d = self.center - origin;
        let on_axis = (d - axial * d.dot(&axial)).length() <= self.gate;
        if !(aligned && on_axis && (chart_radius - self.radius).abs() <= self.gate) {
            return None;
        }
        let (e1, e2) = orthonormal_basis(axial);
        let start = curve.point_at(range.0) - origin;
        let u0 = start.dot(&e2).patan2(start.dot(&e1));
        let sense = if self.normal.dot(&axial) >= 0.0 {
            1.0
        } else {
            -1.0
        };
        Some(linear(
            (u0, station),
            (u0 + sense * (range.1 - range.0), station),
            range,
        ))
    }

    /// `v` runs, `u` fixed: a great circle through both poles, cut
    /// short of a pole crossing (the polar angle is affine in the arc
    /// parameter only on one side of a pole).
    fn sphere_meridian(
        &self,
        curve: &CurveGeometry,
        range: (f64, f64),
        surface: &SurfaceGeometry,
        center: Point3,
        radius: f64,
        pole: Vec3,
    ) -> Option<Pcurve> {
        use crate::geometry::analytic_surfaces::sphere_to_local;
        let concentric = (self.center - center).length() <= self.gate;
        let great = (self.radius - radius).abs() <= self.gate;
        let through_poles = self.normal.dot(&pole).abs() <= GATE;
        if !(concentric && great && through_poles) {
            return None;
        }
        // The polar angle by atan2, not acos: acos(z / r) loses half its
        // digits at a pole (1e-16 in z / r reads 1.4e-8 in angle), so a
        // meridian ENDING on a pole failed the fold check below and fell
        // to the chart fit, whose pole end took an arbitrary azimuth
        // (measured: the cube-corner union's octant arc ran its sphere
        // pcurve to u = 2π at the pole instead of 3.84, PK-930).
        let polar = |t: f64| {
            let d = sphere_to_local(&pole, curve.point_at(t) - center);
            (d.x * d.x + d.y * d.y).sqrt().patan2(d.z)
        };
        let mid = 0.5 * (range.0 + range.1);
        let (v0, v1, vm) = (polar(range.0), polar(range.1), polar(mid));
        if (vm - 0.5 * (v0 + v1)).abs() > GATE {
            return None; // crosses a pole: v folds
        }
        let m = sphere_to_local(&pole, curve.point_at(mid) - center);
        let u = azimuth(m.y.patan2(m.x));
        let sense = sense_along_v(curve, mid, surface, u, vm)?;
        Some(linear(
            (u, v0),
            (u, v0 + sense * (range.1 - range.0)),
            range,
        ))
    }

    /// `v` runs, `u` fixed: a tube circle — radius the minor radius,
    /// centre on the spine circle, plane containing the axis.
    #[allow(clippy::too_many_arguments)]
    fn torus_meridian(
        &self,
        curve: &CurveGeometry,
        range: (f64, f64),
        surface: &SurfaceGeometry,
        center: Point3,
        axial: Vec3,
        major_radius: f64,
        minor_radius: f64,
    ) -> Option<Pcurve> {
        let c = self.center - center;
        let axial_offset = c.dot(&axial);
        let radial = c - axial * axial_offset;
        let on_spine =
            axial_offset.abs() <= self.gate && (radial.length() - major_radius).abs() <= self.gate;
        let tube = (self.radius - minor_radius).abs() <= self.gate;
        let radial_dir = radial.normalized()?;
        let contains_axis =
            self.normal.dot(&axial).abs() <= GATE && self.normal.dot(&radial_dir).abs() <= GATE;
        if !(on_spine && tube && contains_axis) {
            return None;
        }
        let (e1, e2) = orthonormal_basis(axial);
        let u = azimuth(radial.dot(&e2).patan2(radial.dot(&e1)));
        let s = curve.point_at(range.0) - self.center;
        let v0 = s.dot(&axial).patan2(s.dot(&radial_dir));
        let mid = 0.5 * (range.0 + range.1);
        let sm = curve.point_at(mid) - self.center;
        let vm = sm.dot(&axial).patan2(sm.dot(&radial_dir));
        let sense = sense_along_v(curve, mid, surface, u, vm)?;
        Some(linear(
            (u, v0),
            (u, v0 + sense * (range.1 - range.0)),
            range,
        ))
    }
}

/// `+1` when the arc's tangent at `t` runs with the chart's `v`
/// direction at `(u, v)`, `-1` against it; `None` when either is
/// degenerate.
fn sense_along_v(
    curve: &CurveGeometry,
    t: f64,
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
) -> Option<f64> {
    let tangent = curve.tangent_at(t);
    let along = surface.partial_v(u, v);
    let dot = tangent.dot(&along);
    if !(dot.is_finite() && dot != 0.0) {
        return None;
    }
    Some(dot.signum())
}

/// Fold an angle onto `[0, 2π)`.
fn azimuth(angle: f64) -> f64 {
    if angle < 0.0 { angle + TAU } else { angle }
}

fn linear(start: (f64, f64), end: (f64, f64), range: (f64, f64)) -> Pcurve {
    Pcurve {
        geometry: PcurveGeometry::Linear { start, end },
        t_start: range.0,
        t_end: range.1,
    }
}

#[cfg(test)]
#[path = "revolution_iso_pcurve_tests.rs"]
mod tests;
