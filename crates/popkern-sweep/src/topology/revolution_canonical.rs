//! Canonical recognition of revolved profile spans (PK-409): a
//! straight span revolves to a plane (⊥ axis), a cylinder (∥ axis) or
//! a cone; an arc centred on the axis to a sphere (world-`z` axis
//! only — the sphere chart is fixed to world axes); an off-axis arc to
//! a torus. Faces that carry an exact surface are served by the exact
//! boundary-first revolution meshers (seams and poles included) and
//! export exactly; anything unrecognised stays NURBS.
//!
//! Every recognition is VERIFIED by sampling the NURBS surface of
//! revolution against the candidate (projection round trip), never
//! trusted from construction.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use std::f64::consts::TAU;

use crate::surface::canonical::recognize_circle;
use crate::surface::canonical_curves::recognize_line;

/// The native-chart image of the `[0,1]²` revolution chart
/// (u = profile, v = angle) for a recognised surface.
#[derive(Debug, Clone)]
pub(crate) enum RevolutionChart {
    /// Angle-linear surfaces (`Cylinder`, `Cone`, `Torus`, `Sphere`):
    /// `u_native = phi0 + sense · v01 · sweep`, `v_native = v_a +
    /// u01 · (v_b − v_a)` — both linear, so an angle-parameterised
    /// edge and its pcurve stay co-parameterised.
    Angular {
        /// Azimuth of the seam (profile at angle 0) in the surface's frame.
        phi0: f64,
        /// `+1` when increasing v01 increases the native azimuth.
        sense: f64,
        /// Total angle of the sweep.
        sweep: f64,
        /// Native second coordinate at the span start.
        v_a: f64,
        /// Native second coordinate at the span end.
        v_b: f64,
    },
    /// A plane ⊥ the axis: native chart is cartesian in the plane;
    /// the seam is a radial line, circles are rational circles.
    Planar {
        /// Plane-chart coordinates of the axis point.
        center: (f64, f64),
        /// Azimuth of the seam in the plane chart.
        phi0: f64,
        /// `+1` when increasing v01 increases the plane azimuth.
        sense: f64,
        /// Radii at the span start and end.
        r_a: f64,
        r_b: f64,
    },
}

/// A recognised surface of revolution for one span.
#[derive(Debug, Clone)]
pub(crate) struct RevolvedSurface {
    /// The exact surface.
    pub geometry: SurfaceGeometry,
    /// Chart map from the `[0,1]²` revolution chart.
    pub chart: RevolutionChart,
}

fn axial_radial(p: &Point3, origin: &Point3, axis: &Vec3) -> (f64, f64) {
    let d = *p - *origin;
    let z = d.dot(axis);
    let r = (d - *axis * z).length();
    (z, r)
}

/// Azimuth of `p` about the axis in `basis` (the surface evaluator's
/// own frame).
fn azimuth(p: &Point3, origin: &Point3, axis: &Vec3, basis: (Vec3, Vec3)) -> f64 {
    let d = *p - *origin;
    let radial = d - *axis * d.dot(axis);
    radial.dot(&basis.1).atan2(radial.dot(&basis.0))
}

/// Shortest signed angular step from `a` to `b`.
fn angle_step(a: f64, b: f64) -> f64 {
    let mut d = b - a;
    while d > std::f64::consts::PI {
        d -= TAU;
    }
    while d < -std::f64::consts::PI {
        d += TAU;
    }
    d
}

/// Verify that the NURBS revolution surface lies on `candidate`
/// (projection round trip within `gate`).
fn surface_agrees(nurbs: &NurbsSurface, candidate: &SurfaceGeometry, gate: f64) -> bool {
    for i in 0..=4 {
        for j in 0..=8 {
            let p = nurbs.point_at(i as f64 / 4.0, j as f64 / 8.0);
            let Some((u, v)) = project_point_to_uv(candidate, &p) else {
                return false;
            };
            if candidate.point_at(u, v).distance_to(&p) > gate {
                return false;
            }
        }
    }
    true
}

/// The sense of the sweep in the candidate's azimuth: `+1` when the
/// NURBS v direction increases the native azimuth.
fn azimuth_sense(nurbs: &NurbsSurface, origin: &Point3, axis: &Vec3, basis: (Vec3, Vec3)) -> f64 {
    let p0 = nurbs.point_at(0.5, 0.0);
    let p1 = nurbs.point_at(0.5, 0.05);
    let step = angle_step(
        azimuth(&p0, origin, axis, basis),
        azimuth(&p1, origin, axis, basis),
    );
    if step >= 0.0 { 1.0 } else { -1.0 }
}

/// Recognise the surface a span sweeps about the axis over `sweep`
/// radians. `nurbs` is the revolution surface built from the span
/// (u = profile, v = angle in [0, 1]).
pub(crate) fn recognize_revolved_span(
    span: &NurbsCurve,
    nurbs: &NurbsSurface,
    axis_origin: &Point3,
    axis_dir: &Vec3,
    sweep: f64,
    tolerance_linear: f64,
) -> Option<RevolvedSurface> {
    let axis = axis_dir.normalized()?;
    let gate = tolerance_linear * 10.0;
    let (t0, t1) = span.parameter_range();
    let (p_a, p_b) = (span.point_at(t0), span.point_at(t1));
    let (z_a, r_a) = axial_radial(&p_a, axis_origin, &axis);
    let (z_b, r_b) = axial_radial(&p_b, axis_origin, &axis);

    if let Some((_, _, length)) = recognize_line(span) {
        if length < gate {
            return None;
        }
        if (r_a - r_b).abs() < gate && r_a > gate {
            // Cylinder: u = azimuth in orthonormal_basis(axis), v = axial.
            let geometry = SurfaceGeometry::Cylinder {
                origin: *axis_origin,
                axis,
                radius: 0.5 * (r_a + r_b),
            };
            if !surface_agrees(nurbs, &geometry, gate) {
                return None;
            }
            let basis = orthonormal_basis(axis);
            let phi0 = azimuth(&p_a, axis_origin, &axis, basis);
            let sense = azimuth_sense(nurbs, axis_origin, &axis, basis);
            return Some(RevolvedSurface {
                geometry,
                chart: RevolutionChart::Angular {
                    phi0,
                    sense,
                    sweep,
                    v_a: z_a,
                    v_b: z_b,
                },
            });
        }
        if (z_a - z_b).abs() < gate {
            // Plane ⊥ axis through the span's height.
            let origin = *axis_origin + axis * (0.5 * (z_a + z_b));
            let (u_axis, v_axis) = orthonormal_basis(axis);
            let geometry = SurfaceGeometry::Plane {
                origin,
                normal: axis,
                u_axis,
                v_axis,
            };
            if !surface_agrees(nurbs, &geometry, gate) {
                return None;
            }
            let basis = (u_axis, v_axis);
            let phi0 = azimuth(&p_a, axis_origin, &axis, basis);
            let sense = azimuth_sense(nurbs, axis_origin, &axis, basis);
            return Some(RevolvedSurface {
                geometry,
                chart: RevolutionChart::Planar {
                    center: (0.0, 0.0),
                    phi0,
                    sense,
                    r_a,
                    r_b,
                },
            });
        }
        // Cone: apex where the line meets the axis; the cone axis points
        // from the apex towards the widening side.
        let dz = z_b - z_a;
        let dr = r_b - r_a;
        let z_apex = z_a - r_a * dz / dr;
        let apex = *axis_origin + axis * z_apex;
        let cone_axis = if dz * dr > 0.0 { axis } else { -axis };
        let half_angle = (dr / dz).abs().atan();
        let geometry = SurfaceGeometry::Cone {
            apex,
            axis: cone_axis,
            half_angle,
        };
        if !surface_agrees(nurbs, &geometry, gate) {
            return None;
        }
        let basis = orthonormal_basis(cone_axis);
        let phi0 = azimuth(&p_a, &apex, &cone_axis, basis);
        let sense = azimuth_sense(nurbs, &apex, &cone_axis, basis);
        let v_of = |p: &Point3| (*p - apex).dot(&cone_axis);
        return Some(RevolvedSurface {
            geometry,
            chart: RevolutionChart::Angular {
                phi0,
                sense,
                sweep,
                v_a: v_of(&p_a),
                v_b: v_of(&p_b),
            },
        });
    }

    if let Some((center, _normal, radius)) = recognize_arc(span) {
        let (z_c, r_c) = axial_radial(&center, axis_origin, &axis);
        if r_c < gate {
            // Sphere: the chart is fixed to world z, so only an axis
            // along z keeps the seam a meridian.
            if axis.cross(&Vec3::new(0.0, 0.0, 1.0)).length() > 1.0e-9 {
                return None;
            }
            let geometry = SurfaceGeometry::Sphere { center, radius };
            if !surface_agrees(nurbs, &geometry, gate) {
                return None;
            }
            let polar = |p: &Point3| {
                let d = *p - center;
                (d.z / radius).clamp(-1.0, 1.0).acos()
            };
            let world_z = Vec3::new(0.0, 0.0, 1.0);
            let basis = (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
            let phi0 = azimuth(&nurbs.point_at(0.5, 0.0), &center, &world_z, basis);
            let sense = azimuth_sense(nurbs, &center, &world_z, basis);
            return Some(RevolvedSurface {
                geometry,
                chart: RevolutionChart::Angular {
                    phi0,
                    sense,
                    sweep,
                    v_a: polar(&p_a),
                    v_b: polar(&p_b),
                },
            });
        }
        if r_c > radius + gate {
            // Torus: minor angle v measured from the outer equator
            // (v = 0 at +radial, π/2 towards +axis).
            let torus_center = *axis_origin + axis * z_c;
            let geometry = SurfaceGeometry::Torus {
                center: torus_center,
                axis,
                major_radius: r_c,
                minor_radius: radius,
            };
            if !surface_agrees(nurbs, &geometry, gate) {
                return None;
            }
            let basis = orthonormal_basis(axis);
            let phi0 = azimuth(&p_a, axis_origin, &axis, basis);
            let sense = azimuth_sense(nurbs, axis_origin, &axis, basis);
            let minor = |p: &Point3| {
                let (z, r) = axial_radial(p, axis_origin, &axis);
                (z - z_c).atan2(r - r_c)
            };
            // Follow the span's direction through its interior so v_b
            // sits on the correct side of the ±π cut. A FULL circle
            // profile (a torus from a circle sketch) returns to v_a,
            // and two half-turn steps through its midpoint were each
            // ambiguous at exactly ±π and cancelled: the tube seam's
            // pcurve collapsed to a point and a box subtraction found
            // both of its intersection circles "outside" the face
            // (PK-418). Eight steps of at most a quarter turn are
            // unambiguous for any span up to a full circle.
            let v_a = minor(&p_a);
            let steps = 8;
            let mut v_prev = v_a;
            let mut v_b = v_a;
            for k in 1..=steps {
                let t = t0 + (t1 - t0) * k as f64 / steps as f64;
                let v_k = minor(&span.point_at(t));
                v_b += angle_step(v_prev, v_k);
                v_prev = v_k;
            }
            return Some(RevolvedSurface {
                geometry,
                chart: RevolutionChart::Angular {
                    phi0,
                    sense,
                    sweep,
                    v_a,
                    v_b,
                },
            });
        }
    }
    None
}

/// A circular arc span (full or partial): `(center, normal, radius)`.
fn recognize_arc(span: &NurbsCurve) -> Option<(Point3, Vec3, f64)> {
    if let Some(found) = recognize_circle(span) {
        return Some(found);
    }
    // Partial arc: fit a circle through three points and verify.
    let (t0, t1) = span.parameter_range();
    let sample = |f: f64| span.point_at(t0 + (t1 - t0) * f);
    let (a, b, c) = (sample(0.0), sample(0.5), sample(1.0));
    let ab = b - a;
    let ac = c - a;
    let normal = ab.cross(&ac).normalized()?;
    // Circumcentre of a, b, c.
    let ab2 = ab.dot(&ab);
    let ac2 = ac.dot(&ac);
    let cross = ab.cross(&ac);
    let denom = 2.0 * cross.dot(&cross);
    if denom < 1.0e-18 {
        return None;
    }
    let offset = (ac * ab2 - ab * ac2).cross(&cross) * (1.0 / denom);
    let center = a + offset;
    let radius = (a - center).length();
    for i in 0..=12 {
        let p = sample(i as f64 / 12.0);
        let d = p - center;
        if ((d.length() - radius).abs() > 1.0e-9 * radius.max(1.0))
            || d.dot(&normal).abs() > 1.0e-9 * radius.max(1.0)
        {
            return None;
        }
    }
    Some((center, normal, radius))
}
