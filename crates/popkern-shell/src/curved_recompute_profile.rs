//! Coaxial-profile recompute arm (PK-324): transversal intersection
//! circles of two revolution surfaces about a common axis, found as
//! the axial station where the profile radius functions agree (1D
//! Newton seeded from the clone circle's own station). Split from
//! `curved_recompute_similarity.rs` under the 500-line structural cap;
//! shares its verification and similarity mapping.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

use super::{map_similarity_geometry, verify_candidate};

/// Profile radius function of a revolution surface about the common
/// axis `(origin, axis)`: `r(t)` with `t` the axial station. `None`
/// when the surface is not a revolution about that axis or `t` leaves
/// its domain.
struct RevolutionProfile {
    /// Axial station of the surface's own reference point.
    t_ref: f64,
    kind: ProfileKind,
}

enum ProfileKind {
    Cylinder {
        radius: f64,
    },
    /// `dir` is +1 when the cone opens along the common axis, -1 against.
    Cone {
        tan_half: f64,
        dir: f64,
    },
    Sphere {
        radius: f64,
    },
    /// `sheet` is +1 for the outer sheet (r > major), -1 for the inner.
    Torus {
        major: f64,
        minor: f64,
        sheet: f64,
    },
}

impl RevolutionProfile {
    fn radius_at(&self, t: f64) -> Option<f64> {
        match self.kind {
            ProfileKind::Cylinder { radius } => Some(radius),
            ProfileKind::Cone { tan_half, dir } => {
                let s = (t - self.t_ref) * dir;
                (s > 0.0).then_some(s * tan_half)
            }
            ProfileKind::Sphere { radius } => {
                let d = t - self.t_ref;
                let sq = radius * radius - d * d;
                (sq > 0.0).then(|| sq.sqrt())
            }
            ProfileKind::Torus {
                major,
                minor,
                sheet,
            } => {
                let d = t - self.t_ref;
                let sq = minor * minor - d * d;
                (sq > 0.0).then(|| major + sheet * sq.sqrt())
            }
        }
    }
}

/// Build the profile of `surface` about the axis `(origin, axis)`,
/// requiring the surface's own axis (or center) to lie on it. `r0` is
/// the clone circle's radius, used only to pick the torus sheet.
fn revolution_profile(
    surface: &SurfaceGeometry,
    origin: Point3,
    axis: Vec3,
    r0: f64,
    axis_eps: f64,
) -> Option<RevolutionProfile> {
    let on_axis = |p: &Point3| -> Option<f64> {
        let rel = *p - origin;
        let t = rel.dot(&axis);
        ((rel - axis * t).length() <= axis_eps).then_some(t)
    };
    let parallel = |a: &Vec3| -> Option<f64> {
        let a = a.normalized()?;
        let c = a.cross(&axis).length();
        (c <= 1.0e-9).then(|| a.dot(&axis).signum())
    };
    match surface {
        SurfaceGeometry::Cylinder {
            origin: o,
            axis: a,
            radius,
        } => {
            parallel(a)?;
            Some(RevolutionProfile {
                t_ref: on_axis(o)?,
                kind: ProfileKind::Cylinder { radius: *radius },
            })
        }
        SurfaceGeometry::Cone {
            apex,
            axis: a,
            half_angle,
        } => {
            let dir = parallel(a)?;
            Some(RevolutionProfile {
                t_ref: on_axis(apex)?,
                kind: ProfileKind::Cone {
                    tan_half: half_angle.tan(),
                    dir,
                },
            })
        }
        SurfaceGeometry::Sphere { center, radius } => Some(RevolutionProfile {
            t_ref: on_axis(center)?,
            kind: ProfileKind::Sphere { radius: *radius },
        }),
        SurfaceGeometry::Torus {
            center,
            axis: a,
            major_radius,
            minor_radius,
        } => {
            parallel(a)?;
            Some(RevolutionProfile {
                t_ref: on_axis(center)?,
                kind: ProfileKind::Torus {
                    major: *major_radius,
                    minor: *minor_radius,
                    sheet: if r0 >= *major_radius { 1.0 } else { -1.0 },
                },
            })
        }
        _ => None,
    }
}

/// The axis a revolution surface owns, if any (spheres own none — they
/// adopt the partner's).
fn owned_axis(surface: &SurfaceGeometry) -> Option<(Point3, Vec3)> {
    match surface {
        SurfaceGeometry::Cylinder { origin, axis, .. } => Some((*origin, axis.normalized()?)),
        SurfaceGeometry::Cone { apex, axis, .. } => Some((*apex, axis.normalized()?)),
        SurfaceGeometry::Torus { center, axis, .. } => Some((*center, axis.normalized()?)),
        _ => None,
    }
}

/// Coaxial-profile arm: transversal intersection circle of two
/// revolution surfaces about a common axis, found as the axial station
/// where the profile radius functions agree (1D Newton with numeric
/// derivative, seeded from the clone circle's own station — the
/// nearest-branch selection the design mandates).
pub(crate) fn coaxial_profile_edge_curve(
    old: &CurveOnEdge,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let lin = tolerance.linear();
    // Common axis: one surface must own it; the other's own axis (or
    // center) must lie on it. Structural coaxiality guards, not
    // geometric tolerances.
    let (origin, axis) = owned_axis(surf1).or_else(|| owned_axis(surf2))?;
    let axis_eps = 1.0e-9 * (1.0 + origin.x.abs() + origin.y.abs() + origin.z.abs());

    // Clone circle station and radius about the axis, with a 3-sample
    // coherence check (anything else is not an axis-centered circle).
    let (t0, t1) = (old.t_start, old.t_end);
    let samples: Vec<Point3> = [0.0, 0.37, 0.71]
        .iter()
        .map(|f| old.geometry.point_at(t0 + (t1 - t0) * f))
        .collect();
    let heights: Vec<f64> = samples.iter().map(|p| (*p - origin).dot(&axis)).collect();
    let radii: Vec<f64> = samples
        .iter()
        .zip(&heights)
        .map(|(p, h)| (*p - (origin + axis * *h)).length())
        .collect();
    if heights.iter().any(|h| (h - heights[0]).abs() > lin)
        || radii.iter().any(|r| (r - radii[0]).abs() > lin)
        || radii[0] <= lin
    {
        return None;
    }
    let (h0, r0) = (heights[0], radii[0]);

    let p1 = revolution_profile(surf1, origin, axis, r0, axis_eps)?;
    let p2 = revolution_profile(surf2, origin, axis, r0, axis_eps)?;

    // 1D Newton on f(t) = r1(t) − r2(t), numeric derivative, seeded at
    // the clone station.
    let f = |t: f64| -> Option<f64> { Some(p1.radius_at(t)? - p2.radius_at(t)?) };
    let scale = 1.0 + h0.abs() + r0.abs();
    let mut t = h0;
    for _ in 0..64 {
        let ft = f(t)?;
        if ft.abs() < 1.0e-13 * scale {
            break;
        }
        let dt = 1.0e-7 * scale;
        let fp = (f(t + dt)? - f(t - dt)?) / (2.0 * dt);
        if fp.abs() < 1.0e-16 {
            return None;
        }
        let step = ft / fp;
        t -= step.clamp(-0.25 * scale, 0.25 * scale);
    }
    if f(t)?.abs() > 1.0e-10 * scale {
        return None;
    }
    let r_new = p1.radius_at(t)?;
    if r_new <= lin {
        return None;
    }

    // Similarity: scale about the axis + axial translation.
    let c0 = origin + axis * h0;
    let c1 = origin + axis * t;
    let s = r_new / r0;
    let anchor = Vec3::new(c1.x - s * c0.x, c1.y - s * c0.y, c1.z - s * c0.z);
    let geometry = map_similarity_geometry(&old.geometry, anchor, s)?;
    let candidate = CurveOnEdge {
        geometry,
        t_start: t0,
        t_end: t1,
    };
    verify_candidate(&candidate, surf1, surf2, start_pos, end_pos, lin).then_some(candidate)
}
