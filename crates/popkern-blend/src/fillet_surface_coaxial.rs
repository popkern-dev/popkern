//! Exact torus fillets between two supports of revolution about one
//! axis (PK-433).
//!
//! A cylinder, cone, sphere or torus — or a plane perpendicular to
//! the axis — is swept by revolving its meridian section: a line or a
//! circle in the (radius, height) half-plane. A rolling ball touching
//! two coaxial supports keeps its centre at distance `radius` from
//! both sections, so by symmetry the ball-centre locus is a coaxial
//! CIRCLE and the fillet is exactly a torus: its tube centre is a
//! point where the two sections' ±`radius` offsets (again lines and
//! circles) cross. [`super::fillet_surface_analytic`] already serves
//! plane × {sphere, cylinder, cone}; this arm serves every other
//! coaxial pair — a revolved groove's torus beside its wall cylinder
//! (the ball bearing races) fell to the generic arc sweep, whose
//! closed strip exported to OpenCASCADE as a degenerate 3e-7-area
//! patch while popkern's validators read it clean.
//!
//! Branch selection uses the sampled spine (nearest candidate);
//! acceptance uses the walker's solved ball centres when recorded,
//! else the spine samples, at the same `(major + radius) * 1e-4`
//! gate as the plane arms. A pair that is not coaxial, or whose
//! samples miss every candidate, declines to the generic sweep.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use crate::center_curve::CenterCurveResult;
use crate::fillet_surface::FilletSurfaceResult;

/// A support's meridian section in (radius, height) coordinates.
#[derive(Clone, Copy, Debug)]
enum Section {
    /// Through `p` along unit `d`.
    Line { p: (f64, f64), d: (f64, f64) },
    /// Centre `c`, radius `r`.
    Circle { c: (f64, f64), r: f64 },
}

/// Parallelism / coaxiality gate on unit axes and axis-line offsets
/// relative to the support scale (a revolve authors its supports on
/// one exact axis; anything looser is not the coaxial case).
const COAXIAL_REL: f64 = 1e-9;

/// The exact torus fillet between two coaxial supports of revolution,
/// or `None` (see the module docs).
pub(crate) fn try_coaxial_revolution_torus(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
) -> Option<FilletSurfaceResult> {
    let (origin, axis) = shared_axis(surface1, surface2)?;
    let meridian = |p: Point3| -> (f64, f64) {
        let w = p - origin;
        let h = w.dot(&axis);
        ((w - axis * h).length(), h)
    };
    let s1 = section(surface1, origin, axis)?;
    let s2 = section(surface2, origin, axis)?;
    let mut candidates = Vec::new();
    for o1 in offsets(s1, radius) {
        for o2 in offsets(s2, radius) {
            candidates.extend(crossings(o1, o2).into_iter().filter(|q| q.0 > 0.0));
        }
    }
    let (t0, t1) = center_result.center_curve.parameter_range();
    // Interior 90 %: a closed spine is fitted open and wiggles at its
    // end conditions (see `try_circular_center_torus`).
    let spine: Vec<(f64, f64)> = (0..16)
        .map(|i| {
            let t = t0 + (t1 - t0) * (0.05 + 0.9 * i as f64 / 15.0);
            meridian(center_result.center_curve.point_at(t))
        })
        .collect();
    let deviation = |samples: &[(f64, f64)], q: (f64, f64)| {
        samples
            .iter()
            .map(|s| (s.0 - q.0).phypot(s.1 - q.1))
            .fold(0.0_f64, f64::max)
    };
    let best = candidates
        .into_iter()
        .min_by(|a, b| deviation(&spine, *a).total_cmp(&deviation(&spine, *b)))?;
    let verify: Vec<(f64, f64)> = match center_result.contact_solutions.as_deref() {
        Some(solutions) if solutions.len() >= 3 => {
            solutions.iter().map(|s| meridian(s.center)).collect()
        }
        _ => spine,
    };
    let dev = deviation(&verify, best);
    let gate = (best.0 + radius) * 1.0e-4;
    if popkern_core::trace::trace_enabled("analytic-intersection") {
        eprintln!(
            "    [coaxial-torus] tube centre (r {:.6}, h {:.6}) dev {dev:.3e} gate {gate:.3e}",
            best.0, best.1
        );
    }
    if dev > gate {
        return None;
    }
    let trim_curve_1 =
        crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1).ok()?;
    let trim_curve_2 =
        crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2).ok()?;
    Some(FilletSurfaceResult {
        surface: SurfaceGeometry::Torus {
            center: origin + axis * best.1,
            axis,
            major_radius: best.0,
            minor_radius: radius,
        },
        trim_curve_1,
        trim_curve_2,
        continuity_downgrade: None,
    })
}

/// The axis line `(point, unit direction)` both supports revolve
/// about; `None` unless at least one support has a definite axis
/// (cylinder, cone, torus) and the other is coaxial with it.
pub(crate) fn shared_axis(a: &SurfaceGeometry, b: &SurfaceGeometry) -> Option<(Point3, Vec3)> {
    let definite = |s: &SurfaceGeometry| match s {
        SurfaceGeometry::Cylinder { origin, axis, .. } => Some((*origin, *axis)),
        SurfaceGeometry::Cone { apex, axis, .. } => Some((*apex, *axis)),
        SurfaceGeometry::Torus { center, axis, .. } => Some((*center, *axis)),
        _ => None,
    };
    let (origin, axis) = definite(a).or_else(|| definite(b))?;
    let axis = axis.normalized()?;
    let scale = |p: Point3| 1.0 + (p - origin).length();
    let on_axis = |p: Point3| {
        let w = p - origin;
        (w - axis * w.dot(&axis)).length() <= COAXIAL_REL * scale(p)
    };
    let parallel = |d: &Vec3| {
        d.normalized()
            .is_some_and(|d| 1.0 - d.dot(&axis).abs() <= COAXIAL_REL)
    };
    let agrees = |s: &SurfaceGeometry| match s {
        SurfaceGeometry::Plane { normal, .. } => parallel(normal),
        SurfaceGeometry::Sphere { center, .. } => on_axis(*center),
        _ => definite(s).is_some_and(|(p, d)| parallel(&d) && on_axis(p)),
    };
    (agrees(a) && agrees(b)).then_some((origin, axis))
}

/// The support's meridian section about `(origin, axis)`.
fn section(s: &SurfaceGeometry, origin: Point3, axis: Vec3) -> Option<Section> {
    let h = |p: Point3| (p - origin).dot(&axis);
    Some(match s {
        SurfaceGeometry::Plane { origin: o, .. } => Section::Line {
            p: (0.0, h(*o)),
            d: (1.0, 0.0),
        },
        SurfaceGeometry::Cylinder { radius, .. } => Section::Line {
            p: (*radius, 0.0),
            d: (0.0, 1.0),
        },
        SurfaceGeometry::Cone {
            apex,
            axis: a,
            half_angle,
        } => Section::Line {
            p: (0.0, h(*apex)),
            d: (half_angle.psin(), half_angle.pcos() * a.dot(&axis).signum()),
        },
        SurfaceGeometry::Sphere { center, radius, .. } => Section::Circle {
            c: (0.0, h(*center)),
            r: *radius,
        },
        SurfaceGeometry::Torus {
            center,
            major_radius,
            minor_radius,
            ..
        } => Section::Circle {
            c: (*major_radius, h(*center)),
            r: *minor_radius,
        },
        _ => return None,
    })
}

/// The section's two `±distance` offsets (one when an inward circle
/// offset collapses).
fn offsets(s: Section, distance: f64) -> Vec<Section> {
    match s {
        Section::Line { p, d } => [1.0, -1.0]
            .map(|k| Section::Line {
                p: (p.0 - k * distance * d.1, p.1 + k * distance * d.0),
                d,
            })
            .to_vec(),
        Section::Circle { c, r } => [r + distance, r - distance]
            .into_iter()
            .filter(|r| *r > 0.0)
            .map(|r| Section::Circle { c, r })
            .collect(),
    }
}

/// Where two sections cross (tangent contacts count once).
fn crossings(a: Section, b: Section) -> Vec<(f64, f64)> {
    match (a, b) {
        (Section::Line { p, d }, Section::Line { p: q, d: e }) => {
            let det = d.0 * e.1 - d.1 * e.0;
            if det.abs() < 1e-12 {
                return Vec::new();
            }
            let t = ((q.0 - p.0) * e.1 - (q.1 - p.1) * e.0) / det;
            vec![(p.0 + t * d.0, p.1 + t * d.1)]
        }
        (Section::Line { p, d }, Section::Circle { c, r })
        | (Section::Circle { c, r }, Section::Line { p, d }) => {
            // |p + t d - c|² = r², |d| = 1.
            let w = (p.0 - c.0, p.1 - c.1);
            let b = w.0 * d.0 + w.1 * d.1;
            let disc = b * b - (w.0 * w.0 + w.1 * w.1 - r * r);
            if disc < 0.0 {
                return Vec::new();
            }
            let s = disc.sqrt();
            let mut out = vec![(p.0 + (-b - s) * d.0, p.1 + (-b - s) * d.1)];
            if s > 0.0 {
                out.push((p.0 + (-b + s) * d.0, p.1 + (-b + s) * d.1));
            }
            out
        }
        (Section::Circle { c, r }, Section::Circle { c: k, r: q }) => {
            let (dx, dy) = (k.0 - c.0, k.1 - c.1);
            let dist = dx.phypot(dy);
            if dist < 1e-15 || dist > r + q || dist < (r - q).abs() {
                return Vec::new();
            }
            let along = (r * r - q * q + dist * dist) / (2.0 * dist);
            let across = (r * r - along * along).max(0.0).sqrt();
            let (ux, uy) = (dx / dist, dy / dist);
            let m = (c.0 + along * ux, c.1 + along * uy);
            let mut out = vec![(m.0 - across * uy, m.1 + across * ux)];
            if across > 0.0 {
                out.push((m.0 + across * uy, m.1 - across * ux));
            }
            out
        }
    }
}

#[cfg(test)]
#[path = "fillet_surface_coaxial_tests.rs"]
mod tests;
