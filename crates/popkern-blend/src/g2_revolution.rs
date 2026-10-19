//! Exact G2 fillet surfaces for circular spines: a surface of
//! revolution of ONE degree-5 G2 cross-section (PK-389).
//!
//! The generic G2 path lofts a cubic B-spline through per-sample
//! cross-sections along the spine. For a coaxial plane × {cylinder,
//! cone, sphere} pair the spine is a circle and every section is the
//! same meridian curve rotated about the axis; the loft then
//! approximates a circle with a polynomial and its rails leave the
//! exact circular edges by up to 2.6e-4 (measured on the flagship boss
//! plate, r = 0.08: the served body failed the serving-boundary
//! curve–surface check while both rails were exact to 1e-16).
//!
//! This arm revolves the first section with rational-quadratic arcs
//! (≤ 90° per segment, exact), so the surface rows ARE the rails and
//! the curvature-continuous profile is preserved everywhere. It is
//! taken only when a certificate holds: every contact solution's own
//! G2 section equals the first section rotated to that solution's
//! station within a scale-relative tolerance. Any other configuration
//! returns `None` and the loft path stands.

use crate::contact::ContactSolution;
use crate::g2_cross_section::{G2Request, g2_cross_section_bezier};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
use popkern_math::PortableMath;

/// Rotation axis of a circular-spine pair: a point on the axis and
/// the unit direction. `None` when the pair is not coaxial. Beyond the
/// plane arms, any two supports of revolution on one axis qualify
/// (PK-433: a revolved groove's torus beside its wall cylinder — the
/// loft's rails sagged 2.6e-5 off the exact contact circles between
/// stations, and the outer race's G2 strip read 7.4e-6 of volume).
fn revolution_axis(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
) -> Option<(Point3, Vec3)> {
    let coaxial = |normal: &Vec3, axis: &Vec3| -> Option<Vec3> {
        let n = normal.normalized()?;
        let a = axis.normalized()?;
        (n.dot(&a).abs() > 1.0 - 1.0e-9).then_some(a)
    };
    match (surface1, surface2) {
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Cylinder { origin, axis, .. })
        | (SurfaceGeometry::Cylinder { origin, axis, .. }, SurfaceGeometry::Plane { normal, .. }) => {
            Some((*origin, coaxial(normal, axis)?))
        }
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Cone { apex, axis, .. })
        | (SurfaceGeometry::Cone { apex, axis, .. }, SurfaceGeometry::Plane { normal, .. }) => {
            Some((*apex, coaxial(normal, axis)?))
        }
        (SurfaceGeometry::Plane { normal, .. }, SurfaceGeometry::Sphere { center, .. })
        | (SurfaceGeometry::Sphere { center, .. }, SurfaceGeometry::Plane { normal, .. }) => {
            Some((*center, normal.normalized()?))
        }
        _ => crate::fillet_surface_coaxial::shared_axis(surface1, surface2),
    }
}

/// Meridian coordinates of a point about the axis: (axial height,
/// radial distance, unit radial direction).
fn meridian(p: &Point3, origin: &Point3, axis: &Vec3) -> (f64, f64, Option<Vec3>) {
    let r = *p - *origin;
    let h = r.dot(axis);
    let radial = r - *axis * h;
    let rho = radial.length();
    let dir = (rho > 0.0).then(|| radial * (1.0 / rho));
    (h, rho, dir)
}

/// Signed angle of `dir` in the (e1, e2) frame.
fn angle_of(dir: &Vec3, e1: &Vec3, e2: &Vec3) -> f64 {
    dir.dot(e2).patan2(dir.dot(e1))
}

/// Build the exact revolved G2 surface, or `None` when the pair or
/// the contact solutions do not certify as a circular spine.
pub(crate) fn g2_revolution_from_contact_solutions(
    solutions: &[ContactSolution],
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
) -> Option<NurbsSurface> {
    if solutions.len() < 2 {
        return None;
    }
    let (origin, axis) = revolution_axis(surface1, surface2)?;
    let s1: &dyn SurfaceEval = surface1;
    let s2: &dyn SurfaceEval = surface2;
    let section = |solution: &ContactSolution| {
        g2_cross_section_bezier(&G2Request {
            surface_0: s1,
            surface_1: s2,
            solution,
        })
        .ok()
    };
    let first = section(&solutions[0])?;
    let cps = &first.curve.control_points;
    if cps.len() != 6 {
        return None;
    }
    // The meridian frame: e1 through the first ball centre.
    let (_, rho_c, dir_c) = meridian(&solutions[0].center, &origin, &axis);
    let e1 = dir_c?;
    let e2_ccw = axis.cross(&e1);
    let scale = rho_c.max(first.chord_length).max(1.0e-12);
    let coord_tol = 1.0e-9 * scale;
    // Every control point of the first section lies in the meridian
    // half-plane through e1.
    let profile: Vec<(f64, f64)> = cps
        .iter()
        .map(|p| {
            let (h, rho, dir) = meridian(p, &origin, &axis);
            match dir {
                Some(d) if rho > coord_tol && (d - e1).length() > 1.0e-7 => None,
                _ => Some((h, rho)),
            }
        })
        .collect::<Option<Vec<_>>>()?;
    // Rotation sense from the second station; stations unwrap
    // monotonically in that sense.
    let (_, _, dir_1) = meridian(&solutions[1].center, &origin, &axis);
    let ccw = angle_of(&dir_1?, &e1, &e2_ccw) >= 0.0;
    let e2 = if ccw { e2_ccw } else { e2_ccw * -1.0 };
    let mut sweep = 0.0f64;
    let mut prev = 0.0f64;
    for solution in solutions.iter().skip(1) {
        let (_, _, dir) = meridian(&solution.center, &origin, &axis);
        let mut theta = angle_of(&dir?, &e1, &e2);
        while theta < prev - 1.0e-9 {
            theta += std::f64::consts::TAU;
        }
        sweep = theta;
        prev = theta;
        // Certificate: this station's own G2 section is the first
        // section rotated here.
        let own = section(solution)?;
        if own.curve.control_points.len() != 6 {
            return None;
        }
        let (c, s) = (theta.pcos(), theta.psin());
        for (k, (h, rho)) in profile.iter().enumerate() {
            let expected = origin + axis * *h + (e1 * c + e2 * s) * *rho;
            if own.curve.control_points[k].distance_to(&expected) > 1.0e-7 * scale {
                return None;
            }
        }
    }
    if sweep <= 1.0e-9 {
        return None;
    }
    // Closed spine: the last station returns to the first.
    let closed = solutions[0]
        .center
        .distance_to(&solutions[solutions.len() - 1].center)
        < 1.0e-7 * scale;
    let sweep = if closed { std::f64::consts::TAU } else { sweep };
    // Rational-quadratic arcs, at most 90° each.
    let segments = ((sweep / std::f64::consts::FRAC_PI_2) - 1.0e-9)
        .ceil()
        .max(1.0) as usize;
    let delta = sweep / segments as f64;
    let w_mid = (delta / 2.0).pcos();
    let num_u = 2 * segments + 1;
    let num_v = 6;
    let mut control_points = Vec::with_capacity(num_u * num_v);
    let mut weights = Vec::with_capacity(num_u * num_v);
    let place = |theta: f64, h: f64, rho: f64| {
        origin + axis * h + (e1 * theta.pcos() + e2 * theta.psin()) * rho
    };
    for iu in 0..num_u {
        let seg = iu / 2;
        let is_mid = iu % 2 == 1;
        let theta = if is_mid {
            (seg as f64 + 0.5) * delta
        } else {
            seg as f64 * delta
        };
        let (radial_scale, w) = if is_mid {
            (1.0 / w_mid, w_mid)
        } else {
            (1.0, 1.0)
        };
        for (h, rho) in &profile {
            control_points.push(place(theta, *h, *rho * radial_scale));
            weights.push(w);
        }
    }
    let mut knots_u = vec![0.0, 0.0, 0.0];
    for seg in 1..segments {
        let t = seg as f64 / segments as f64;
        knots_u.push(t);
        knots_u.push(t);
    }
    knots_u.extend([1.0, 1.0, 1.0]);
    let surface = NurbsSurface {
        degree_u: 2,
        degree_v: 5,
        knots_u,
        knots_v: first.curve.knots.clone(),
        control_points,
        num_u,
        num_v,
        weights,
    };
    // A closed spine's strip is a ring, smooth across its seam (the
    // section rotated by the exact circle): it is written in periodic
    // form, a chart with no seam boundary like a torus (PK-893). The
    // seam row sits at the first station to the coordinate tolerance;
    // past it the clamped ring stands, as before.
    if closed && let Some(ring) = surface.periodic_u_from_closed_clamped(coord_tol) {
        return Some(ring);
    }
    Some(surface)
}

#[cfg(test)]
#[path = "g2_revolution_tests.rs"]
mod tests;
