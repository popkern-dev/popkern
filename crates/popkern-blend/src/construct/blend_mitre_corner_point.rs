//! Mitre point and chart helpers for the mitre-corner settlement.
//!
//! The mitre point is the point near the corner at distance r from
//! every strip axis. It is solved by Gauss–Newton on the residuals
//! `|radial_i(P)| − r`, starting from the stale corner (which sits in
//! the solution's basin: √2·r from each axis on a right-angled
//! corner). The solve is accepted only when every residual meets
//! `pos_tol`. A least-squares compromise is never a mitre point.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use super::StripCylinder;

/// Iteration cap. Gauss–Newton converges quadratically from the
/// corner basin in a handful of steps. The cap only bounds a
/// non-converging (no common point) solve.
const MAX_ITERATIONS: usize = 60;

fn radial(p: Point3, cyl: &StripCylinder) -> Option<(f64, Vec3)> {
    let axis = cyl.axis.normalized()?;
    let d = p - cyl.origin;
    let r = d - axis * d.dot(&axis);
    let len = r.length();
    (len > 1e-300).then(|| (len, r * (1.0 / len)))
}

/// The point near `start` on every cylinder, or `None` when no such
/// point exists within `pos_tol`.
pub(super) fn solve_mitre_point(
    cylinders: &[&StripCylinder],
    start: Point3,
    pos_tol: f64,
) -> Option<Point3> {
    let radius = cylinders.first()?.radius;
    let mut p = start;
    for _ in 0..MAX_ITERATIONS {
        // Normal equations JᵀJ δ = −Jᵀf, rows = radial unit vectors.
        let mut jtj = [[0.0_f64; 3]; 3];
        let mut jtf = [0.0_f64; 3];
        let mut worst = 0.0_f64;
        for cyl in cylinders {
            let (len, n) = radial(p, cyl)?;
            let f = len - cyl.radius;
            worst = worst.max(f.abs());
            let row = [n.x, n.y, n.z];
            for i in 0..3 {
                jtf[i] += row[i] * f;
                for j in 0..3 {
                    jtj[i][j] += row[i] * row[j];
                }
            }
        }
        if worst <= pos_tol * 1e-3 {
            break;
        }
        let delta = super::super::surface_pair_march::solve3(jtj, [-jtf[0], -jtf[1], -jtf[2]])?;
        p = p + Vec3::new(delta[0], delta[1], delta[2]);
        if (p - start).length() > 4.0 * radius {
            return None;
        }
    }
    let on_all = cylinders
        .iter()
        .all(|cyl| radial(p, cyl).is_some_and(|(len, _)| (len - cyl.radius).abs() <= pos_tol));
    on_all.then_some(p)
}

/// Chart point of `p` on a cylinder strip, its angle unwrapped onto
/// the branch of `near` (a chart point of the same face).
pub(super) fn chart_point_near(
    surface: &SurfaceGeometry,
    p: Point3,
    near: (f64, f64),
) -> Option<(f64, f64)> {
    let SurfaceGeometry::Cylinder { origin, axis, .. } = surface else {
        return None;
    };
    let (perp1, perp2) = orthonormal_basis(*axis);
    let d = p - *origin;
    let radial = d - *axis * d.dot(axis);
    let mut u = radial.dot(&perp2).patan2(radial.dot(&perp1));
    while u - near.0 > std::f64::consts::PI {
        u -= std::f64::consts::TAU;
    }
    while near.0 - u > std::f64::consts::PI {
        u += std::f64::consts::TAU;
    }
    Some((u, d.dot(axis)))
}
