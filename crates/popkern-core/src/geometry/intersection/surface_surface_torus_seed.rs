//! Closest-approach seeding fallback for torus pairs (PK-310;
//! extracted from `surface_surface.rs` under the file-size ratchet).
//!
//! A shallow outer-outer graze intersects in a lens narrower than the
//! patch scan's refinement floor, and an empty seed set let the
//! boolean serve BOTH OPERANDS OVERLAID UNCUT with validate_all
//! silent (the self-intersection probe consults the same blind
//! seeding). The tube-centre circles' closest approach is cheap and
//! exact per probe; when the tubes penetrate beyond tolerance, Newton
//! from the closest pair recovers real seeds.

use crate::geometry::intersection::marching::{MarchSample, march_project};
use crate::geometry::intersection::starting_points::StartingPoint;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Push closest-approach seeds when the scan found none and the two
/// tori's tubes measurably penetrate.
pub(super) fn push_torus_closest_approach_seeds(
    starting_points: &mut Vec<StartingPoint>,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    tolerance: &Tolerance,
) {
    if !starting_points.is_empty() {
        return;
    }
    let (
        SurfaceGeometry::Torus {
            minor_radius: r1, ..
        },
        SurfaceGeometry::Torus {
            minor_radius: r2, ..
        },
    ) = (surface1, surface2)
    else {
        return;
    };
    let Some((u1c, u2c, centre_dist)) = torus_tube_centre_closest_approach(surface1, surface2)
    else {
        return;
    };
    if centre_dist >= r1 + r2 - tolerance.linear() * 10.0 {
        return;
    }
    let p1 = surface1.point_at(u1c, 0.0);
    let p2 = surface2.point_at(u2c, 0.0);
    let mid = Point3::new(
        0.5 * (p1.x + p2.x),
        0.5 * (p1.y + p2.y),
        0.5 * (p1.z + p2.z),
    );
    let n1 = surface1.normal_at(u1c, 0.0);
    let n2 = surface2.normal_at(u2c, 0.0);
    let dir = {
        let c = n1.cross(&n2);
        c.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
    };
    // Two seeds nudged to either side along the curve tangent so the
    // marcher can trace both lens loops.
    for sign in [1.0, -1.0] {
        let target = mid + dir * (tolerance.linear() * 100.0 * sign);
        if let Some(sample) = march_project(
            surface1, u1c, 0.0, surface2, u2c, 0.0, &target, &dir, tolerance,
        ) {
            let _: &MarchSample = &sample;
            starting_points.push(StartingPoint {
                sample,
                is_tangent: false,
            });
        }
    }
}

/// Closest approach between two tori's tube-centre circles: returns
/// `(u1, u2, distance)` — coarse scan + golden-section refinement.
fn torus_tube_centre_closest_approach(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
) -> Option<(f64, f64, f64)> {
    let centre_at = |s: &SurfaceGeometry, u: f64| -> Option<Point3> {
        if let SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            ..
        } = s
        {
            let (p1, p2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
            Some(*center + (p1 * u.pcos() + p2 * u.psin()) * *major_radius)
        } else {
            None
        }
    };
    let tau = std::f64::consts::TAU;
    let n = 64;
    let mut best = (f64::MAX, 0.0, 0.0);
    for i in 0..n {
        let u1 = tau * i as f64 / n as f64;
        let c1 = centre_at(surface1, u1)?;
        for j in 0..n {
            let u2 = tau * j as f64 / n as f64;
            let c2 = centre_at(surface2, u2)?;
            let d = c1.distance_to(&c2);
            if d < best.0 {
                best = (d, u1, u2);
            }
        }
    }
    // Coordinate-wise golden refinement.
    let (mut u1, mut u2) = (best.1, best.2);
    let h = tau / n as f64;
    for _ in 0..4 {
        for which in 0..2 {
            let eval = |x: f64| -> f64 {
                let (a, b) = if which == 0 { (x, u2) } else { (u1, x) };
                match (centre_at(surface1, a), centre_at(surface2, b)) {
                    (Some(c1), Some(c2)) => c1.distance_to(&c2),
                    _ => f64::MAX,
                }
            };
            let x0 = if which == 0 { u1 } else { u2 };
            let (mut lo, mut hi) = (x0 - h, x0 + h);
            for _ in 0..30 {
                let m1 = lo + (hi - lo) / 3.0;
                let m2 = hi - (hi - lo) / 3.0;
                if eval(m1) < eval(m2) {
                    hi = m2;
                } else {
                    lo = m1;
                }
            }
            if which == 0 {
                u1 = 0.5 * (lo + hi);
            } else {
                u2 = 0.5 * (lo + hi);
            }
        }
    }
    let d = centre_at(surface1, u1)?.distance_to(&centre_at(surface2, u2)?);
    Some((u1, u2, d))
}
