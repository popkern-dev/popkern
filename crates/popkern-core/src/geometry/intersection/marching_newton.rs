//! Newton projection solvers used by the differential marching algorithm.

use super::marching_linalg::solve_4x4;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

use super::MarchSample;

/// Project onto both surfaces using Newton 4×4.
///
/// Solves the system:
/// - `S1(u1,v1) - S2(u2,v2) = 0`  (3 scalar equations)
/// - `(S1(u1,v1) - target) · march_dir = 0`  (advancing front constraint)
///
/// Returns `None` if the solver fails to converge or the Jacobian is singular.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn march_project(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    target: &Point3,
    march_dir: &Vec3,
    tolerance: &Tolerance,
) -> Option<MarchSample> {
    let max_iter = 20;
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;

    for _ in 0..max_iter {
        let p1 = surface1.point_at(u1, v1);
        let p2 = surface2.point_at(u2, v2);

        // Rows 0-2 residual: S1 - S2
        let f = p1 - p2;
        // Row 3 residual: (S1 - target) · march_dir
        let f4 = (p1 - *target).dot(march_dir);

        if crate::adaptive_precision::adaptive_dist_sq_lt(
            (p1.x, p1.y, p1.z),
            (p2.x, p2.y, p2.z),
            tolerance.linear() * tolerance.linear(),
        ) && f4.abs() < tolerance.linear()
        {
            return Some(MarchSample {
                u1,
                v1,
                u2,
                v2,
                point: p1,
            });
        }

        let s1u = surface1.partial_u(u1, v1);
        let s1v = surface1.partial_v(u1, v1);
        let s2u = surface2.partial_u(u2, v2);
        let s2v = surface2.partial_v(u2, v2);

        // Build 4×4 Jacobian J such that J * [du1, dv1, du2, dv2]^T = [f.x, f.y, f.z, f4]:
        //
        // Row 0: [S1_u.x,         S1_v.x,         -S2_u.x,  -S2_v.x ]
        // Row 1: [S1_u.y,         S1_v.y,         -S2_u.y,  -S2_v.y ]
        // Row 2: [S1_u.z,         S1_v.z,         -S2_u.z,  -S2_v.z ]
        // Row 3: [T·S1_u,         T·S1_v,          0,        0       ]
        //
        // Derivation of row 3:
        //   F4 = (S1(u1,v1) - target) · T
        //   ∂F4/∂u1 = S1_u · T,  ∂F4/∂v1 = S1_v · T,  ∂F4/∂u2 = 0,  ∂F4/∂v2 = 0
        let j = [
            [s1u.x, s1v.x, -s2u.x, -s2v.x],
            [s1u.y, s1v.y, -s2u.y, -s2v.y],
            [s1u.z, s1v.z, -s2u.z, -s2v.z],
            [march_dir.dot(&s1u), march_dir.dot(&s1v), 0.0, 0.0],
        ];
        let rhs = [f.x, f.y, f.z, f4];

        let delta = solve_4x4(j, rhs)?;

        u1 -= delta[0];
        v1 -= delta[1];
        u2 -= delta[2];
        v2 -= delta[3];
    }

    None
}

/// Like [`march_project`], but also returns the number of Newton iterations used.
///
/// Returns `None` if the solver fails to converge.
#[allow(clippy::too_many_arguments)]
pub(crate) fn march_project_with_iters(
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
    surface2: &SurfaceGeometry,
    u2: f64,
    v2: f64,
    target: &Point3,
    march_dir: &Vec3,
    tolerance: &Tolerance,
) -> Option<(MarchSample, usize)> {
    let max_iter = 20;
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;

    for iter in 0..max_iter {
        let p1 = surface1.point_at(u1, v1);
        let p2 = surface2.point_at(u2, v2);

        let f = p1 - p2;
        let f4 = (p1 - *target).dot(march_dir);

        if crate::adaptive_precision::adaptive_dist_sq_lt(
            (p1.x, p1.y, p1.z),
            (p2.x, p2.y, p2.z),
            tolerance.linear() * tolerance.linear(),
        ) && f4.abs() < tolerance.linear()
        {
            return Some((
                MarchSample {
                    u1,
                    v1,
                    u2,
                    v2,
                    point: p1,
                },
                iter + 1,
            ));
        }

        let s1u = surface1.partial_u(u1, v1);
        let s1v = surface1.partial_v(u1, v1);
        let s2u = surface2.partial_u(u2, v2);
        let s2v = surface2.partial_v(u2, v2);

        let j = [
            [s1u.x, s1v.x, -s2u.x, -s2v.x],
            [s1u.y, s1v.y, -s2u.y, -s2v.y],
            [s1u.z, s1v.z, -s2u.z, -s2v.z],
            [march_dir.dot(&s1u), march_dir.dot(&s1v), 0.0, 0.0],
        ];
        let rhs = [f.x, f.y, f.z, f4];

        let delta = solve_4x4(j, rhs)?;

        u1 -= delta[0];
        v1 -= delta[1];
        u2 -= delta[2];
        v2 -= delta[3];
    }

    None
}
