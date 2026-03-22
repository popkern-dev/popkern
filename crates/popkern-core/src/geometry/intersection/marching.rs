//! Differential marching for surface-surface intersection curves.

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// A sample point on a surface-surface intersection curve.
#[derive(Debug, Clone)]
pub struct MarchSample {
    /// Parameter on surface 1, u direction.
    pub u1: f64,
    /// Parameter on surface 1, v direction.
    pub v1: f64,
    /// Parameter on surface 2, u direction.
    pub u2: f64,
    /// Parameter on surface 2, v direction.
    pub v2: f64,
    /// The 3D point (equal to surface1.point_at(u1, v1) at convergence).
    pub point: Point3,
}

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

        if f.length() < tolerance.linear && f4.abs() < tolerance.linear {
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

/// Solve a 4×4 linear system `A * x = b` via Gaussian elimination with partial pivoting.
///
/// Returns `None` if the matrix is singular (pivot smaller than `1e-30`).
#[allow(dead_code, clippy::needless_range_loop)]
fn solve_4x4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    // Forward elimination with partial pivoting.
    for col in 0..4 {
        // Find the row with the largest absolute value in this column.
        let mut max_row = col;
        let mut max_val = a[col][col].abs();
        for row in (col + 1)..4 {
            if a[row][col].abs() > max_val {
                max_val = a[row][col].abs();
                max_row = row;
            }
        }

        if max_val < 1e-30 {
            return None;
        }

        // Swap the pivot row into position.
        if max_row != col {
            a.swap(col, max_row);
            b.swap(col, max_row);
        }

        // Eliminate entries below the pivot.
        for row in (col + 1)..4 {
            let factor = a[row][col] / a[col][col];
            for k in col..4 {
                a[row][k] -= factor * a[col][k];
            }
            b[row] -= factor * b[col];
        }
    }

    // Back substitution.
    let mut x = [0.0_f64; 4];
    for i in (0..4).rev() {
        x[i] = b[i];
        for j in (i + 1)..4 {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }

    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    fn xz_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        }
    }

    fn z5_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 5.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    /// The solver finds the intersection of the XY plane (z=0) and the XZ
    /// plane (y=0), which is the X axis.  The advancing-front constraint
    /// `(p - target) · X = 0` fixes x along the march, pinning the solution
    /// to the unique point on the X axis with x=0.5, y=0, z=0.
    #[test]
    fn project_onto_two_perpendicular_planes() {
        let s1 = xy_plane(); // z = 0
        let s2 = xz_plane(); // y = 0

        // Use march_dir = X so the front constraint has non-zero Jacobian entries
        // for the XY plane (whose partials are in the XY plane).
        let march_dir = Vec3::new(1.0, 0.0, 0.0);
        let target = Point3::new(0.5, 0.0, 0.0);

        let result = march_project(
            &s1,
            0.1,
            0.1,
            &s2,
            0.1,
            0.4,
            &target,
            &march_dir,
            &Tolerance::default(),
        );

        assert!(result.is_some(), "expected convergence");
        let s = result.unwrap();
        assert!(s.point.y.abs() < 1e-8, "point should lie on XZ plane (y=0)");
        assert!(s.point.z.abs() < 1e-8, "point should lie on XY plane (z=0)");
    }

    /// Two parallel planes (z=0 and z=5) never intersect; the Jacobian for
    /// the surface-difference rows is singular.
    #[test]
    fn fails_for_parallel_planes() {
        let s1 = xy_plane(); // z = 0
        let s2 = z5_plane(); // z = 5

        let result = march_project(
            &s1,
            0.0,
            0.0,
            &s2,
            0.0,
            0.0,
            &Point3::new(0.0, 0.0, 2.5),
            &Vec3::new(1.0, 0.0, 0.0),
            &Tolerance::default(),
        );

        assert!(result.is_none(), "parallel planes should not converge");
    }

    /// Basic sanity check: solving against the identity matrix returns the
    /// right-hand side unchanged.
    #[test]
    fn solve_4x4_identity() {
        let a = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let b = [1.0, 2.0, 3.0, 4.0];
        let x = solve_4x4(a, b).expect("identity system must be solvable");
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 2.0).abs() < 1e-10);
        assert!((x[2] - 3.0).abs() < 1e-10);
        assert!((x[3] - 4.0).abs() < 1e-10);
    }

    /// A truly singular matrix must return None.
    #[test]
    fn solve_4x4_singular_returns_none() {
        let a = [
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 6.0, 8.0], // row 1 = 2 * row 0
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ];
        let b = [1.0, 2.0, 3.0, 4.0];
        assert!(solve_4x4(a, b).is_none());
    }

    /// Verify that the solver produces a correct answer for a known 4×4 system.
    #[test]
    fn solve_4x4_known_solution() {
        // System: x0=1, x1=2, x2=3, x3=4 encoded as a non-trivial matrix.
        let a = [
            [2.0, 1.0, 0.0, 0.0],
            [1.0, 3.0, 1.0, 0.0],
            [0.0, 1.0, 4.0, 1.0],
            [0.0, 0.0, 1.0, 3.0],
        ];
        // b = A * [1, 2, 3, 4]^T
        let b = [
            2.0 * 1.0 + 1.0 * 2.0,
            1.0 * 1.0 + 3.0 * 2.0 + 1.0 * 3.0,
            1.0 * 2.0 + 4.0 * 3.0 + 1.0 * 4.0,
            1.0 * 3.0 + 3.0 * 4.0,
        ];
        let x = solve_4x4(a, b).expect("system must be solvable");
        assert!((x[0] - 1.0).abs() < 1e-10, "x[0] = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-10, "x[1] = {}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-10, "x[2] = {}", x[2]);
        assert!((x[3] - 4.0).abs() < 1e-10, "x[3] = {}", x[3]);
    }
}
