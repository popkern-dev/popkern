//! Differential marching for surface-surface intersection curves.

use crate::geometry::intersection::degeneracy::DegeneracyHandler;
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

/// Like `march_project`, but also returns the number of Newton iterations used.
///
/// Returns `None` if the solver fails to converge.
#[allow(clippy::too_many_arguments)]
fn march_project_with_iters(
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

        if f.length() < tolerance.linear && f4.abs() < tolerance.linear {
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

        let delta = match solve_4x4(j, rhs) {
            Some(d) => d,
            None => return None,
        };

        u1 -= delta[0];
        v1 -= delta[1];
        u2 -= delta[2];
        v2 -= delta[3];
    }

    None
}

/// March along a surface-surface intersection curve from a starting point.
///
/// Marches bidirectionally (+T and -T), merges results.
/// Returns `(samples, is_closed)`.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn march(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    start: &MarchSample,
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> (Vec<MarchSample>, bool) {
    // March in +T direction.
    let forward = march_one_direction(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, 1.0, tolerance,
        handler,
    );

    // March in -T direction.
    let mut backward = march_one_direction(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, -1.0, tolerance,
        handler,
    );

    // Reverse backward so it goes from far-negative to start.
    backward.reverse();

    // Merge: backward (without duplicate start) + forward.
    // backward ends at start, forward starts at start, so drop last of backward.
    if !backward.is_empty() {
        backward.pop(); // remove duplicate start point
    }
    backward.extend(forward);
    let samples = backward;

    // Check if closed: last sample within tolerance of first, and enough steps.
    let is_closed = if samples.len() > 3 {
        let first = &samples[0];
        let last = &samples[samples.len() - 1];
        first.point.distance_to(&last.point) < tolerance.linear
    } else {
        false
    };

    (samples, is_closed)
}

/// March in one direction along the intersection curve.
#[allow(clippy::too_many_arguments)]
fn march_one_direction(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    start: &MarchSample,
    direction_sign: f64,
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> Vec<MarchSample> {
    let max_steps = 10000;
    let epsilon = 1e-12;
    let min_step = tolerance.linear;
    let max_step = tolerance.linear * 1000.0;
    let mut step_size = tolerance.linear * 50.0;

    let mut samples = vec![start.clone()];

    for _ in 0..max_steps {
        let current = &samples[samples.len() - 1];

        // Compute surface normals at current parameters.
        let n1 = surface1.normal_at(current.u1, current.v1);
        let n2 = surface2.normal_at(current.u2, current.v2);

        // Check singularity: degenerate normal.
        if n1.length() < epsilon {
            if handler.on_singularity(&current.point) {
                break;
            }
        }
        if n2.length() < epsilon {
            if handler.on_singularity(&current.point) {
                break;
            }
        }

        // Tangent direction = N1 x N2.
        let t_raw = n1.cross(&n2);
        let t_len = t_raw.length();

        // Check tangency: surfaces are tangent when cross product is near zero.
        if t_len < epsilon {
            let angle = n1.angle_to(&n2);
            if handler.on_tangent_during_march(&current.point, angle) {
                break;
            }
        }

        // Normalize tangent direction and apply direction sign.
        let t_normalized = match t_raw.normalized() {
            Some(t) => t * direction_sign,
            None => break, // cannot normalize zero vector
        };

        // Try Newton projection, with retry on failure.
        let mut succeeded = false;
        let mut retry_step = step_size;
        for _retry in 0..3 {
            let retry_target = current.point + t_normalized * retry_step;
            match march_project_with_iters(
                surface1,
                current.u1,
                current.v1,
                surface2,
                current.u2,
                current.v2,
                &retry_target,
                &t_normalized,
                tolerance,
            ) {
                Some((sample, iters)) => {
                    // Adaptive step based on iteration count.
                    if iters < 3 {
                        step_size = (retry_step * 1.5).min(max_step);
                    } else if iters > 8 {
                        step_size = (retry_step * 0.5).max(min_step);
                    } else {
                        step_size = retry_step;
                    }

                    // Check domain bounds.
                    if sample.u1 < u1_range.0
                        || sample.u1 > u1_range.1
                        || sample.v1 < v1_range.0
                        || sample.v1 > v1_range.1
                        || sample.u2 < u2_range.0
                        || sample.u2 > u2_range.1
                        || sample.v2 < v2_range.0
                        || sample.v2 > v2_range.1
                    {
                        // Exited domain; stop marching in this direction.
                        return samples;
                    }

                    samples.push(sample);
                    succeeded = true;
                    break;
                }
                None => {
                    // Halve step and retry.
                    retry_step *= 0.5;
                    if retry_step < min_step {
                        break;
                    }
                }
            }
        }

        if !succeeded {
            break;
        }
    }

    samples
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
    use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;

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

    /// Two perpendicular planes (z=0 normal and x=0 normal) intersect along the
    /// y-axis.  Marching from the origin should produce samples along y with
    /// x ~ 0 and z ~ 0.
    #[test]
    fn march_two_perpendicular_planes() {
        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(1.0, 0.0, 0.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        };
        let start = MarchSample {
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            point: Point3::new(0.0, 0.0, 0.0),
        };
        let (samples, is_closed) = march(
            &s1,
            (-1.0, 1.0),
            (-1.0, 1.0),
            &s2,
            (-1.0, 1.0),
            (-1.0, 1.0),
            &start,
            &Tolerance::default(),
            &DefaultDegeneracyHandler,
        );
        assert!(samples.len() > 2, "Should produce multiple samples");
        assert!(!is_closed, "Planes produce open intersection");
        // All samples should have x ~ 0 and z ~ 0 (intersection is y-axis)
        for s in &samples {
            assert!(
                (s.point.x).abs() < 1e-4,
                "x should be ~0, got {}",
                s.point.x
            );
            assert!(
                (s.point.z).abs() < 1e-4,
                "z should be ~0, got {}",
                s.point.z
            );
        }
    }
}
