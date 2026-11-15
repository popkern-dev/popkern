//! The exact offset of a NURBS surface that is a sphere (PK-1215).
//!
//! A sphere written as a NURBS net (a NURBS ball, a revolved semicircle)
//! offsets to the concentric sphere, and since every point of the net's
//! surface lies at the radius from the centre, scaling its control
//! points about the centre by the new radius over the old is that offset
//! exactly: the same knots and weights, the homothety commuting with the
//! rational form. The general path fitted it with C2 cubics across the
//! net's C1 double-knot rows and stopped at 8.6e-4 on a ball of radius 5
//! offset by 2, against a tolerance of 2.5e-7.
//!
//! Recognition is measured: the centre is the least-squares sphere
//! through a grid of samples, and every sample of a denser grid, the
//! knot rows among them, must lie on it within a quarter of the linear
//! tolerance (the approximation path's own stage tolerance). Anything
//! else keeps the general path.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;

/// Grid lines per direction for the fit and the check.
const GRID: usize = 17;

/// The offset of `n` by `distance` (`inward` against its normal), when
/// `n` is a sphere and the offset radius stays positive; `None` keeps
/// the general path.
pub(crate) fn sphere_offset(
    n: &NurbsSurface,
    distance: f64,
    inward: bool,
    tolerance: &Tolerance,
) -> Option<NurbsSurface> {
    let ((u0, u1), (v0, v1)) = n.parameter_range();
    let grid = |k: usize, lo: f64, hi: f64| lo + (hi - lo) * k as f64 / (GRID - 1) as f64;
    let mut samples = Vec::with_capacity(GRID * GRID);
    for i in 0..GRID {
        for j in 0..GRID {
            samples.push(n.point_at(grid(i, u0, u1), grid(j, v0, v1)));
        }
    }
    // Knot rows: where a C1 net could leave a sphere unseen by the grid.
    for &u in distinct_interior(&n.knots_u, n.degree_u).iter() {
        for j in 0..GRID {
            samples.push(n.point_at(u, grid(j, v0, v1)));
        }
    }
    for &v in distinct_interior(&n.knots_v, n.degree_v).iter() {
        for i in 0..GRID {
            samples.push(n.point_at(grid(i, u0, u1), v));
        }
    }
    let centre = least_squares_centre(&samples)?;
    let radius = samples.iter().map(|p| p.distance_to(&centre)).sum::<f64>() / samples.len() as f64;
    let gate = 0.25 * tolerance.linear();
    if radius.is_nan()
        || radius <= gate
        || samples
            .iter()
            .any(|p| (p.distance_to(&centre) - radius).abs() > gate)
    {
        return None;
    }
    // Which way the surface's normal faces the centre, read where the
    // chart is regular (the grid's middle).
    let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
    let normal = n.normal_at(um, vm).normalized()?;
    let radial = (n.point_at(um, vm) - centre).normalized()?;
    let facing = normal.dot(&radial);
    if facing.abs() < 0.5 {
        return None;
    }
    let sign = if inward { -1.0 } else { 1.0 };
    let new_radius = radius + sign * distance * facing.signum();
    if new_radius.is_nan() || new_radius <= gate {
        return None;
    }
    let k = new_radius / radius;
    let mut out = n.clone();
    for p in &mut out.control_points {
        *p = centre + (*p - centre) * k;
    }
    Some(out)
}

/// The distinct interior knots of a clamped vector of `degree`.
fn distinct_interior(knots: &[f64], degree: usize) -> Vec<f64> {
    let (lo, hi) = (knots[degree], knots[knots.len() - 1 - degree]);
    let mut out: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|&k| k > lo && k < hi)
        .collect();
    out.dedup();
    out
}

/// The centre of the least-squares sphere through `points`:
/// `2 (p_i − p_0) · c = |p_i|² − |p_0|²`, solved by its normal equations.
fn least_squares_centre(points: &[Point3]) -> Option<Point3> {
    let p0 = points.first()?;
    let mut ata = [[0.0_f64; 3]; 3];
    let mut atb = [0.0_f64; 3];
    let sq = |p: &Point3| p.x * p.x + p.y * p.y + p.z * p.z;
    for p in &points[1..] {
        let row = [2.0 * (p.x - p0.x), 2.0 * (p.y - p0.y), 2.0 * (p.z - p0.z)];
        let rhs = sq(p) - sq(p0);
        for i in 0..3 {
            for j in 0..3 {
                ata[i][j] += row[i] * row[j];
            }
            atb[i] += row[i] * rhs;
        }
    }
    let det = |m: &[[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(&ata);
    let scale = ata[0][0].abs() + ata[1][1].abs() + ata[2][2].abs();
    if d.is_nan() || d.abs() <= 1e-12 * scale.powi(3) {
        return None;
    }
    // Cramer's rule, one column replaced by the right-hand side.
    let solve = |col: usize| {
        let mut m = ata;
        for (r, row) in m.iter_mut().enumerate() {
            row[col] = atb[r];
        }
        det(&m) / d
    };
    let c = Vec3::new(solve(0), solve(1), solve(2));
    Some(Point3::new(c.x, c.y, c.z))
}
