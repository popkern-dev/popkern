//! Damped boundary-row fitting for KP20 support-cell reconciliation.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

use super::BoundarySide;

#[derive(Debug, Clone, Copy)]
pub(super) struct BoundaryRowFitSample {
    pub(super) parameter: f64,
    pub(super) target: Point3,
}

pub(super) fn eval_boundary(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    t: f64,
) -> Point3 {
    match side {
        BoundarySide::UMin => patch.point_at(0.0, t),
        BoundarySide::UMax => patch.point_at(1.0, t),
        BoundarySide::VMin => patch.point_at(t, 0.0),
        BoundarySide::VMax => patch.point_at(t, 1.0),
    }
}

pub(super) fn max_boundary_displacement(
    before: &RectangularGregoryBiNonic,
    after: &RectangularGregoryBiNonic,
    side: BoundarySide,
) -> f64 {
    let mut worst = 0.0_f64;
    for index in 0..10 {
        let t = index as f64 / 9.0;
        worst =
            worst.max(eval_boundary(before, side, t).distance_to(&eval_boundary(after, side, t)));
    }
    worst
}

pub(super) fn fit_boundary_samples(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[BoundaryRowFitSample],
) {
    fit_boundary_samples_with_regularization(patch, side, samples, 1.0e-4);
}

pub(super) fn fit_boundary_samples_with_regularization(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[BoundaryRowFitSample],
    regularization: f64,
) {
    if samples.is_empty() {
        return;
    }
    let prior = boundary_controls(patch, side);
    let lambda = regularization.max(0.0);
    let mut normal = [[0.0; 10]; 10];
    let mut rhs_x = [0.0; 10];
    let mut rhs_y = [0.0; 10];
    let mut rhs_z = [0.0; 10];
    for sample in samples {
        let basis = bernstein9(sample.parameter);
        for i in 0..10 {
            rhs_x[i] += basis[i] * sample.target.x;
            rhs_y[i] += basis[i] * sample.target.y;
            rhs_z[i] += basis[i] * sample.target.z;
            for j in 0..10 {
                normal[i][j] += basis[i] * basis[j];
            }
        }
    }
    for i in 0..10 {
        normal[i][i] += lambda;
        rhs_x[i] += lambda * prior[i].x;
        rhs_y[i] += lambda * prior[i].y;
        rhs_z[i] += lambda * prior[i].z;
    }
    let xs = solve_10(normal, rhs_x);
    let ys = solve_10(normal, rhs_y);
    let zs = solve_10(normal, rhs_z);
    for index in 0..10 {
        set_boundary_control(
            patch,
            side,
            index,
            Point3::new(xs[index], ys[index], zs[index]),
        );
    }
}

fn boundary_controls(patch: &RectangularGregoryBiNonic, side: BoundarySide) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = match side {
            BoundarySide::UMin => patch.grid[0][index],
            BoundarySide::UMax => patch.grid[9][index],
            BoundarySide::VMin => patch.grid[index][0],
            BoundarySide::VMax => patch.grid[index][9],
        };
    }
    out
}

fn set_boundary_control(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    index: usize,
    point: Point3,
) {
    match side {
        BoundarySide::UMin => patch.grid[0][index] = point,
        BoundarySide::UMax => patch.grid[9][index] = point,
        BoundarySide::VMin => patch.grid[index][0] = point,
        BoundarySide::VMax => patch.grid[index][9] = point,
    }
}

pub(super) fn bernstein9(t: f64) -> [f64; 10] {
    let omt = 1.0 - t;
    let mut out = [0.0; 10];
    for (i, value) in out.iter_mut().enumerate() {
        *value = binomial9(i) * t.powi(i as i32) * omt.powi((9 - i) as i32);
    }
    out
}

pub(super) fn eval_row(row: &[Point3; 10], t: f64) -> Point3 {
    let basis = bernstein9(t);
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for (index, weight) in basis.iter().enumerate() {
        x += row[index].x * weight;
        y += row[index].y * weight;
        z += row[index].z * weight;
    }
    Point3::new(x, y, z)
}

fn binomial9(i: usize) -> f64 {
    [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0][i]
}

#[allow(clippy::needless_range_loop)] // Gaussian elimination is explicit row/column indexing.
fn solve_10(mut matrix: [[f64; 10]; 10], mut rhs: [f64; 10]) -> [f64; 10] {
    for pivot in 0..10 {
        let mut best = pivot;
        for row in pivot + 1..10 {
            if matrix[row][pivot].abs() > matrix[best][pivot].abs() {
                best = row;
            }
        }
        matrix.swap(pivot, best);
        rhs.swap(pivot, best);
        let scale = matrix[pivot][pivot];
        if scale.abs() <= 1.0e-14 {
            continue;
        }
        for column in pivot..10 {
            matrix[pivot][column] /= scale;
        }
        rhs[pivot] /= scale;
        for row in 0..10 {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            for column in pivot..10 {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    rhs
}
