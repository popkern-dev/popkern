//! Damped boundary-row fitting for KP20 support-cell reconciliation.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;

use crate::construct::kp20_polar_support_cell_boundary_trim::BoundarySide;

#[derive(Debug, Clone, Copy)]
pub(super) struct BoundaryRowFitSample {
    pub(super) parameter: f64,
    pub(super) target: Point3,
}

#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub(super) struct BoundaryRowFitConstraint {
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

#[allow(dead_code)]
pub(super) fn fit_boundary_samples_with_constraints(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[BoundaryRowFitSample],
    constraints: &[BoundaryRowFitConstraint],
    regularization: f64,
) -> bool {
    if samples.is_empty() && constraints.is_empty() {
        return true;
    }
    let prior = boundary_controls(patch, side);
    let lambda = regularization.max(0.0);
    let constraint_count = constraints.len();
    let system_size = 10 + constraint_count;
    let mut normal = vec![vec![0.0_f64; system_size]; system_size];
    let mut rhs_x = vec![0.0_f64; system_size];
    let mut rhs_y = vec![0.0_f64; system_size];
    let mut rhs_z = vec![0.0_f64; system_size];
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
    for (constraint_index, constraint) in constraints.iter().enumerate() {
        let row = 10 + constraint_index;
        let basis = bernstein9(constraint.parameter);
        for i in 0..10 {
            normal[row][i] = basis[i];
            normal[i][row] = basis[i];
        }
        rhs_x[row] = constraint.target.x;
        rhs_y[row] = constraint.target.y;
        rhs_z[row] = constraint.target.z;
    }
    let Some(xs) = solve_dense(normal.clone(), rhs_x) else {
        return false;
    };
    let Some(ys) = solve_dense(normal.clone(), rhs_y) else {
        return false;
    };
    let Some(zs) = solve_dense(normal, rhs_z) else {
        return false;
    };
    for index in 0..10 {
        set_boundary_control(
            patch,
            side,
            index,
            Point3::new(xs[index], ys[index], zs[index]),
        );
    }
    true
}

#[allow(dead_code)]
pub(super) fn max_boundary_constraint_residual(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
    constraints: &[BoundaryRowFitConstraint],
) -> f64 {
    constraints.iter().fold(0.0_f64, |worst, constraint| {
        worst.max(eval_boundary(patch, side, constraint.parameter).distance_to(&constraint.target))
    })
}

#[allow(dead_code)]
pub(super) fn snap_boundary_endpoints_to_samples(
    patch: &mut RectangularGregoryBiNonic,
    side: BoundarySide,
    samples: &[BoundaryRowFitSample],
    tolerance: f64,
) {
    for sample in samples {
        if sample.parameter.abs() <= tolerance {
            set_boundary_control(patch, side, 0, sample.target);
        } else if (sample.parameter - 1.0).abs() <= tolerance {
            set_boundary_control(patch, side, 9, sample.target);
        }
    }
}

#[allow(dead_code)]
pub(super) fn scale_boundary_delta_to_displacement(
    before: &RectangularGregoryBiNonic,
    after: &RectangularGregoryBiNonic,
    side: BoundarySide,
    max_displacement: f64,
) -> RectangularGregoryBiNonic {
    let displacement = max_boundary_displacement(before, after, side);
    if displacement <= max_displacement || displacement <= 1.0e-12 {
        return after.clone();
    }
    let scale = (max_displacement / displacement).clamp(0.0, 1.0);
    let before_controls = boundary_controls(before, side);
    let after_controls = boundary_controls(after, side);
    let mut scaled = before.clone();
    for index in 0..10 {
        set_boundary_control(
            &mut scaled,
            side,
            index,
            before_controls[index] + (after_controls[index] - before_controls[index]) * scale,
        );
    }
    scaled
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

#[allow(dead_code)]
#[allow(clippy::needless_range_loop)] // Gaussian elimination is explicit row/column indexing.
fn solve_dense(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    for pivot in 0..n {
        let mut best = pivot;
        for row in pivot + 1..n {
            if matrix[row][pivot].abs() > matrix[best][pivot].abs() {
                best = row;
            }
        }
        matrix.swap(pivot, best);
        rhs.swap(pivot, best);
        let scale = matrix[pivot][pivot];
        if scale.abs() <= 1.0e-13 {
            return None;
        }
        for column in pivot..n {
            matrix[pivot][column] /= scale;
        }
        rhs[pivot] /= scale;
        for row in 0..n {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            if factor == 0.0 {
                continue;
            }
            for column in pivot..n {
                matrix[row][column] -= factor * matrix[pivot][column];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    Some(rhs)
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
