//! Cubic/quartic Bezier evaluation, splitting, elevation and the
//! 5x5 linear solve used by the setback-cap candidate fit. Pure
//! geometry math with no BRep or candidate state — split from
//! `pns_setback_candidate.rs`, which `use`s it wholesale.

use popkern_core::geometry::{Point3, Vec3};

/// De Casteljau split of a cubic Bezier at `t`.
pub(crate) fn split_cubic(c: &[Point3; 4], t: f64) -> ([Point3; 4], [Point3; 4]) {
    let l = |a: Point3, b: Point3| {
        Point3::new(
            a.x + (b.x - a.x) * t,
            a.y + (b.y - a.y) * t,
            a.z + (b.z - a.z) * t,
        )
    };
    let q0 = l(c[0], c[1]);
    let q1 = l(c[1], c[2]);
    let q2 = l(c[2], c[3]);
    let r0 = l(q0, q1);
    let r1 = l(q1, q2);
    let s0 = l(r0, r1);
    ([c[0], q0, r0, s0], [s0, r1, q2, c[3]])
}

/// Exact degree elevation of a cubic Bezier to quartic.
pub(crate) fn elevate_cubic_to_quartic(c: &[Point3; 4]) -> [Point3; 5] {
    let lerp = |a: Point3, b: Point3, f: f64| {
        Point3::new(
            a.x * f + b.x * (1.0 - f),
            a.y * f + b.y * (1.0 - f),
            a.z * f + b.z * (1.0 - f),
        )
    };
    [
        c[0],
        lerp(c[0], c[1], 1.0 / 4.0),
        lerp(c[1], c[2], 2.0 / 4.0),
        lerp(c[2], c[3], 3.0 / 4.0),
        c[3],
    ]
}

pub(crate) fn bezier3_eval(cps: &[Point3; 4], t: f64) -> Point3 {
    let mut w = *cps;
    for level in (1..4).rev() {
        for i in 0..level {
            w[i] = Point3::new(
                w[i].x + (w[i + 1].x - w[i].x) * t,
                w[i].y + (w[i + 1].y - w[i].y) * t,
                w[i].z + (w[i + 1].z - w[i].z) * t,
            );
        }
    }
    w[0]
}

pub(crate) fn bezier3_tangent(cps: &[Point3; 4], t: f64) -> Vec3 {
    // Derivative of a cubic Bezier: 3·Σ (P_{i+1}−P_i) B_i^2(t).
    let d: [Vec3; 3] = core::array::from_fn(|i| {
        Vec3::new(
            cps[i + 1].x - cps[i].x,
            cps[i + 1].y - cps[i].y,
            cps[i + 1].z - cps[i].z,
        )
    });
    let u = 1.0 - t;
    (d[0] * (u * u) + d[1] * (2.0 * u * t) + d[2] * (t * t)) * 3.0
}

pub(crate) fn bezier4_eval(cps: &[Point3; 5], t: f64) -> Point3 {
    let mut w = *cps;
    for level in (1..5).rev() {
        for i in 0..level {
            w[i] = Point3::new(
                w[i].x + (w[i + 1].x - w[i].x) * t,
                w[i].y + (w[i + 1].y - w[i].y) * t,
                w[i].z + (w[i + 1].z - w[i].z) * t,
            );
        }
    }
    w[0]
}

pub(crate) fn bernstein4_all(t: f64) -> [f64; 5] {
    const C: [f64; 5] = [1.0, 4.0, 6.0, 4.0, 1.0];
    let u = 1.0 - t;
    let mut out = [0.0_f64; 5];
    let mut tp = 1.0;
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = C[i] * tp * u.powi(4 - i as i32);
        tp *= t;
    }
    out
}

pub(crate) fn bernstein4_deriv_all(t: f64) -> [f64; 5] {
    // d/dt B_i^4 = 4 (B_{i-1}^3 − B_i^3).
    let u = 1.0 - t;
    let b3 = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
    [
        -4.0 * b3[0],
        4.0 * (b3[0] - b3[1]),
        4.0 * (b3[1] - b3[2]),
        4.0 * (b3[2] - b3[3]),
        4.0 * b3[3],
    ]
}

/// Gaussian elimination for the 5×5 Bernstein interpolation system.
// Column loops index two rows of the same matrix (pivot and target); the
// iterator form would need split borrows and obscure the elimination math.
#[allow(clippy::needless_range_loop)]
pub(crate) fn solve_5(mat: &[[f64; 5]; 5], rhs: &[[f64; 3]; 5]) -> Option<[[f64; 3]; 5]> {
    let mut a = *mat;
    let mut b = *rhs;
    for col in 0..5 {
        let mut piv = col;
        for r in (col + 1)..5 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for r in (col + 1)..5 {
            let f = a[r][col] / a[col][col];
            for c in col..5 {
                a[r][c] -= f * a[col][c];
            }
            for c in 0..3 {
                b[r][c] -= f * b[col][c];
            }
        }
    }
    let mut x = [[0.0_f64; 3]; 5];
    for row in (0..5).rev() {
        for c in 0..3 {
            let mut acc = b[row][c];
            for col in (row + 1)..5 {
                acc -= a[row][col] * x[col][c];
            }
            x[row][c] = acc / a[row][row];
        }
    }
    Some(x)
}
