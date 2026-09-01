//! Degree bookkeeping for the KP20 paper border rows (PK-398).
//!
//! The paper border carries its rim rows at native degree 6 (seven
//! control points along each piece) so the setback fan's exact strip
//! section reaches the band machinery without a cubic detour. This
//! module owns the row-level degree conversions the border pipeline
//! needs:
//!
//! - exact elevation of bi-3 rows (the characteristic border, the
//!   production tensor-border) to degree 6;
//! - end-pinned least-squares reduction of degree-6 rows to cubics
//!   (the KP15 TBdata backstop consumes bi-3 rows; exact when the
//!   content is an elevated cubic);
//! - the eq (14) product `d · (row1 − row0)` with a cubic scalar
//!   profile and a degree-6 vector row — degree 9 — reduced to degree
//!   6 while preserving BOTH end 2-jets exactly. Every spoke/corner
//!   relation of the assembly consumes only those end jets (the three
//!   end coefficients of the band rows), so the reduction changes
//!   nothing the strip projection or corner overlap can see; the one
//!   free interior coefficient is fitted in least squares. Exact when
//!   the product has true degree ≤ 6 (elevated-cubic inputs).
//! - de Casteljau halving and reversal of degree-6 rows.

use crate::geometry::{Point3, Vec3};

/// Seven-coefficient row: one degree-6 Bernstein row of a border piece.
pub type SexticRow = [Point3; 7];

/// Elevate one cubic Bernstein row to degree 6 exactly.
pub fn elevate_cubic_row_to_sextic(row: [Point3; 4]) -> SexticRow {
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let binom3 = |i: usize| -> f64 {
        match i {
            0 | 3 => 1.0,
            1 | 2 => 3.0,
            _ => 0.0,
        }
    };
    let mut out = [Point3::new(0.0, 0.0, 0.0); 7];
    for (i, slot) in out.iter_mut().enumerate() {
        let j_min = i.saturating_sub(3);
        let j_max = i.min(3);
        let mut acc = [0.0_f64; 3];
        for (j, cp) in row.iter().enumerate().take(j_max + 1).skip(j_min) {
            let w = C3[j] * binom3(i - j) / C6[i];
            acc[0] += w * cp.x;
            acc[1] += w * cp.y;
            acc[2] += w * cp.z;
        }
        *slot = Point3::new(acc[0], acc[1], acc[2]);
    }
    out
}

/// Elevate a bi-3 border piece (three depth rows of four) to degree 6
/// along, exactly.
pub fn elevate_bi3_rows_to_sextic(rows: [[Point3; 4]; 3]) -> [SexticRow; 3] {
    [
        elevate_cubic_row_to_sextic(rows[0]),
        elevate_cubic_row_to_sextic(rows[1]),
        elevate_cubic_row_to_sextic(rows[2]),
    ]
}

/// Degree-6 Bernstein basis at `t`.
pub fn bernstein6(t: f64) -> [f64; 7] {
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let u = 1.0 - t;
    let mut out = [0.0_f64; 7];
    let mut tp = 1.0;
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = C6[i] * tp * u.powi(6 - i as i32);
        tp *= t;
    }
    out
}

/// Evaluate a degree-6 row at `t`.
pub fn eval_sextic_row(row: &SexticRow, t: f64) -> Point3 {
    let b = bernstein6(t);
    let mut acc = [0.0_f64; 3];
    for (w, cp) in b.iter().zip(row) {
        acc[0] += w * cp.x;
        acc[1] += w * cp.y;
        acc[2] += w * cp.z;
    }
    Point3::new(acc[0], acc[1], acc[2])
}

/// Reduce one degree-6 row to a cubic with both endpoints pinned and
/// the two interior coefficients fitted in least squares over a
/// uniform sample grid. Exact (to roundoff) for elevated cubics.
pub fn reduce_sextic_row_to_cubic(row: &SexticRow) -> [Point3; 4] {
    const M: usize = 33;
    let mut gram = [[0.0_f64; 2]; 2];
    let mut moment = [[0.0_f64; 3]; 2];
    for k in 0..M {
        let t = k as f64 / (M - 1) as f64;
        let u = 1.0 - t;
        let b = [u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t];
        let target = eval_sextic_row(row, t);
        let r = [
            target.x - b[0] * row[0].x - b[3] * row[6].x,
            target.y - b[0] * row[0].y - b[3] * row[6].y,
            target.z - b[0] * row[0].z - b[3] * row[6].z,
        ];
        for i in 0..2 {
            for j in 0..2 {
                gram[i][j] += b[i + 1] * b[j + 1];
            }
            for c in 0..3 {
                moment[i][c] += b[i + 1] * r[c];
            }
        }
    }
    let det = gram[0][0] * gram[1][1] - gram[0][1] * gram[1][0];
    let mut p1 = [0.0_f64; 3];
    let mut p2 = [0.0_f64; 3];
    for c in 0..3 {
        p1[c] = (moment[0][c] * gram[1][1] - moment[1][c] * gram[0][1]) / det;
        p2[c] = (gram[0][0] * moment[1][c] - gram[1][0] * moment[0][c]) / det;
    }
    [
        row[0],
        Point3::new(p1[0], p1[1], p1[2]),
        Point3::new(p2[0], p2[1], p2[2]),
        row[6],
    ]
}

/// Reduce a degree-6 border piece to bi-3 rows (see
/// [`reduce_sextic_row_to_cubic`]).
pub fn reduce_sextic_rows_to_bi3(rows: &[SexticRow; 3]) -> [[Point3; 4]; 3] {
    [
        reduce_sextic_row_to_cubic(&rows[0]),
        reduce_sextic_row_to_cubic(&rows[1]),
        reduce_sextic_row_to_cubic(&rows[2]),
    ]
}

/// De Casteljau halves of one degree-6 row at `t = ½`.
pub fn split_sextic_row(row: &SexticRow) -> (SexticRow, SexticRow) {
    let mid =
        |a: Point3, b: Point3| Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
    let mut levels: [[Point3; 7]; 7] = [[Point3::new(0.0, 0.0, 0.0); 7]; 7];
    levels[0] = *row;
    for level in 1..7 {
        for i in 0..(7 - level) {
            levels[level][i] = mid(levels[level - 1][i], levels[level - 1][i + 1]);
        }
    }
    let first: SexticRow = core::array::from_fn(|i| levels[i][0]);
    let second: SexticRow = core::array::from_fn(|i| levels[6 - i][i]);
    (first, second)
}

/// De Casteljau halves of a three-row piece.
pub fn split_sextic_rows(rows: &[SexticRow; 3]) -> ([SexticRow; 3], [SexticRow; 3]) {
    let (a0, b0) = split_sextic_row(&rows[0]);
    let (a1, b1) = split_sextic_row(&rows[1]);
    let (a2, b2) = split_sextic_row(&rows[2]);
    ([a0, a1, a2], [b0, b1, b2])
}

/// Reverse the along direction of a three-row piece.
pub fn reverse_sextic_rows(rows: &[SexticRow; 3]) -> [SexticRow; 3] {
    let mut out = *rows;
    for row in out.iter_mut() {
        row.reverse();
    }
    out
}

/// Product of a cubic scalar profile `d` (Bernstein coefficients) and
/// a degree-6 vector row, reduced from degree 9 to degree 6 with both
/// end 2-jets preserved exactly and the single remaining coefficient
/// fitted in least squares. Exact when the product's true degree is
/// at most 6.
pub fn scalar_cubic_times_sextic_reduced(d: [f64; 4], v: [Vec3; 7]) -> [Vec3; 7] {
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    const C9: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    // Degree-9 product coefficients.
    let mut q = [Vec3::zero(); 10];
    for (i, slot) in q.iter_mut().enumerate() {
        let j_min = i.saturating_sub(6);
        let j_max = i.min(3);
        for j in j_min..=j_max {
            let l = i - j;
            *slot = *slot + v[l] * (d[j] * C3[j] * C6[l] / C9[i]);
        }
    }
    reduce_nonic_to_sextic_end_jets(&q)
}

/// Degree-9 → degree-6 reduction preserving the end 2-jets exactly:
/// `p₀ = q₀`, `6(p₁ − p₀) = 9(q₁ − q₀)`,
/// `30(p₂ − 2p₁ + p₀) = 72(q₂ − 2q₁ + q₀)` and mirrored at the far
/// end; `p₃` least squares over a uniform sample grid.
pub fn reduce_nonic_to_sextic_end_jets(q: &[Vec3; 10]) -> [Vec3; 7] {
    let p0 = q[0];
    let p1 = p0 + (q[1] - q[0]) * 1.5;
    let p2 = p1 * 2.0 - p0 + (q[2] - q[1] * 2.0 + q[0]) * 2.4;
    let p6 = q[9];
    let p5 = p6 + (q[8] - q[9]) * 1.5;
    let p4 = p5 * 2.0 - p6 + (q[7] - q[8] * 2.0 + q[9]) * 2.4;
    const M: usize = 33;
    let mut gram = 0.0_f64;
    let mut moment = Vec3::zero();
    for k in 0..M {
        let t = k as f64 / (M - 1) as f64;
        let b6 = bernstein6(t);
        let b9 = bernstein9(t);
        let mut target = Vec3::zero();
        for (w, c) in b9.iter().zip(q) {
            target = target + *c * *w;
        }
        let partial = p0 * b6[0] + p1 * b6[1] + p2 * b6[2] + p4 * b6[4] + p5 * b6[5] + p6 * b6[6];
        gram += b6[3] * b6[3];
        moment = moment + (target - partial) * b6[3];
    }
    let p3 = if gram > 0.0 {
        moment * (1.0 / gram)
    } else {
        Vec3::zero()
    };
    [p0, p1, p2, p3, p4, p5, p6]
}

fn bernstein9(t: f64) -> [f64; 10] {
    const C9: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    let u = 1.0 - t;
    let mut out = [0.0_f64; 10];
    let mut tp = 1.0;
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = C9[i] * tp * u.powi(9 - i as i32);
        tp *= t;
    }
    out
}

#[cfg(test)]
#[path = "kp20_border_degree_tests.rs"]
mod tests;
