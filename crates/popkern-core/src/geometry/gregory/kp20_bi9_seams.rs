//! Seam fits for the KP20 bi-9 resampling (PK-599): the band middles
//! (depth 0–2, along 3–6) of every internal lattice seam and every
//! sector-separating seam, re-derived from the exact composition while
//! the corner blocks stay at the feasible start.
//!
//! - Internal seams are parametric C² (the τ quadrant lattice): one
//!   jet fit — value `V`, first difference `F`, second difference `S`
//!   along the seam — writes both pieces' bands, so C² is exact by
//!   construction.
//! - Sector seams carry C⁰ and the paper's eq (8)/(9) with the eq (10)
//!   reparameterisation `b(u)`; the pair's band middles solve a small
//!   least-squares problem over the null space of those rows,
//!   collocated above their polynomial degree (eq (9) ≤ 11 at bi-9).

use faer::Mat;

use super::kp20_bi9_exact::ExactPiece;
use crate::geometry::Point3;
use crate::geometry::Vec3;
use crate::geometry::variational_cap::linalg;

/// A bi-9 control net, `[i][j]` with `i ↔ u`, `j ↔ v`.
pub(super) type Net = [[Point3; 10]; 10];

/// Which side of a piece a band belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Side {
    /// `u = 0` (depth along `i`, along `j`).
    I0,
    /// `u = 1`.
    I9,
    /// `v = 0` (depth along `j`, along `i`).
    J0,
    /// `v = 1`.
    J9,
}

/// Net index of the band coefficient at `depth` from `side`, `along`.
pub(super) fn slot(side: Side, depth: usize, along: usize) -> (usize, usize) {
    match side {
        Side::I0 => (depth, along),
        Side::I9 => (9 - depth, along),
        Side::J0 => (along, depth),
        Side::J9 => (along, 9 - depth),
    }
}

/// The free middle positions of a band row.
pub(super) const MIDDLE: [usize; 4] = [3, 4, 5, 6];

fn bernstein(n: usize, i: isize, t: f64) -> f64 {
    if i < 0 || i as usize > n {
        return 0.0;
    }
    let i = i as usize;
    let c = (0..i).fold(1.0, |c, k| c * (n - k) as f64 / (k + 1) as f64);
    c * t.powi(i as i32) * (1.0 - t).powi((n - i) as i32)
}

/// Degree-9 basis value, first and second derivative at `t`, index `j`.
pub(super) fn b9(j: usize, t: f64) -> [f64; 3] {
    let k = j as isize;
    [
        bernstein(9, k, t),
        9.0 * (bernstein(8, k - 1, t) - bernstein(8, k, t)),
        72.0 * (bernstein(7, k - 2, t) - 2.0 * bernstein(7, k - 1, t) + bernstein(7, k, t)),
    ]
}

fn v3(p: Point3) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}

fn p3(v: Vec3) -> Point3 {
    Point3::new(v.x, v.y, v.z)
}

/// Gauss–Legendre-like interior samples on [0, 1] (uniform midpoints;
/// enough of them to overdetermine every fit).
pub(super) fn samples(count: usize) -> Vec<f64> {
    (0..count)
        .map(|k| (k as f64 + 0.5) / count as f64)
        .collect()
}

/// Solve `min |A x − b|` per coordinate for 4 middle unknowns.
fn fit_middles(rows: &[[f64; 4]], rhs: &[Vec3]) -> Option<[Vec3; 4]> {
    let mut a = Mat::<f64>::zeros(rows.len(), 4);
    let mut b = Mat::<f64>::zeros(rows.len(), 3);
    for (r, (row, t)) in rows.iter().zip(rhs).enumerate() {
        for (c, w) in row.iter().enumerate() {
            a.write(r, c, *w);
        }
        b.write(r, 0, t.x);
        b.write(r, 1, t.y);
        b.write(r, 2, t.z);
    }
    let x = linalg::lstsq(a, b, 1e-12)?;
    Some(core::array::from_fn(|c| {
        Vec3::new(x.read(c, 0), x.read(c, 1), x.read(c, 2))
    }))
}

/// Fit the band middles of an internal (parametric C²) seam between
/// `minus` (seam on its `side_minus`, an `I9`/`J9` side) and `plus` (on
/// `side_plus`, the matching `I0`/`J0` side), and write both bands.
pub(super) fn fit_internal_seam(
    minus: &mut Net,
    plus: &mut Net,
    side_minus: Side,
    side_plus: Side,
    exact_minus: &ExactPiece<'_>,
    exact_plus: &ExactPiece<'_>,
) -> Option<()> {
    let row = |net: &Net, side: Side, depth: usize, along: usize| {
        let (i, j) = slot(side, depth, along);
        v3(net[i][j])
    };
    // Held jet coefficients from the plus side.
    let held = |along: usize| {
        let (p0, p1, p2) = (
            row(plus, side_plus, 0, along),
            row(plus, side_plus, 1, along),
            row(plus, side_plus, 2, along),
        );
        (p0, p1 - p0, p2 - p1 * 2.0 + p0)
    };
    let ts = samples(24);
    let mut targets = [Vec::new(), Vec::new(), Vec::new()];
    let mut rows = Vec::new();
    for &t in &ts {
        let (jp, jm) = match side_plus {
            Side::I0 => (exact_plus.jet(0.0, t), exact_minus.jet(1.0, t)),
            _ => (exact_plus.jet(t, 0.0), exact_minus.jet(t, 1.0)),
        };
        let (d1, d2) = match side_plus {
            Side::I0 => ((jp.u + jm.u) * 0.5, (jp.uu + jm.uu) * 0.5),
            _ => ((jp.v + jm.v) * 0.5, (jp.vv + jm.vv) * 0.5),
        };
        let goal = [(jp.p + jm.p) * 0.5, d1 * (1.0 / 9.0), d2 * (1.0 / 72.0)];
        let mut known = [Vec3::new(0.0, 0.0, 0.0); 3];
        let mut basis = [0.0; 4];
        for along in 0..10 {
            let w = b9(along, t)[0];
            match MIDDLE.iter().position(|&m| m == along) {
                Some(k) => basis[k] = w,
                None => {
                    let (v, f, s) = held(along);
                    known[0] = known[0] + v * w;
                    known[1] = known[1] + f * w;
                    known[2] = known[2] + s * w;
                }
            }
        }
        rows.push(basis);
        for q in 0..3 {
            targets[q].push(goal[q] - known[q]);
        }
    }
    let v = fit_middles(&rows, &targets[0])?;
    let f = fit_middles(&rows, &targets[1])?;
    let s = fit_middles(&rows, &targets[2])?;
    for (k, &along) in MIDDLE.iter().enumerate() {
        let writes = [
            (side_plus, 0, v[k]),
            (side_plus, 1, v[k] + f[k]),
            (side_plus, 2, v[k] + f[k] * 2.0 + s[k]),
            (side_minus, 0, v[k]),
            (side_minus, 1, v[k] - f[k]),
            (side_minus, 2, v[k] - f[k] * 2.0 + s[k]),
        ];
        for (side, depth, value) in writes {
            let (i, j) = slot(side, depth, along);
            let net = if side == side_plus {
                &mut *plus
            } else {
                &mut *minus
            };
            net[i][j] = p3(value);
        }
    }
    Some(())
}

/// Fit one side's band middles (rows 0–2) to the exact composition's
/// cross-derivatives along that side alone.
pub(super) fn fit_band_side(net: &mut Net, side: Side, exact: &ExactPiece<'_>) -> Option<()> {
    let row = |net: &Net, depth: usize, along: usize| {
        let (i, j) = slot(side, depth, along);
        v3(net[i][j])
    };
    let held = |along: usize| {
        let (p0, p1, p2) = (row(net, 0, along), row(net, 1, along), row(net, 2, along));
        (p0, p1 - p0, p2 - p1 * 2.0 + p0)
    };
    let ts = samples(24);
    let mut targets = [Vec::new(), Vec::new(), Vec::new()];
    let mut rows = Vec::new();
    for &t in &ts {
        // Inward first derivative (sign flips on the far sides) and
        // second derivative along the depth direction.
        let (jet, d1, d2) = match side {
            Side::I0 => {
                let j = exact.jet(0.0, t);
                (j, j.u, j.uu)
            }
            Side::I9 => {
                let j = exact.jet(1.0, t);
                (j, j.u * -1.0, j.uu)
            }
            Side::J0 => {
                let j = exact.jet(t, 0.0);
                (j, j.v, j.vv)
            }
            Side::J9 => {
                let j = exact.jet(t, 1.0);
                (j, j.v * -1.0, j.vv)
            }
        };
        let goal = [jet.p, d1 * (1.0 / 9.0), d2 * (1.0 / 72.0)];
        let mut known = [Vec3::new(0.0, 0.0, 0.0); 3];
        let mut basis = [0.0; 4];
        for along in 0..10 {
            let w = b9(along, t)[0];
            match MIDDLE.iter().position(|&m| m == along) {
                Some(k) => basis[k] = w,
                None => {
                    let (v, f, s) = held(along);
                    known[0] = known[0] + v * w;
                    known[1] = known[1] + f * w;
                    known[2] = known[2] + s * w;
                }
            }
        }
        rows.push(basis);
        for q in 0..3 {
            targets[q].push(goal[q] - known[q]);
        }
    }
    let v = fit_middles(&rows, &targets[0])?;
    let f = fit_middles(&rows, &targets[1])?;
    let s = fit_middles(&rows, &targets[2])?;
    for (k, &along) in MIDDLE.iter().enumerate() {
        for (depth, value) in [(0, v[k]), (1, v[k] + f[k]), (2, v[k] + f[k] * 2.0 + s[k])] {
            let (i, j) = slot(side, depth, along);
            net[i][j] = p3(value);
        }
    }
    Some(())
}

#[path = "kp20_bi9_sector_seam.rs"]
mod kp20_bi9_sector_seam;
pub(super) use kp20_bi9_sector_seam::{fit_sector_seam, sector_seam_residual};
