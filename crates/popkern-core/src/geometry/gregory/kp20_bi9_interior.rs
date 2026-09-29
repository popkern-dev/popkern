//! Interior fit of one KP20 bi-9 piece (PK-599): the 4×4 coefficients
//! away from every band, least squares to the exact composition with
//! all band and corner coefficients held.

// Coefficient nets are addressed by (i, j) index pairs throughout.
#![allow(clippy::needless_range_loop)]

use faer::Mat;

use super::kp20_bi9_exact::ExactPiece;
use super::kp20_bi9_seams::{Net, b9, samples};
use crate::geometry::Point3;
use crate::geometry::variational_cap::linalg;

/// Re-fit `net[3..=6][3..=6]` to `exact`.
pub(super) fn fit_interior(net: &mut Net, exact: &ExactPiece<'_>) -> Option<()> {
    let free = |i: usize, j: usize| (3..=6).contains(&i) && (3..=6).contains(&j);
    let ts = samples(14);
    let rows = ts.len() * ts.len();
    let mut a = Mat::<f64>::zeros(rows, 16);
    let mut b = Mat::<f64>::zeros(rows, 3);
    let mut r = 0;
    for &u in &ts {
        for &v in &ts {
            let bu: Vec<f64> = (0..10).map(|i| b9(i, u)[0]).collect();
            let bv: Vec<f64> = (0..10).map(|j| b9(j, v)[0]).collect();
            let target = exact.point(u, v);
            let (mut kx, mut ky, mut kz) = (0.0, 0.0, 0.0);
            for i in 0..10 {
                for j in 0..10 {
                    let w = bu[i] * bv[j];
                    if free(i, j) {
                        a.write(r, (i - 3) * 4 + (j - 3), w);
                    } else {
                        kx += w * net[i][j].x;
                        ky += w * net[i][j].y;
                        kz += w * net[i][j].z;
                    }
                }
            }
            b.write(r, 0, target.x - kx);
            b.write(r, 1, target.y - ky);
            b.write(r, 2, target.z - kz);
            r += 1;
        }
    }
    let x = linalg::lstsq(a, b, 1e-12)?;
    for i in 3..=6 {
        for j in 3..=6 {
            let k = (i - 3) * 4 + (j - 3);
            net[i][j] = Point3::new(x.read(k, 0), x.read(k, 1), x.read(k, 2));
        }
    }
    Some(())
}
