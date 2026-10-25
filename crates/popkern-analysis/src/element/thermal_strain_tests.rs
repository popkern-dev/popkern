use super::*;
use crate::element::VectorTables;
use crate::element::vector::lame;
use crate::grid::CELL_CORNERS;
use crate::moments::moment_index;

/// Moments (physical volume) of the box [lo, hi] (local coordinates)
/// in a cell of size h.
fn box_moments(lo: [f64; 3], hi: [f64; 3], h: f64) -> [f64; 27] {
    let one = |a: usize, k: usize| {
        let p = k as i32 + 1;
        (hi[a].powi(p) - lo[a].powi(p)) / p as f64
    };
    let jac = (0.5 * h).powi(3);
    let mut m = [0.0; 27];
    for c in 0..3 {
        for b in 0..3 {
            for a in 0..3 {
                m[moment_index(a, b, c)] = jac * one(0, a) * one(1, b) * one(2, c);
            }
        }
    }
    m
}

/// A uniform rise's free expansion is stress-free, so the cell's
/// stiffness times it is exactly its thermal load (full and cut cells,
/// with and without enhanced modes).
#[test]
fn free_expansion_balances_the_thermal_load() {
    let h = 0.5;
    let (lambda, mu) = lame(1000.0, 0.3);
    let (alpha, theta) = (2.0e-5, 40.0);
    let vt = VectorTables::new();
    for m in [
        box_moments([-1.0; 3], [1.0; 3], h),
        box_moments([-1.0, -0.4, -1.0], [0.3, 1.0, 0.7], h),
    ] {
        for (op, _) in [
            vt.cell_operator(&m, h, lambda, mu),
            Ok((
                vt.cell(&m, h, lambda, mu, crate::element::vector::Enhancement::None)
                    .unwrap_or_else(|e| panic!("{e}")),
                false,
            )),
        ]
        .into_iter()
        .map(|r| r.unwrap_or_else(|e| panic!("{e}")))
        {
            let t = cell_thermal(&op, &m, h, (lambda, mu), alpha, &[theta; 8]);
            // u = αθ x from the cell centre, at the corners.
            let u: [f64; NODAL_DOFS] = std::array::from_fn(|r| {
                alpha * theta * 0.5 * h * (2.0 * CELL_CORNERS[r / 3][r % 3] as f64 - 1.0)
            });
            let scale = t.load.iter().fold(0.0f64, |a, x| a.max(x.abs()));
            for r in 0..NODAL_DOFS {
                let ku: f64 = (0..NODAL_DOFS).map(|c| op.k[r][c] * u[c]).sum();
                assert!(
                    (ku - t.load[r]).abs() <= 1e-10 * scale,
                    "{r}: {ku} vs {}",
                    t.load[r]
                );
            }
            assert!(t.alpha0.iter().all(|a| a.abs() <= 1e-12 * alpha * theta));
        }
    }
}
