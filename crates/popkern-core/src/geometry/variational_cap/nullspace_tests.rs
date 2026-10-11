//! The structured null space and solves against their dense
//! definitions (PK-620).

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use faer::Mat;

use super::chart::CharacteristicMap;
use super::nullspace::hard_system;
use super::space::SectorSpace;
use super::spoke::spoke_rows;

/// Dense `Z` (free rows × k) from the basis rows.
fn dense_z(hs: &super::nullspace::HardSystem) -> Mat<f64> {
    let k = hs.dim();
    let mut z = Mat::<f64>::zeros(hs.free.len(), k);
    for r in 0..hs.free.len() {
        let mut row = vec![0.0; k];
        hs.add_row(r, 1.0, &mut row);
        for (j, v) in row.into_iter().enumerate() {
            z.write(r, j, v);
        }
    }
    z
}

fn spaces() -> Vec<SectorSpace> {
    vec![
        SectorSpace::bezier(3, 7),
        SectorSpace::bezier(5, 7),
        SectorSpace::new(
            4,
            5,
            vec![
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            ],
        ),
        SectorSpace::bezier(8, 9).with_corner_cells(),
        SectorSpace::bezier(16, 9).with_corner_cells(),
    ]
}

#[test]
fn the_nested_basis_is_an_orthonormal_null_space_basis() {
    for space in spaces() {
        let hs = hard_system(&space).expect("hard system");
        let z = dense_z(&hs);
        let k = hs.dim();
        // Orthonormal.
        let mut worst_orth = 0.0_f64;
        for a in 0..k {
            for b in a..k {
                let dot: f64 = (0..z.nrows()).map(|r| z.read(r, a) * z.read(r, b)).sum();
                let target = if a == b { 1.0 } else { 0.0 };
                worst_orth = worst_orth.max((dot - target).abs());
            }
        }
        // In the spoke rows' null space.
        let rows = spoke_rows(&space, CharacteristicMap::new(space.n).c);
        let mut worst_null = 0.0_f64;
        for row in &rows {
            let scale = row.iter().map(|(_, w)| w.abs()).fold(0.0, f64::max);
            for j in 0..k {
                let mut acc = 0.0;
                for &(i, w) in row {
                    if let (false, r) = hs.slots[i] {
                        acc += w * z.read(r, j);
                    }
                }
                worst_null = worst_null.max(acc.abs() / scale.max(1.0));
            }
        }
        eprintln!(
            "n {} degree {} cells {}: k {} orth {worst_orth:.2e} null {worst_null:.2e}",
            space.n, space.degree, space.corner_cells, k
        );
        assert!(worst_orth < 1e-10, "orthonormality {worst_orth:.2e}");
        assert!(worst_null < 1e-10, "null space {worst_null:.2e}");
    }
}

/// A level of the 5-sided spherical cap, its characteristic energy,
/// and the dense reference reduction `Zᵀ E Z / mean + 1e-12 I` with
/// `h = Zᵀ (E x_c − b) / mean`.
fn reference_level(
    space: SectorSpace,
) -> (
    super::level::Level,
    super::factor::EnergyFactor,
    Mat<f64>,
    Vec<f64>,
) {
    use super::energy::{Integrand, energy_gram};
    use super::level::{Frame, build_level};
    let sides = super::tests::spherical_ngon(space.n, 0.6);
    let frame = Frame {
        origin: crate::geometry::Point3::new(0.0, 0.0, 0.0),
        scale: 1.0,
    };
    let splits = vec![0.5; space.n];
    let level = build_level(space, &sides, &frame, &splits).expect("level");
    let phi = CharacteristicMap::new(level.space.n);
    let metric = |_s: usize, u: f64, v: f64| phi.metric(u, v);
    let integrand = Integrand::Parametric {
        metric: &metric,
        eps: 1e-2,
    };
    let factor = super::factor::energy_factor(&level, &integrand, None).expect("factor");
    let gram = energy_gram(&level.space, &integrand, None);
    let hs = &level.hs;
    let z = dense_z(hs);
    let (k, per, size) = (hs.dim(), level.space.per(), gram.size);
    let mut g = Mat::<f64>::zeros(3 * k, 3 * k);
    let mut h = vec![0.0; 3 * k];
    for s in 0..level.space.n {
        let block = &gram.blocks[s];
        let base = s * per;
        for c in 0..3 {
            for d in 0..3 {
                for la in 0..per {
                    let (false, ra) = hs.slots[base + la] else {
                        continue;
                    };
                    for lb in 0..per {
                        let e = block[(c * per + la) * size + d * per + lb];
                        if e == 0.0 {
                            continue;
                        }
                        match hs.slots[base + lb] {
                            (false, rb) => {
                                for i in 0..k {
                                    let zi = z.read(ra, i) * e;
                                    if zi == 0.0 {
                                        continue;
                                    }
                                    for j in 0..k {
                                        g.write(
                                            c * k + i,
                                            d * k + j,
                                            g.read(c * k + i, d * k + j) + zi * z.read(rb, j),
                                        );
                                    }
                                    h[c * k + i] += zi * level.xp[d][rb];
                                }
                            }
                            (true, rb) => {
                                for i in 0..k {
                                    h[c * k + i] += z.read(ra, i) * e * level.fixed[d][rb];
                                }
                            }
                        }
                    }
                }
            }
            for la in 0..per {
                if let (false, ra) = hs.slots[base + la] {
                    for i in 0..k {
                        h[c * k + i] -= z.read(ra, i) * gram.linear[s][c * per + la];
                    }
                }
            }
        }
    }
    let mean = (0..3 * k).map(|i| g.read(i, i)).sum::<f64>() / (3 * k) as f64;
    for i in 0..3 * k {
        for j in 0..3 * k {
            g.write(i, j, g.read(i, j) / mean);
        }
        g.write(i, i, g.read(i, i) + 1e-12);
    }
    for v in &mut h {
        *v /= mean;
    }
    (level, factor, g, h)
}

#[test]
fn the_sector_energy_rows_reproduce_the_dense_reduction() {
    for space in [
        SectorSpace::bezier(5, 7),
        SectorSpace::bezier(8, 9).with_corner_cells(),
    ] {
        let (_, factor, g, h) = reference_level(space);
        let dense = factor.dense_gram();
        let scale = (0..g.nrows()).map(|i| g.read(i, i)).fold(0.0, f64::max);
        let mut worst_g = 0.0_f64;
        for i in 0..g.nrows() {
            for j in 0..g.ncols() {
                worst_g = worst_g.max((dense.read(i, j) - g.read(i, j)).abs());
            }
        }
        let hscale = h.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let worst_h = h
            .iter()
            .zip(&factor.h)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        eprintln!(
            "gram {:.2e} (scale {scale:.2e}) h {:.2e} (scale {hscale:.2e})",
            worst_g, worst_h
        );
        assert!(worst_g < 1e-9 * scale, "gram {worst_g:.2e}");
        assert!(worst_h < 1e-9 * hscale.max(1.0), "h {worst_h:.2e}");
    }
}

#[test]
fn the_banded_solve_is_the_dense_stacked_solve() {
    for space in [
        SectorSpace::bezier(5, 7),
        SectorSpace::bezier(8, 9).with_corner_cells(),
    ] {
        let (_, factor, g, h) = reference_level(space);
        let hs = factor.hs.clone();
        let k = hs.dim();
        // The minimiser solves G y = −h.
        let y = factor.minimiser().expect("minimiser");
        let mut worst = 0.0_f64;
        for (i, &hi) in h.iter().enumerate() {
            let r: f64 = (0..3 * k)
                .map(|j| g.read(i, j) * y[j / k][j % k])
                .sum::<f64>()
                + hi;
            worst = worst.max(r.abs());
        }
        eprintln!("minimiser residual {worst:.2e}");
        assert!(worst < 1e-8, "minimiser residual {worst:.2e}");
        // Rim-like rows: each over one sector's columns, weighted 1e3.
        let n = hs.layout.interior.len();
        let per_sector = 6;
        let mut rows = Mat::<f64>::zeros(n * per_sector, 3 * k);
        let mut rhs = vec![0.0; n * per_sector];
        for s in 0..n {
            for t in 0..per_sector {
                let i = s * per_sector + t;
                for (a, &j) in hs.sector_cols[s].iter().enumerate() {
                    let c = (t + a) % 3;
                    rows.write(i, c * k + j, 1e3 * ((i * 31 + a * 7) as f64 * 0.13).sin());
                }
                rhs[i] = 1e3 * (i as f64 * 0.7).cos();
            }
        }
        let lambda = 1e-2;
        let banded = super::banded::solve(&factor, lambda, &rows, &rhs).expect("banded");
        let (er, ec) = factor.dense_rows();
        let total = er.nrows() + rows.nrows();
        let sl = lambda.sqrt();
        let a = Mat::from_fn(total, 3 * k, |i, j| {
            if i < er.nrows() {
                sl * er.read(i, j)
            } else {
                rows.read(i - er.nrows(), j)
            }
        });
        let b = Mat::from_fn(total, 1, |i, _| {
            if i < er.nrows() {
                -sl * ec[i]
            } else {
                rhs[i - er.nrows()]
            }
        });
        let dense = super::linalg::lstsq(a, b, 1e-14).expect("dense");
        let scale = (0..3 * k)
            .map(|j| dense.read(j, 0).abs())
            .fold(0.0, f64::max);
        let worst = (0..3 * k)
            .map(|j| (dense.read(j, 0) - banded[j / k][j % k]).abs())
            .fold(0.0, f64::max);
        eprintln!("banded vs dense {worst:.2e} (scale {scale:.2e})");
        assert!(worst < 1e-8 * scale.max(1.0), "banded vs dense {worst:.2e}");
    }
}
