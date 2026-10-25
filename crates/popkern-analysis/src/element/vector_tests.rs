use super::vector::*;
use crate::grid::CELL_CORNERS;
use crate::linalg::Cholesky;

fn full_moments(h: f64) -> [f64; 27] {
    let axis = [2.0, 0.0, 2.0 / 3.0];
    let s = 0.125 * h * h * h;
    std::array::from_fn(|i| s * axis[i % 3] * axis[(i / 3) % 3] * axis[i / 9])
}

/// Corner positions of a cell of size h centred at the origin.
fn corners(h: f64) -> [[f64; 3]; 8] {
    CELL_CORNERS.map(|c| c.map(|o| if o == 0 { -0.5 * h } else { 0.5 * h }))
}

fn apply(k: &Matrix24, u: &[f64; NODAL_DOFS]) -> [f64; NODAL_DOFS] {
    std::array::from_fn(|r| (0..NODAL_DOFS).map(|c| k[r][c] * u[c]).sum())
}

fn cell(h: f64, e: f64, nu: f64, enhanced: bool) -> VectorCell {
    let (l, m) = lame(e, nu);
    VectorTables::new()
        .cell(
            &full_moments(h),
            h,
            l,
            m,
            if enhanced {
                Enhancement::FragmentMean
            } else {
                Enhancement::None
            },
        )
        .unwrap_or_else(|e| panic!("{e}"))
}

/// Nodal values of the field `f` at the corners.
fn nodal(h: f64, f: impl Fn([f64; 3]) -> [f64; 3]) -> [f64; NODAL_DOFS] {
    let x = corners(h);
    std::array::from_fn(|d| f(x[d / 3])[d % 3])
}

#[test]
fn x1_the_condensed_cell_has_exactly_the_six_rigid_modes() {
    let h = 0.7;
    let c = cell(h, 3.0, 0.3, true);
    let scale = c.k.iter().flatten().fold(0.0f64, |a, x| a.max(x.abs()));
    let mut modes: Vec<[f64; NODAL_DOFS]> = Vec::new();
    for t in 0..3 {
        modes.push(nodal(h, |_| {
            std::array::from_fn(|i| if i == t { 1.0 } else { 0.0 })
        }));
    }
    for axis in 0..3 {
        modes.push(nodal(h, |x| {
            let mut w = [0.0; 3];
            w[axis] = 1.0;
            [
                w[1] * x[2] - w[2] * x[1],
                w[2] * x[0] - w[0] * x[2],
                w[0] * x[1] - w[1] * x[0],
            ]
        }));
    }
    for (i, r) in modes.iter().enumerate() {
        let kr = apply(&c.k, r);
        let worst = kr.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        assert!(worst <= 1e-13 * scale, "mode {i}: {worst}");
    }
    for r in 0..NODAL_DOFS {
        for col in 0..NODAL_DOFS {
            assert!((c.k[r][col] - c.k[col][r]).abs() <= 1e-14 * scale);
        }
    }
    // Rank 18: adding the rigid modes makes it positive definite.
    let mut flat = vec![0.0; NODAL_DOFS * NODAL_DOFS];
    for r in 0..NODAL_DOFS {
        for col in 0..NODAL_DOFS {
            let rr: f64 = modes.iter().map(|m| m[r] * m[col]).sum();
            flat[r * NODAL_DOFS + col] = c.k[r][col] + scale * rr;
        }
    }
    assert!(Cholesky::factor(&flat, NODAL_DOFS, 1e-10).is_ok());
}

/// Boundary load ∮ N_a σ·n dA of a constant stress on a full cell.
fn traction_load(h: f64, sigma: [[f64; 3]; 3]) -> [f64; NODAL_DOFS] {
    let x = corners(h);
    let mut f = [0.0; NODAL_DOFS];
    for (a, xa) in x.iter().enumerate() {
        for axis in 0..3 {
            // each node lies on one wall per axis; ∫ N_a over it is h²/4
            let n = if xa[axis] > 0.0 { 1.0 } else { -1.0 };
            for i in 0..3 {
                f[3 * a + i] += sigma[i][axis] * n * 0.25 * h * h;
            }
        }
    }
    f
}

#[test]
fn x2_linear_fields_pass_the_patch_test_with_zero_enhancement() {
    let (h, e, nu) = (0.6, 2.0, 0.3);
    let (l, m) = lame(e, nu);
    let c = cell(h, e, nu, true);
    for grad in [
        [[1e-3, 0.0, 0.0], [0.0; 3], [0.0; 3]],
        [[0.0, 2e-3, 0.0], [0.0; 3], [0.0; 3]],
    ] {
        let u = nodal(h, |x| {
            std::array::from_fn(|i| (0..3).map(|j| grad[i][j] * x[j]).sum())
        });
        let eps: [[f64; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (grad[i][j] + grad[j][i])));
        let tr = eps[0][0] + eps[1][1] + eps[2][2];
        let sigma = std::array::from_fn(|i| {
            std::array::from_fn(|j| 2.0 * m * eps[i][j] + if i == j { l * tr } else { 0.0 })
        });
        let ku = apply(&c.k, &u);
        let f = traction_load(h, sigma);
        let scale = f.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        for d in 0..NODAL_DOFS {
            assert!((ku[d] - f[d]).abs() <= 1e-13 * scale, "dof {d}");
        }
        let recovery = c.recovery.unwrap_or_else(|| panic!("no recovery"));
        for row in recovery {
            let alpha: f64 = row.iter().zip(&u).map(|(r, x)| r * x).sum();
            assert!(
                alpha.abs() <= 1e-15,
                "enhancement {alpha} for a linear field"
            );
        }
    }
}

#[test]
fn c3_a_single_cell_bends_exactly() {
    let (h, e, nu, kappa) = (0.8, 5.0, 0.3, 1e-2);
    // Pure bending about y: σ_xx = E κ z, other stresses zero.
    let field = |x: [f64; 3]| {
        [
            kappa * x[0] * x[2],
            -nu * kappa * x[1] * x[2],
            -0.5 * kappa * (x[0] * x[0] + nu * (x[2] * x[2] - x[1] * x[1])),
        ]
    };
    let u = nodal(h, field);
    let exact = e * kappa * kappa * h.powi(5) / 24.0;
    let energy = |k: &Matrix24| 0.5 * u.iter().zip(apply(k, &u)).map(|(a, b)| a * b).sum::<f64>();
    let enhanced = energy(&cell(h, e, nu, true).k);
    assert!(
        (enhanced - exact).abs() <= 1e-12 * exact,
        "{enhanced} vs {exact}"
    );
    // Negative control: plain Q1 is too stiff in bending.
    let plain = energy(&cell(h, e, nu, false).k);
    assert!(plain > 1.5 * exact, "plain Q1 {plain} vs {exact}");
}

/// Tip deflection of a cantilever L × b × t (b = t = 1, along x), E = 1,
/// clamped at x = 0, end shear P in −z, on full cells of size h.
fn cantilever_tip(nu: f64, h: f64, enhanced: bool) -> f64 {
    let (l, p) = (10.0, 1e-3);
    let n = [
        (l / h).round() as usize,
        (1.0 / h).round() as usize,
        (1.0 / h).round() as usize,
    ];
    let nodes = [n[0] + 1, n[1] + 1, n[2] + 1];
    let node = |i: usize, j: usize, k: usize| i + nodes[0] * (j + nodes[1] * k);
    let dofs = 3 * nodes[0] * nodes[1] * nodes[2];
    let (lam, mu) = lame(1.0, nu);
    let c = VectorTables::new()
        .cell(
            &full_moments(h),
            h,
            lam,
            mu,
            if enhanced {
                Enhancement::FragmentMean
            } else {
                Enhancement::None
            },
        )
        .unwrap_or_else(|e| panic!("{e}"));
    let mut k = vec![0.0; dofs * dofs];
    let mut f = vec![0.0; dofs];
    for ck in 0..n[2] {
        for cj in 0..n[1] {
            for ci in 0..n[0] {
                let ids: [usize; 8] = CELL_CORNERS.map(|o| node(ci + o[0], cj + o[1], ck + o[2]));
                for a in 0..8 {
                    for b in 0..8 {
                        for i in 0..3 {
                            for j in 0..3 {
                                k[(3 * ids[a] + i) * dofs + 3 * ids[b] + j] +=
                                    c.k[3 * a + i][3 * b + j];
                            }
                        }
                    }
                }
            }
        }
    }
    // Uniform shear traction on the end face: h²/4 per corner per cell.
    let traction = p / 1.0;
    for ck in 0..n[2] {
        for cj in 0..n[1] {
            for (dj, dk) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                f[3 * node(n[0], cj + dj, ck + dk) + 2] -= traction * 0.25 * h * h;
            }
        }
    }
    // Beam-theory root: u_x = 0 over the x = 0 node plane, u_z = 0 on
    // its mid-height line and u_y = 0 on its mid-width line. This removes
    // the rigid modes but leaves the root section free to contract
    // laterally, as beam theory assumes (a fully clamped root restrains
    // Poisson contraction and is ~2.6% stiffer at ν = 0.49).
    let (mid_j, mid_k) = (n[1] / 2, n[2] / 2);
    let fixed = |d: usize| {
        let nd = d / 3;
        let (i, j, k) = (
            nd % nodes[0],
            (nd / nodes[0]) % nodes[1],
            nd / (nodes[0] * nodes[1]),
        );
        let component = d % 3;
        i == 0
            && (component == 0 || (component == 2 && k == mid_k) || (component == 1 && j == mid_j))
    };
    let free: Vec<usize> = (0..dofs).filter(|&d| !fixed(d)).collect();
    let m = free.len();
    let kf: Vec<f64> = (0..m * m)
        .map(|i| k[free[i / m] * dofs + free[i % m]])
        .collect();
    let ff: Vec<f64> = free.iter().map(|&d| f[d]).collect();
    let u = Cholesky::factor(&kf, m, 0.0)
        .and_then(|c| c.solve(&ff))
        .unwrap_or_else(|e| panic!("{e}"));
    let mut tip = 0.0;
    let mut count = 0.0;
    for (slot, &d) in free.iter().enumerate() {
        if d % 3 == 2 && (d / 3) % nodes[0] == n[0] {
            tip += u[slot];
            count += 1.0;
        }
    }
    -tip / count
}

#[test]
fn c4_a_nearly_incompressible_cantilever_does_not_lock() {
    let nu = 0.49;
    let (p, l, e) = (1e-3, 10.0, 1.0);
    let g = e / (2.0 * (1.0 + nu));
    let exact = p * l * l * l / (3.0 * e / 12.0) + p * l / (5.0 / 6.0 * g);
    let eas = cantilever_tip(nu, 0.5, true);
    assert!((eas - exact).abs() <= 0.05 * exact, "EAS {eas} vs {exact}");
    // Negative control: plain Q1 locks.
    let plain = cantilever_tip(nu, 0.5, false);
    assert!(plain < 0.8 * exact, "plain Q1 {plain} vs {exact}");
}
