use std::collections::BTreeMap;

use super::poly::{eval, shape, shape_derivative};
use super::scalar::{Matrix8, ScalarTables, surface_moments};
use crate::clip::{ClippedSolid, clip_solid};
use crate::grid::{CELL_CORNERS, Grid};
use crate::linalg::Cholesky;
use crate::mesh_test_helpers::{bounds, box_mesh, rotated};
use crate::model::{SolidId, SolidMesh};
use crate::moments::{SolidMoments, moment_index, rasterise_solid};

/// Moments of a full cell of size h.
fn full_moments(h: f64) -> [f64; 27] {
    let axis = [2.0, 0.0, 2.0 / 3.0];
    let s = 0.125 * h * h * h;
    std::array::from_fn(|i| s * axis[i % 3] * axis[(i / 3) % 3] * axis[i / 9])
}

#[test]
fn the_full_cell_conduction_matrix_has_only_the_constant_null_mode() {
    let t = ScalarTables::new();
    let h = 0.7;
    let k = 2.5;
    let km = t.conduction(&full_moments(h), k, h);
    let scale = km.iter().flatten().fold(0.0f64, |a, x| a.max(x.abs()));
    for (a, row) in km.iter().enumerate() {
        assert!(row.iter().sum::<f64>().abs() <= 1e-14 * scale, "row {a}");
        for (b, x) in row.iter().enumerate() {
            assert_eq!(*x, km[b][a]);
        }
        // Hex8 Laplacian diagonal: k h / 3.
        assert!((row[a] - k * h / 3.0).abs() <= 1e-15 * scale);
    }
    // K + c 1·1ᵀ is SPD exactly when the constant is the only null mode.
    let flat: Vec<f64> = (0..64).map(|i| km[i / 8][i % 8] + scale).collect();
    assert!(Cholesky::factor(&flat, 8, 1e-12).is_ok());
}

#[test]
fn shape_functions_interpolate_the_corners() {
    for a in 0..8 {
        for (b, corner) in CELL_CORNERS.iter().enumerate() {
            let x = corner.map(|o| if o == 0 { -1.0 } else { 1.0 });
            let want = if a == b { 1.0 } else { 0.0 };
            assert_eq!(eval(&shape(a), x), want);
        }
    }
}

fn setup(mesh: &SolidMesh, h: f64, shift: f64) -> (Grid, ClippedSolid, SolidMoments) {
    let (lo, hi) = bounds(mesh);
    let g = Grid::covering(lo.map(|c| c - shift), hi, h).unwrap_or_else(|e| panic!("{e}"));
    let c = clip_solid(&g, mesh).unwrap_or_else(|e| panic!("{e}"));
    let r = rasterise_solid(&g, SolidId(0), &c, 2).unwrap_or_else(|e| panic!("{e}"));
    (g, c, r)
}

/// ∫ over the physical box [lo, hi] (inside one cell starting at `cmin`)
/// of f(local ξ) by 2-point Gauss per axis: exact for degree ≤ 3.
fn gauss(lo: [f64; 3], hi: [f64; 3], cmin: [f64; 3], h: f64, f: impl Fn([f64; 3]) -> f64) -> f64 {
    let r = 1.0 / 3f64.sqrt();
    let mut sum = 0.0;
    for gz in [-r, r] {
        for gy in [-r, r] {
            for gx in [-r, r] {
                let g = [gx, gy, gz];
                let x = [0, 1, 2].map(|a| 0.5 * (lo[a] + hi[a]) + 0.5 * (hi[a] - lo[a]) * g[a]);
                sum += f([0, 1, 2].map(|a| 2.0 * (x[a] - cmin[a]) / h - 1.0));
            }
        }
    }
    sum * (0..3).map(|a| 0.5 * (hi[a] - lo[a])).product::<f64>()
}

#[test]
fn moment_tables_match_gauss_on_box_cut_cells() {
    let t = ScalarTables::new();
    let mut rng = crate::linalg::spd_test_helpers::Lcg(77);
    let mut checked = 0;
    for _ in 0..4 {
        let min = [0, 1, 2].map(|_| rng.next_signed());
        let max = min.map(|m| m + 1.3 + 0.8 * (rng.next_signed() + 1.0));
        let mesh = box_mesh(min, max);
        let h = 0.55;
        let (g, _, r) = setup(&mesh, h, 0.17);
        for cell in r
            .cells
            .iter()
            .filter(|c| c.kind == crate::moments::CellKind::Cut)
            .take(5)
        {
            let ijk = g.cell_ijk(cell.cell).unwrap_or_default();
            let cmin = g.cell_min(ijk);
            let lo = [0, 1, 2].map(|a| cmin[a].max(min[a]));
            let hi = [0, 1, 2].map(|a| (cmin[a] + h).min(max[a]));
            let km = t.conduction(&cell.m, 1.0, h);
            let mm = t.mass(&cell.m);
            let fv = t.volume_load(&cell.m);
            for a in 0..8 {
                let want_f = gauss(lo, hi, cmin, h, |x| eval(&shape(a), x));
                assert!((fv[a] - want_f).abs() <= 1e-12 * h * h * h);
                for b in 0..8 {
                    let want_k = gauss(lo, hi, cmin, h, |x| {
                        (0..3)
                            .map(|d| {
                                eval(&shape_derivative(a, d), x) * eval(&shape_derivative(b, d), x)
                            })
                            .sum::<f64>()
                            * 4.0
                            / (h * h)
                    });
                    let want_m =
                        gauss(lo, hi, cmin, h, |x| eval(&shape(a), x) * eval(&shape(b), x));
                    assert!((km[a][b] - want_k).abs() <= 1e-12 * h, "K {a}{b}");
                    assert!((mm[a][b] - want_m).abs() <= 1e-12 * h * h * h, "M {a}{b}");
                }
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 20);
}

/// Global residual K·T − ∮ N k∇T·n over all active cells, for linear
/// T = t0 + g·x, using `matrix(cell moments)` as the cell operator.
fn patch_residual(
    mesh: &SolidMesh,
    h: f64,
    k: f64,
    matrix: &dyn Fn(&ScalarTables, &[f64; 27]) -> Matrix8,
) -> (BTreeMap<usize, f64>, f64, f64) {
    let t = ScalarTables::new();
    let (g, c, r) = setup(mesh, h, 0.11);
    let grad = [0.3, -1.2, 0.7];
    let field = |x: [f64; 3]| 2.0 + grad[0] * x[0] + grad[1] * x[1] + grad[2] * x[2];
    let mut residual: BTreeMap<usize, f64> = BTreeMap::new();
    let mut energy = 0.0;
    for cell in &r.cells {
        let ijk = g.cell_ijk(cell.cell).unwrap_or_default();
        let nodes = g.cell_nodes(ijk).unwrap_or_default();
        let temps: [f64; 8] = std::array::from_fn(|a| {
            let n = g.node_ijk(nodes[a]).unwrap_or_default();
            field(g.node_position(n))
        });
        let km = matrix(&t, &cell.m);
        let start = c.pieces.partition_point(|p| p.cell < cell.cell);
        let end = c.pieces.partition_point(|p| p.cell <= cell.cell);
        let s = surface_moments(&c.pieces[start..end], g.cell_min(ijk), h);
        let flux = t.normal_load(&s, grad.map(|x| k * x));
        for a in 0..8 {
            let kt: f64 = (0..8).map(|b| km[a][b] * temps[b]).sum();
            *residual.entry(nodes[a]).or_insert(0.0) += kt - flux[a];
            energy += temps[a] * kt;
        }
    }
    let exact_energy = k * grad.iter().map(|x| x * x).sum::<f64>() * r.volume;
    (residual, energy, exact_energy)
}

#[test]
fn a_linear_temperature_passes_the_cut_cell_patch_test() {
    let mesh = rotated(box_mesh([0.0; 3], [3.0; 3]), 17.0, [1.5; 3]);
    let k = 1.7;
    let h = 0.6;
    let exact = |t: &ScalarTables, m: &[f64; 27]| t.conduction(m, k, h);
    let (res, energy, want) = patch_residual(&mesh, h, k, &exact);
    let worst = res.values().fold(0.0f64, |a, r| a.max(r.abs()));
    // Nodal fluxes are of order k |g| h²; demand round-off.
    assert!(worst <= 1e-12 * k * 1.4 * h * h, "worst residual {worst}");
    assert!((energy - want).abs() <= 1e-12 * want, "{energy} vs {want}");
    // Negative control: the volume-fraction-scaled full-cell matrix
    // (a smeared boundary) fails at the cut nodes.
    let full = full_moments(h);
    let smeared = |t: &ScalarTables, m: &[f64; 27]| {
        let eta = m[moment_index(0, 0, 0)] / (h * h * h);
        t.conduction(&full, k, h).map(|row| row.map(|x| eta * x))
    };
    let (res, _, _) = patch_residual(&mesh, h, k, &smeared);
    let worst = res.values().fold(0.0f64, |a, r| a.max(r.abs()));
    assert!(worst > 1e-3 * k * h * h, "smeared operator passed: {worst}");
}
