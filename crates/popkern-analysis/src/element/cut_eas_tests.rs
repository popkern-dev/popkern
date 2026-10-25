use std::collections::BTreeMap;

use super::scalar::{ScalarTables, surface_moments};
use super::vector::*;
use crate::clip::clip_solid;
use crate::grid::Grid;
use crate::mesh_test_helpers::{box_mesh, rotated};
use crate::model::SolidId;
use crate::moments::{CellKind, SolidMoments, rasterise_solid};

const E: f64 = 2.5;
const NU: f64 = 0.3;

/// A 17°-rotated 3-mm block on cells of 0.45: every surface cell is cut.
fn block() -> (Grid, crate::clip::ClippedSolid, SolidMoments) {
    let mesh = rotated(box_mesh([0.0; 3], [3.0; 3]), 17.0, [1.5; 3]);
    let g = Grid::covering([-0.7; 3], [3.7; 3], 0.45).unwrap_or_else(|e| panic!("{e}"));
    let c = clip_solid(&g, &mesh).unwrap_or_else(|e| panic!("{e}"));
    let r = rasterise_solid(&g, SolidId(0), &c, 2).unwrap_or_else(|e| panic!("{e}"));
    (g, c, r)
}

/// A linear displacement field's gradient, strain and stress.
fn linear_state() -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let grad = [
        [1.0e-3, 2.0e-4, -3.0e-4],
        [5.0e-4, -8.0e-4, 1.0e-4],
        [-2.0e-4, 6.0e-4, 4.0e-4],
    ];
    let (l, m) = lame(E, NU);
    let eps: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| 0.5 * (grad[i][j] + grad[j][i])));
    let tr = eps[0][0] + eps[1][1] + eps[2][2];
    let sigma = std::array::from_fn(|i| {
        std::array::from_fn(|j| 2.0 * m * eps[i][j] + if i == j { l * tr } else { 0.0 })
    });
    (grad, sigma)
}

/// Global residual K u_lin − ∮ N σ·n over every dof, the strain energy
/// ½ uᵀKu, and the number of cells that skipped EAS.
fn patch(enhancement: Option<Enhancement>) -> (f64, f64, f64, usize, usize) {
    let (g, c, r) = block();
    let (grad, sigma) = linear_state();
    let (l, m) = lame(E, NU);
    let vt = VectorTables::new();
    let st = ScalarTables::new();
    let h = g.cell_size();
    let mut residual: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    let mut energy = 0.0;
    let mut skipped = 0;
    for cell in &r.cells {
        // The control differs from production only in the mean shift:
        // it keeps the same eta_EAS skip.
        let (op, used) = match enhancement {
            Some(e) if cell.volume_fraction(h) >= ETA_EAS => (
                vt.cell(&cell.m, h, l, m, e)
                    .unwrap_or_else(|e| panic!("{e}")),
                true,
            ),
            _ => vt
                .cell_operator(&cell.m, h, l, m)
                .unwrap_or_else(|e| panic!("{e}")),
        };
        if !used {
            skipped += 1;
        }
        let ijk = g.cell_ijk(cell.cell).unwrap_or_default();
        let nodes = g.cell_nodes(ijk).unwrap_or_default();
        let u: [f64; NODAL_DOFS] = std::array::from_fn(|d| {
            let x = g.node_position(g.node_ijk(nodes[d / 3]).unwrap_or_default());
            (0..3).map(|j| grad[d % 3][j] * x[j]).sum()
        });
        let start = c.pieces.partition_point(|p| p.cell < cell.cell);
        let end = c.pieces.partition_point(|p| p.cell <= cell.cell);
        let s = surface_moments(&c.pieces[start..end], g.cell_min(ijk), h);
        let traction: [[f64; 8]; 3] = std::array::from_fn(|i| st.normal_load(&s, sigma[i]));
        for d in 0..NODAL_DOFS {
            let ku: f64 = (0..NODAL_DOFS).map(|e| op.k[d][e] * u[e]).sum();
            *residual.entry((nodes[d / 3], d % 3)).or_insert(0.0) += ku - traction[d % 3][d / 3];
            energy += 0.5 * u[d] * ku;
        }
    }
    let worst = residual.values().fold(0.0f64, |a, x| a.max(x.abs()));
    let sigma_eps: f64 = (0..3)
        .flat_map(|i| (0..3).map(move |j| (i, j)))
        .map(|(i, j)| sigma[i][j] * 0.5 * (grad[i][j] + grad[j][i]))
        .sum();
    let exact_energy = 0.5 * sigma_eps * r.volume;
    let cut_small = r
        .cells
        .iter()
        .filter(|c| c.kind == CellKind::Cut && c.volume_fraction(h) < ETA_EAS)
        .count();
    (worst, energy, exact_energy, skipped, cut_small)
}

/// Scale of nodal forces: |σ| h².
fn force_scale() -> f64 {
    let (_, sigma) = linear_state();
    let s = sigma.iter().flatten().fold(0.0f64, |a, x| a.max(x.abs()));
    s * 0.45 * 0.45
}

#[test]
fn x3_cut_cells_pass_the_patch_test_with_the_fragment_mean() {
    let (worst, _, _, skipped, cut_small) = patch(None);
    assert!(worst <= 1e-10 * force_scale(), "worst residual {worst}");
    assert_eq!(
        skipped, cut_small,
        "EAS skipped exactly on cells with eta < {ETA_EAS}"
    );
    assert!(cut_small > 0, "the block should produce small cut cells");
    // Negative control: without the mean shift the enhanced modes pick
    // up the linear field and the patch test fails.
    let (raw, _, _, _, _) = patch(Some(Enhancement::Raw));
    assert!(raw > 1e-4 * force_scale(), "raw EAS passed: {raw}");
}

#[test]
fn x4_the_strain_energy_of_a_linear_field_is_exact() {
    let (_, energy, exact, _, _) = patch(None);
    assert!(
        (energy - exact).abs() <= 1e-12 * exact,
        "{energy} vs {exact}"
    );
}
