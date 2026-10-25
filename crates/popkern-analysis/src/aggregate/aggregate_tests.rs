use super::*;
use crate::clip::clip_solid;
use crate::element::ScalarTables;
use crate::linalg::Cholesky;
use crate::linalg::spd_test_helpers::Lcg;
use crate::mesh_test_helpers::box_mesh;
use crate::model::SolidMesh;
use crate::moments::rasterise_solid;

const H: f64 = 0.25;

fn moments_on(grid: &Grid, mesh: &SolidMesh, solid: u32) -> (SolidMoments, f64) {
    let c = clip_solid(grid, mesh).unwrap_or_else(|e| panic!("{e}"));
    let area = c.pieces.iter().map(|p| p.area).sum();
    let r = rasterise_solid(grid, SolidId(solid), &c, 2).unwrap_or_else(|e| panic!("{e}"));
    (r, area)
}

fn grid() -> Grid {
    Grid::covering([0.0; 3], [2.5, 1.0, 1.0], H).unwrap_or_else(|e| panic!("{e}"))
}

#[test]
fn the_extension_reproduces_linear_fields_and_is_adjoint() {
    // A bar ending in a 2% sliver layer, so some nodes are ill-posed.
    let mesh = box_mesh([0.05, 0.1, 0.12], [1.0 + 0.02 * H, 0.6, 0.62]);
    let g = grid();
    let (m, area) = moments_on(&g, &mesh, 0);
    let agg = aggregate(&g, &m, 0.25, area).unwrap_or_else(|e| panic!("{e}"));
    let e = &agg.extension;
    assert!(
        e.rows.iter().any(|r| matches!(r, Row::Extended(_))),
        "no ill-posed node"
    );
    assert!(agg.cells_aggregated > 0);
    // Free values of a linear field, extended, equal the field.
    let field = |x: [f64; 3]| 1.5 - 0.7 * x[0] + 2.0 * x[1] + 0.3 * x[2];
    let at = |node: usize| {
        field(
            g.node_position(
                g.node_ijk(node / crate::moments::fragments::NODE_COPIES)
                    .unwrap_or_default(),
            ),
        )
    };
    let free: Vec<f64> = e
        .nodes
        .iter()
        .zip(&e.rows)
        .filter(|(_, r)| matches!(r, Row::Free(_)))
        .map(|(n, _)| at(*n))
        .collect();
    assert_eq!(free.len(), e.free);
    for (node, v) in e.nodes.iter().zip(e.expand(&free)) {
        assert!((v - at(*node)).abs() < 1e-12, "node {node}");
    }
    // <E x, y> = <x, Eᵀ y>
    let mut rng = Lcg(5);
    let x: Vec<f64> = (0..e.free).map(|_| rng.next_signed()).collect();
    let y: Vec<f64> = (0..e.nodes.len()).map(|_| rng.next_signed()).collect();
    let lhs: f64 = e.expand(&x).iter().zip(&y).map(|(a, b)| a * b).sum();
    let rhs: f64 = x.iter().zip(e.restrict(&y)).map(|(a, b)| a * b).sum();
    assert!((lhs - rhs).abs() < 1e-12);
}

#[test]
fn aggregation_is_deterministic_and_stays_within_its_solid() {
    // Two boxes 0.3h apart sharing a column of cells: B's small cells
    // must root in B even where A's large cells are nearer.
    let a = box_mesh([0.05, 0.1, 0.1], [1.0, 0.9, 0.9]);
    let b = box_mesh([1.0 + 0.3 * H, 0.1, 0.1], [1.0 + 0.3 * H + 0.3, 0.9, 0.9]);
    let g = grid();
    let (ma, area_a) = moments_on(&g, &a, 0);
    let (mb, area_b) = moments_on(&g, &b, 1);
    let first = aggregate(&g, &mb, 0.25, area_b).unwrap_or_else(|e| panic!("{e}"));
    let second = aggregate(&g, &mb, 0.25, area_b).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(first, second);
    let b_large: Vec<usize> = mb
        .cells
        .iter()
        .filter(|c| c.volume_fraction(H) >= 0.25)
        .map(|c| c.cell)
        .collect();
    // Roots are entries of the solid's own moments.
    assert!(
        first
            .root
            .iter()
            .all(|&r| b_large.contains(&mb.cells[r].cell))
    );
    let agg_a = aggregate(&g, &ma, 0.25, area_a).unwrap_or_else(|e| panic!("{e}"));
    assert!(agg_a.root.iter().all(|&r| r < ma.cells.len()));
    assert_eq!(first.unreached, 0);
}

/// Dense K + M over active nodes, plain or reduced by E.
fn system(g: &Grid, m: &SolidMoments, agg: &Aggregation, reduced: bool) -> (Vec<f64>, usize) {
    let t = ScalarTables::new();
    let e = &agg.extension;
    let n = e.nodes.len();
    let mut a = vec![0.0; n * n];
    for (entry, cell) in m.cells.iter().enumerate() {
        let nodes = m
            .node_keys(g, entry)
            .unwrap_or_default()
            .map(|v| e.position(v).unwrap_or(0));
        let k = t.conduction(&cell.m, 1.0, H);
        let mm = t.mass(&cell.m);
        for i in 0..8 {
            for j in 0..8 {
                a[nodes[i] * n + nodes[j]] += k[i][j] + mm[i][j];
            }
        }
    }
    if !reduced {
        return (a, n);
    }
    // Eᵀ A E column by column.
    let f = e.free;
    let mut r = vec![0.0; f * f];
    for col in 0..f {
        let mut unit = vec![0.0; f];
        unit[col] = 1.0;
        let ex = e.expand(&unit);
        let aex: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| a[i * n + j] * ex[j]).sum())
            .collect();
        for (row, v) in e.restrict(&aex).into_iter().enumerate() {
            r[row * f + col] = v;
        }
    }
    (r, f)
}

/// Condition number by power and inverse iteration.
fn condition(a: &[f64], n: usize) -> f64 {
    let mul = |x: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
            .collect()
    };
    let norm = |x: &[f64]| x.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut rng = Lcg(9);
    let start: Vec<f64> = (0..n).map(|_| 1.0 + 0.1 * rng.next_signed()).collect();
    let mut x = start.clone();
    let mut lmax = 0.0;
    for _ in 0..400 {
        let y = mul(&x);
        lmax = norm(&y) / norm(&x);
        let s = norm(&y);
        x = y.into_iter().map(|v| v / s).collect();
    }
    let chol = Cholesky::factor(a, n, 0.0).unwrap_or_else(|e| panic!("{e}"));
    let mut x = start;
    let mut lmin = 0.0;
    for _ in 0..200 {
        let y = chol.solve(&x).unwrap_or_else(|e| panic!("{e}"));
        lmin = norm(&x) / norm(&y);
        let s = norm(&y);
        x = y.into_iter().map(|v| v / s).collect();
    }
    lmax / lmin
}

#[test]
fn agfem_keeps_the_condition_number_bounded_over_sliver_offsets() {
    let g = grid();
    let mut reduced = Vec::new();
    let mut plain = Vec::new();
    for i in 0..50 {
        // Sliver depths from 1e-1 to 1e-8 of a cell, log-spaced.
        let depth = 10f64.powf(-1.0 - 7.0 * i as f64 / 49.0);
        let mesh = box_mesh([0.05, 0.1, 0.12], [1.0 + depth * H, 0.6, 0.62]);
        let (m, area) = moments_on(&g, &mesh, 0);
        let agg = aggregate(&g, &m, 0.25, area).unwrap_or_else(|e| panic!("{e}"));
        let (a, n) = system(&g, &m, &agg, true);
        reduced.push(condition(&a, n));
        let (a, n) = system(&g, &m, &agg, false);
        plain.push(condition(&a, n));
    }
    let (lo, hi) = reduced
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(l, h), &k| (l.min(k), h.max(k)));
    assert!(hi / lo < 2.0, "AgFEM condition varies {lo}..{hi}");
    // Negative control: without aggregation the deepest sliver blows up.
    assert!(plain[49] > 1e4 * hi, "plain {} vs reduced {hi}", plain[49]);
}

#[test]
fn a_wall_thinner_than_a_cell_is_under_resolved() {
    let g = grid();
    let thin = box_mesh([0.1, 0.1, 0.1], [2.1, 0.9, 0.15]);
    let (m, area) = moments_on(&g, &thin, 0);
    let agg = aggregate(&g, &m, 0.25, area).unwrap_or_else(|e| panic!("{e}"));
    let Some(Warning::UnderResolved { wall, cell, .. }) = agg.warnings.first() else {
        panic!("no UnderResolved warning: {:?}", agg.warnings);
    };
    assert_eq!(*cell, H);
    assert!((wall - 0.05).abs() < 0.1 * 0.05, "wall {wall}");
    assert!(agg.gamma < 0.25);
    // Negative control: a wall of two cells is resolved.
    let thick = box_mesh([0.1, 0.1, 0.1], [2.1, 0.9, 0.6]);
    let (m, area) = moments_on(&g, &thick, 0);
    let agg = aggregate(&g, &m, 0.25, area).unwrap_or_else(|e| panic!("{e}"));
    assert!(agg.warnings.is_empty(), "{:?}", agg.warnings);
}
