use super::pcg::pcg;
use super::*;
use crate::aggregate::aggregate;
use crate::clip::clip_solid;
use crate::element::{ScalarTables, surface_moments};
use crate::grid::Grid;
use crate::linalg::spd_test_helpers::bits_hash;
use crate::mesh_test_helpers::{box_mesh, rotated};
use crate::model::{SolidId, SolidMesh};
use crate::moments::rasterise_solid;

/// A reduced scalar system K + Robin(c = 1) with unit volume source,
/// assembled over `meshes` (one solid each, block diagonal), with the
/// solid of each free dof.
struct System {
    a: Csr,
    b: Vec<f64>,
    positions: Vec<[usize; 3]>,
    solid: Vec<usize>,
}

fn system(grid: &Grid, meshes: &[SolidMesh]) -> System {
    let t = ScalarTables::new();
    let h = grid.cell_size();
    let mut triplets = Vec::new();
    let mut b = Vec::new();
    let mut positions = Vec::new();
    let mut solid_of = Vec::new();
    let mut offset = 0;
    for (s, mesh) in meshes.iter().enumerate() {
        let c = clip_solid(grid, mesh).unwrap_or_else(|e| panic!("{e}"));
        let area = c.pieces.iter().map(|p| p.area).sum();
        let m = rasterise_solid(grid, SolidId(s as u32), &c, 2).unwrap_or_else(|e| panic!("{e}"));
        let agg = aggregate(grid, &m, 0.25, area).unwrap_or_else(|e| panic!("{e}"));
        let e = &agg.extension;
        let mut node_triplets = Vec::new();
        let mut node_load = vec![0.0; e.nodes.len()];
        for (entry, cell) in m.cells.iter().enumerate() {
            let ijk = grid.cell_ijk(cell.cell).unwrap_or_default();
            let nodes = m
                .node_keys(grid, entry)
                .unwrap_or_default()
                .map(|n| e.position(n).unwrap_or(0));
            let start = c.pieces.partition_point(|p| p.cell < cell.cell);
            let end = c.pieces.partition_point(|p| p.cell <= cell.cell);
            let sm = surface_moments(&c.pieces[start..end], grid.cell_min(ijk), h);
            let k = t.conduction(&cell.m, 1.0, h);
            let robin = t.facet_mass(&sm);
            let f = t.volume_load(&cell.m);
            for i in 0..8 {
                node_load[nodes[i]] += f[i];
                for j in 0..8 {
                    node_triplets.push((nodes[i], nodes[j], k[i][j] + robin[i][j]));
                }
            }
        }
        let a_nodes = Csr::from_triplets(e.nodes.len(), e.nodes.len(), &node_triplets)
            .unwrap_or_else(|e| panic!("{e}"));
        let ecsr = e.to_csr().unwrap_or_else(|e| panic!("{e}"));
        let reduced = ecsr
            .transpose()
            .mul(&a_nodes.mul(&ecsr).unwrap_or_else(|e| panic!("{e}")))
            .unwrap_or_else(|e| panic!("{e}"));
        for r in 0..reduced.rows() {
            for (cc, v) in reduced.row(r) {
                triplets.push((offset + r, offset + cc, v));
            }
        }
        b.extend(e.restrict(&node_load));
        positions.extend(e.free_positions(grid));
        solid_of.extend(std::iter::repeat_n(s, e.free));
        offset += e.free;
    }
    System {
        a: Csr::from_triplets(offset, offset, &triplets).unwrap_or_else(|e| panic!("{e}")),
        b,
        positions,
        solid: solid_of,
    }
}

fn n_of(b: &[f64]) -> usize {
    b.len()
}

fn options(split: bool) -> HierarchyOptions {
    HierarchyOptions {
        coarsest_dofs: 500,
        split_components: split,
    }
}

fn cube_grid(h: f64) -> (Grid, SolidMesh) {
    let mesh = rotated(box_mesh([0.0; 3], [3.0; 3]), 17.0, [1.5; 3]);
    let g = Grid::covering([-0.6; 3], [3.6; 3], h).unwrap_or_else(|e| panic!("{e}"));
    (g, mesh)
}

#[test]
fn iterations_are_mesh_independent() {
    let mut iterations = Vec::new();
    for h in [0.3, 0.2, 0.12] {
        let (g, mesh) = cube_grid(h);
        let s = system(&g, &[mesh]);
        let hier = Hierarchy::build(s.a, s.positions, vec![0; n_of(&s.b)], &options(true))
            .unwrap_or_else(|e| panic!("{e}"));
        assert!(hier.levels.len() >= 2, "h {h}: no coarsening");
        assert_eq!(hier.jacobi_fallbacks, 0);
        let r = pcg(&hier, &s.b, 1e-8, 200).unwrap_or_else(|e| panic!("{e}"));
        assert!(r.converged, "h {h}: {:?}", r.history.last());
        iterations.push(r.iterations);
    }
    assert!(iterations.iter().all(|&i| i < 60), "{iterations:?}");
    assert!(
        iterations[2] as f64 <= 1.2 * iterations[0] as f64,
        "iterations grow with refinement: {iterations:?}"
    );
}

/// FNV-1a of the solution bits on the h = 0.3 cube.
const PINNED_SOLUTION_HASH: u64 = 6364974095852511351;

#[test]
fn the_solution_bits_are_pinned() {
    let (g, mesh) = cube_grid(0.3);
    let solve = || {
        let s = system(&g, std::slice::from_ref(&mesh));
        let hier = Hierarchy::build(s.a, s.positions, vec![0; n_of(&s.b)], &options(true))
            .unwrap_or_else(|e| panic!("{e}"));
        pcg(&hier, &s.b, 1e-8, 200)
            .unwrap_or_else(|e| panic!("{e}"))
            .x
    };
    let x = solve();
    assert_eq!(bits_hash(&x), bits_hash(&solve()));
    assert_eq!(bits_hash(&x), PINNED_SOLUTION_HASH);
}

/// Solids of the finest dofs each coarse dof of `level` reaches.
fn solids_per_dof(hier: &Hierarchy, solid: &[usize], level: usize) -> Vec<Vec<usize>> {
    let mut sets: Vec<Vec<usize>> = solid.iter().map(|&s| vec![s]).collect();
    for l in 0..level {
        let p = hier.levels[l]
            .p
            .as_ref()
            .unwrap_or_else(|| panic!("level {l} has no P"));
        let mut next = vec![Vec::new(); p.cols()];
        for (fine, set) in sets.iter().enumerate() {
            for (c, _) in p.row(fine) {
                next[c].extend(set);
            }
        }
        for s in &mut next {
            s.sort_unstable();
            s.dedup();
        }
        sets = next;
    }
    sets
}

#[test]
fn separate_solids_stay_decoupled_on_coarse_levels() {
    let h = 0.2;
    let g = Grid::covering([0.0; 3], [3.0, 1.2, 1.2], h).unwrap_or_else(|e| panic!("{e}"));
    let a = box_mesh([0.05, 0.1, 0.1], [1.5, 1.1, 1.1]);
    let b = box_mesh([1.5 + 0.3 * h, 0.1, 0.1], [2.9, 1.1, 1.1]);
    let s = system(&g, &[a, b]);
    let solid = s.solid.clone();
    let split = Hierarchy::build(
        s.a.clone(),
        s.positions.clone(),
        vec![0; n_of(&s.b)],
        &options(true),
    )
    .unwrap_or_else(|e| panic!("{e}"));
    for level in 1..split.levels.len() {
        let sets = solids_per_dof(&split, &solid, level);
        assert!(
            sets.iter().all(|s| s.len() == 1),
            "level {level} mixes solids"
        );
    }
    let joined = Hierarchy::build(s.a, s.positions, vec![0; n_of(&s.b)], &options(false))
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(joined.levels.len() >= 2);
    let sets = solids_per_dof(&joined, &solid, 1);
    assert!(
        sets.iter().any(|s| s.len() == 2),
        "without splitting a coarse dof joins the solids"
    );
}

#[test]
fn an_unconverged_solve_reports_its_history_and_best_iterate() {
    let (g, mesh) = cube_grid(0.3);
    let s = system(&g, &[mesh]);
    let hier = Hierarchy::build(s.a, s.positions, vec![0; n_of(&s.b)], &options(true))
        .unwrap_or_else(|e| panic!("{e}"));
    let r = pcg(&hier, &s.b, 1e-14, 2).unwrap_or_else(|e| panic!("{e}"));
    assert!(!r.converged);
    assert_eq!(r.iterations, 2);
    assert_eq!(r.history.len(), 3);
    assert!(r.history[1..].iter().all(|&v| v > 1e-14));
}

#[test]
fn a_stalled_coarsening_declines_only_for_a_connected_block_above_the_cap() {
    // Dofs at isolated nodes far apart: no two share a coarse node, so
    // every level keeps all dofs and coarsening stalls.
    let n = COARSEST_HARD_CAP + 10;
    let positions: Vec<[usize; 3]> = (0..n).map(|i| [64 * i, 0, 0]).collect();
    // Uncoupled dofs (separate solids): one-dof blocks, factored apart.
    let diagonal: Vec<(usize, usize, f64)> = (0..n).map(|i| (i, i, 1.0)).collect();
    let a = Csr::from_triplets(n, n, &diagonal).unwrap_or_else(|e| panic!("{e}"));
    let h = Hierarchy::build(a, positions.clone(), vec![0; n], &options(true))
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(h.coarsest.block_count(), n);
    // A chain couples them into one block above the cap: declined.
    let mut chain = diagonal.clone();
    for t in &mut chain {
        t.2 = 3.0;
    }
    for i in 1..n {
        chain.push((i - 1, i, -1.0));
        chain.push((i, i - 1, -1.0));
    }
    let a = Csr::from_triplets(n, n, &chain).unwrap_or_else(|e| panic!("{e}"));
    let r = Hierarchy::build(a, positions, vec![0; n], &options(true));
    assert_eq!(
        r.err(),
        Some(MultigridError::CoarsestTooLarge {
            dofs: n,
            cap: COARSEST_HARD_CAP
        })
    );
}

/// A rescaled operator (a design update) rebuilt against the previous
/// hierarchy is the hierarchy built afresh, level by level, and solves
/// to the same bits; with other labels (the negative control) the
/// prolongations differ and are not reused.
#[test]
fn a_rebuilt_hierarchy_is_the_fresh_one() {
    let h = 0.2;
    let g = Grid::covering([0.0; 3], [3.0, 1.2, 1.2], h).unwrap_or_else(|e| panic!("{e}"));
    let a = box_mesh([0.05, 0.1, 0.1], [1.5, 1.1, 1.1]);
    let b = box_mesh([1.5 + 0.3 * h, 0.1, 0.1], [2.9, 1.1, 1.1]);
    let s = system(&g, &[a, b]);
    let n = n_of(&s.b);
    let zeros = vec![0u8; n];
    let opts = HierarchyOptions {
        coarsest_dofs: 40,
        split_components: true,
    };
    let first = Hierarchy::build(s.a.clone(), s.positions.clone(), zeros.clone(), &opts)
        .unwrap_or_else(|e| panic!("{e}"));
    let mut scaled = s.a.clone();
    for v in scaled.values_mut() {
        *v *= 1.7;
    }
    let build = |labels: &[u8], previous: Option<&Hierarchy>| {
        Hierarchy::rebuild(
            scaled.clone(),
            s.positions.clone(),
            labels.to_vec(),
            &opts,
            previous,
        )
        .unwrap_or_else(|e| panic!("{e}"))
    };
    let fresh = build(&zeros, None);
    let rebuilt = build(&zeros, Some(&first));
    assert!(fresh.levels.len() >= 3);
    assert_eq!(fresh.levels.len(), rebuilt.levels.len());
    for (f, r) in fresh.levels.iter().zip(&rebuilt.levels) {
        assert_eq!(f.a, r.a);
        assert_eq!(f.p, r.p);
        assert_eq!(f.r, r.r);
        assert_eq!(f.smoother, r.smoother);
    }
    let x = |h: &Hierarchy| pcg(h, &s.b, 1e-8, 200).unwrap_or_else(|e| panic!("{e}")).x;
    assert_eq!(bits_hash(&x(&fresh)), bits_hash(&x(&rebuilt)));
    let other: Vec<u8> = (0..n).map(|i| (i % 3) as u8).collect();
    let relabelled = build(&other, Some(&first));
    assert_ne!(relabelled.levels[0].p, first.levels[0].p);
    assert_eq!(relabelled.levels[0].p, build(&other, None).levels[0].p);
}
