use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

use super::{
    Kp20SupportCellConstructionMode, Kp20SupportCellPrimitiveTrial, gregory_from_bi6,
    support_cell_topology_path_g0, support_cells,
};

#[test]
fn topology_path_g0_measures_rectangular_support_grid_adjacency() {
    let base = flat_bi6_patch();
    let patches = support_cells(&base, 2, 4)
        .into_iter()
        .map(gregory_from_bi6)
        .collect::<Vec<_>>();

    let g0 = support_cell_topology_path_g0(&patches, 2, 4);

    assert!(
        g0 <= 1.0e-10,
        "2x4 support-cell grid should measure exact internal adjacency, got {g0:.3e}"
    );
}

#[test]
fn topology_path_g0_declines_incomplete_rectangular_inventory() {
    let base = flat_bi6_patch();
    let mut patches = support_cells(&base, 2, 4)
        .into_iter()
        .map(gregory_from_bi6)
        .collect::<Vec<_>>();
    patches.pop();

    assert!(support_cell_topology_path_g0(&patches, 2, 4).is_infinite());
}

#[test]
fn support_cell_trial_prefers_commit_supported_inventory() {
    let supported = trial(2, 4, 48, 0, 1.0e-6);
    let unsupported = trial(3, 3, 54, 0, 0.0);

    assert!(supported.is_better_than(&unsupported));
}

#[test]
fn support_cell_trial_keeps_lower_supported_resolution_on_quality_tie() {
    let eight_cell = trial(2, 4, 48, 0, 0.0);
    let four_cell = trial(1, 4, 24, 0, 0.0);

    assert!(four_cell.is_better_than(&eight_cell));
}

#[test]
fn support_cell_trial_promotes_higher_resolution_on_measured_g0_improvement() {
    let eight_cell = trial(2, 4, 48, 0, 1.0e-10);
    let four_cell = trial(1, 4, 24, 0, 2.0e-9);

    assert!(eight_cell.is_better_than(&four_cell));
}

#[test]
fn support_cell_trial_promotes_sixteen_cell_resolution_on_measured_g0_improvement() {
    let sixteen_cell = trial(4, 4, 96, 0, 1.0e-10);
    let eight_cell = trial(2, 4, 48, 0, 2.0e-9);

    assert!(sixteen_cell.is_better_than(&eight_cell));
}

fn trial(
    cells_u: usize,
    cells_v: usize,
    patch_count: usize,
    folded_patch_count: usize,
    topology_path_g0: f64,
) -> Kp20SupportCellPrimitiveTrial {
    Kp20SupportCellPrimitiveTrial {
        patch_count,
        folded_patch_count,
        cells_u,
        cells_v,
        topology_path_g0,
        construction_mode: Kp20SupportCellConstructionMode::InheritedCellBoundary,
        emitted_patches: Vec::new(),
        emitted_provenance: Vec::new(),
    }
}

fn flat_bi6_patch() -> Kp20Bi6Patch {
    let mut patch = Kp20Bi6Patch::zero();
    for u in 0..7 {
        for v in 0..7 {
            patch.cps[u][v] = Point3::new(u as f64 / 6.0, v as f64 / 6.0, 0.0);
        }
    }
    patch
}
