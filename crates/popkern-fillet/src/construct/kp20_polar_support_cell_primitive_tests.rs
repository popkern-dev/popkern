use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

use super::{gregory_from_bi6, support_cell_topology_path_g0, support_cells};

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

fn flat_bi6_patch() -> Kp20Bi6Patch {
    let mut patch = Kp20Bi6Patch::zero();
    for u in 0..7 {
        for v in 0..7 {
            patch.cps[u][v] = Point3::new(u as f64 / 6.0, v as f64 / 6.0, 0.0);
        }
    }
    patch
}
