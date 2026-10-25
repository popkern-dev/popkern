use super::lattice::{contains, lattice_moments};
use super::*;
use crate::clip::clip_solid;
use crate::mesh_test_helpers::{bounds, box_mesh, icosphere, u_prism};
use crate::model::SolidMesh;

fn setup(mesh: &SolidMesh, h: f64, shift: f64) -> (Grid, ClippedSolid) {
    let (lo, hi) = bounds(mesh);
    let g = Grid::covering(lo.map(|c| c - shift), hi, h).unwrap_or_else(|e| panic!("{e}"));
    let c = clip_solid(&g, mesh).unwrap_or_else(|e| panic!("{e}"));
    (g, c)
}

fn query(g: &Grid, c: &ClippedSolid, p: [f64; 3]) -> bool {
    let ijk = g
        .cell_of_point(p)
        .unwrap_or_else(|| panic!("{p:?} outside grid"));
    contains(g, c, ijk, p)
}

#[test]
fn containment_matches_the_geometry() {
    let mesh = box_mesh([0.3, -1.1, 2.0], [4.1, 2.2, 3.7]);
    let (g, c) = setup(&mesh, 0.7, 0.13);
    let mut rng = crate::linalg::spd_test_helpers::Lcg(3);
    for _ in 0..2000 {
        let p = [
            2.2 + 2.5 * rng.next_signed(),
            0.55 + 2.2 * rng.next_signed(),
            2.85 + 1.4 * rng.next_signed(),
        ];
        let inside =
            p[0] > 0.3 && p[0] < 4.1 && p[1] > -1.1 && p[1] < 2.2 && p[2] > 2.0 && p[2] < 3.7;
        assert_eq!(query(&g, &c, p), inside, "box {p:?}");
    }
    // Icosphere of radius 2 (inscribed): decide only away from the
    // faceting band.
    let sphere = icosphere([0.2, 0.1, -0.3], 2.0, 3);
    let (g, c) = setup(&sphere, 0.45, 0.07);
    let mut checked = 0;
    for _ in 0..3000 {
        let p = [
            0.2 + 2.3 * rng.next_signed(),
            0.1 + 2.3 * rng.next_signed(),
            -0.3 + 2.3 * rng.next_signed(),
        ];
        let r = ((p[0] - 0.2).powi(2) + (p[1] - 0.1).powi(2) + (p[2] + 0.3).powi(2)).sqrt();
        if (1.95..=2.0).contains(&r) || g.cell_of_point(p).is_none() {
            continue;
        }
        assert_eq!(query(&g, &c, p), r < 1.95, "sphere {p:?} r {r}");
        checked += 1;
    }
    assert!(checked > 2000);
}

#[test]
fn lattice_moments_converge_to_the_exact_moments() {
    let sphere = icosphere([0.2, 0.1, -0.3], 2.0, 3);
    let (g, c) = setup(&sphere, 0.6, 0.0);
    let exact = rasterise_solid(&g, SolidId(0), &c, 2).unwrap_or_else(|e| panic!("{e}"));
    let h = g.cell_size();
    for depth in [2u8, 3] {
        let spacing = h / f64::from(1u32 << depth);
        for cell in exact.cells.iter().filter(|c| c.kind == CellKind::Cut) {
            let ijk = g.cell_ijk(cell.cell).unwrap_or_default();
            let area: f64 = c
                .pieces
                .iter()
                .filter(|p| p.cell == cell.cell)
                .map(|p| p.area)
                .sum();
            let approx = lattice_moments(&g, &c, ijk, depth);
            // The midpoint rule misclassifies at most the sub-cells the
            // surface passes through: volume ≤ √3 · spacing · area. The
            // exact value carries round-off at 1e-12 h³ (an extreme vertex
            // an ulp past a grid plane leaves a ~1e-29 sliver piece).
            let bound = 3f64.sqrt() * spacing * area + 1e-12 * h * h * h;
            assert!(
                (approx[0] - cell.m[0]).abs() <= bound,
                "depth {depth} cell {ijk:?}: lattice {} exact {} bound {bound} area {area}",
                approx[0],
                cell.m[0]
            );
        }
    }
    // A full cell is exact for the volume.
    let full = exact
        .cells
        .iter()
        .find(|c| c.kind == CellKind::Full)
        .unwrap_or_else(|| panic!("no full cell"));
    let ijk = g.cell_ijk(full.cell).unwrap_or_default();
    assert!((lattice_moments(&g, &c, ijk, 2)[0] - h * h * h).abs() < 1e-15);
}

#[test]
fn forced_fallbacks_are_counted_and_stay_within_the_bound() {
    let sphere = icosphere([0.2, 0.1, -0.3], 2.0, 3);
    let (g, c) = setup(&sphere, 0.6, 0.0);
    let exact = rasterise_solid(&g, SolidId(0), &c, 2).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(exact.fallback_cells, 0, "exact geometry needs no fallback");
    let forced = rasterise_with(
        &g,
        SolidId(0),
        &c,
        3,
        &|_| true,
        crate::moments::Splitting::Exact,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(forced.fallback_cells, exact.cut);
    let spacing = g.cell_size() / 8.0;
    let area: f64 = c.pieces.iter().map(|p| p.area).sum();
    assert!((forced.volume - exact.volume).abs() <= 3f64.sqrt() * spacing * area);
    assert!(forced.volume != exact.volume, "the fallback really ran");
}

#[test]
fn inconsistent_volumes_are_detected() {
    let h = 0.5;
    let mut m = [0.0; 27];
    m[0] = h * h * h;
    assert!(consistent(&m, h));
    m[0] = 2.0 * h * h * h;
    assert!(!consistent(&m, h));
    m[0] = -0.01 * h * h * h;
    assert!(!consistent(&m, h));
    m[0] = f64::NAN;
    assert!(!consistent(&m, h));
}

#[test]
fn a_slot_narrower_than_a_cell_is_detected() {
    let slotted = u_prism(0.3);
    assert!((slotted.enclosed_volume() - (6.0 - 0.3 * 1.5)).abs() < 1e-12);
    // Unit cells on integer planes: the slot (x 1.35..1.65, y > 0.5)
    // splits the cells of column x in [1, 2], y in [1, 2].
    let g = Grid::covering([0.0; 3], [3.0, 2.0, 1.0], 1.0).unwrap_or_else(|e| panic!("{e}"));
    let c = clip_solid(&g, &slotted).unwrap_or_else(|e| panic!("{e}"));
    let r = rasterise_solid(&g, SolidId(0), &c, 2).unwrap_or_else(|e| panic!("{e}"));
    let slits = slit_cells(&g, &c, &r, 2);
    let ijk: Vec<[usize; 3]> = slits.iter().filter_map(|&s| g.cell_ijk(s)).collect();
    assert_eq!(
        ijk,
        [[2, 2, 1]],
        "grid has one padding cell, so x 1..2 is i = 2"
    );
    // Negative control: the same outline without a slot has none.
    let solid = box_mesh([0.0; 3], [3.0, 2.0, 1.0]);
    let c = clip_solid(&g, &solid).unwrap_or_else(|e| panic!("{e}"));
    let r = rasterise_solid(&g, SolidId(0), &c, 2).unwrap_or_else(|e| panic!("{e}"));
    assert!(slit_cells(&g, &c, &r, 2).is_empty());
}
