//! PH-054: exact fragment moments where every edge line of every cell
//! wall is reached by two or more fragments: eight octant blocks a slot
//! apart in all three axes, meeting inside one cell. Each wall of that
//! cell has four fragments and each of its edges two, so neither a
//! single-owner wall nor a single-owner edge line closes it.

use super::wall_split_tests::fragment_errors;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::moments::Splitting;

/// Slot width and the point the three slots cross: every slot face at
/// least 0.02 from a cell wall at all five grid shifts, since material
/// thinner than the lattice spacing beside a slot is a separate defect
/// (PH-060: at (1.09, 0.9125, 1.1625) a 0.005 sliver welds across).
const W: f64 = 0.05;
const CENTRE: [f64; 3] = [1.06, 0.95, 1.07];
const SIZE: f64 = 2.0;

/// The eight octants' (lower, upper) corners for slots crossing at
/// `centre`.
fn octants(centre: [f64; 3]) -> Vec<([f64; 3], [f64; 3])> {
    let d = 0.5 * W;
    (0..8)
        .map(|k: usize| {
            let lo = [0, 1, 2].map(|a| if k >> a & 1 == 0 { 0.0 } else { centre[a] + d });
            let hi = [0, 1, 2].map(|a| if k >> a & 1 == 0 { centre[a] - d } else { SIZE });
            (lo, hi)
        })
        .collect()
}

/// The octants as one solid (eight closed shells).
fn octant_model(centre: [f64; 3]) -> crate::model::AnalysisModel {
    let mut mesh = box_mesh([0.0; 3], [0.0; 3]);
    mesh.vertices.clear();
    mesh.triangles.clear();
    mesh.triangle_face.clear();
    let mut volume = 0.0;
    for (lo, hi) in octants(centre) {
        let b = box_mesh(lo, hi);
        let base = mesh.vertices.len() as u32;
        mesh.vertices.extend(b.vertices);
        mesh.triangles
            .extend(b.triangles.iter().map(|t| t.map(|i| i + base)));
        mesh.triangle_face.extend(b.triangle_face);
        volume += b.exact_volume.unwrap_or(0.0);
    }
    mesh.exact_volume = Some(volume);
    model_of(mesh, 6)
}

/// Measured 2026-10-11 over the five shifts: 1816, 420 and 32 fragments
/// of 2-, 4- and 8-fragment cells, exact to 5.6e-16 h³; before PH-054
/// the four 8-fragment cells kept lattice moments.
#[test]
fn the_octants_at_a_three_way_slot_crossing_get_exact_moments() {
    let model = octant_model(CENTRE);
    let forms = octants(CENTRE);
    let everywhere = ([-9.0; 3], [9.0; 3]);
    let hi = [SIZE; 3];
    for fragments in [2, 4, 8] {
        let (exact, checked, lattice) =
            fragment_errors(&model, hi, &forms, fragments, everywhere, Splitting::Exact);
        assert_eq!(lattice, 0, "{fragments} fragments");
        assert!(exact <= 1e-12, "{fragments} fragments: {exact}");
        assert!(checked > 0, "no cell of {fragments} fragments");
    }
    // Negative control: lattice moments miss the eight-fragment cells.
    let control = fragment_errors(&model, hi, &forms, 8, everywhere, Splitting::LatticeMoments);
    assert!(control.0 > 1e-6, "control {}", control.0);
}

/// PH-060: at this crossing and grid shift (0.13, 0.2, 0.09) the +x
/// octants' faces lie 0.005 (0.02 h) inside cells [5, 5, 1..3], under
/// the lattice spacing h/16. Before, those cells formed two fragments
/// for four octants (the slivers welded across the slot, moments off by
/// 9e-3 h³); measured now: every fragment exact to 1e-15 h³.
#[test]
fn slivers_beside_a_slot_are_their_own_fragments() {
    let centre = [1.09, 0.9125, 1.1625];
    let model = octant_model(centre);
    let forms = octants(centre);
    let everywhere = ([-9.0; 3], [9.0; 3]);
    let hi = [SIZE; 3];
    for fragments in [2, 4, 8] {
        let (exact, _, lattice) =
            fragment_errors(&model, hi, &forms, fragments, everywhere, Splitting::Exact);
        assert_eq!(lattice, 0, "{fragments} fragments");
        assert!(exact <= 1e-12, "{fragments} fragments: {exact}");
    }
    // The cells that welded hold four fragments now, two of them slivers.
    use crate::clip::clip_solid;
    use crate::grid::Grid;
    use crate::moments::rasterise_with;
    use crate::options::AnalysisOptions;
    let options = AnalysisOptions {
        target_cell_size: Some(super::wall_split_tests::H),
        grid_shift: [0.13, 0.2, 0.09],
        ..AnalysisOptions::default()
    };
    let grid = Grid::for_model([0.0; 3], hi, &options).unwrap_or_else(|e| panic!("{e}"));
    let mesh = &model.solids[0];
    let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
    let m = rasterise_with(&grid, mesh.id, &clipped, 2, &|_| false, Splitting::Exact)
        .unwrap_or_else(|e| panic!("{e}"));
    for k in 1..=3 {
        let cell = grid.cell_index([5, 5, k]).unwrap_or_else(|| panic!("cell"));
        let sc = m
            .fragments
            .split
            .iter()
            .find(|s| s.cell == cell)
            .unwrap_or_else(|| panic!("cell [5, 5, {k}] not split"));
        assert_eq!(m.entries(cell).len(), 4, "cell [5, 5, {k}]");
        assert_eq!(sc.slivers.len(), 2, "cell [5, 5, {k}]");
    }
}
