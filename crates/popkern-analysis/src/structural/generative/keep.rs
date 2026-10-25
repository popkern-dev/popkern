//! Material a design must keep, near faces the caller names (mounting
//! faces, bolt holes), as a density floor per entry: the fraction of
//! the entry's volume within the keep thickness of such a face, from an
//! 8³ lattice. A floor rather than a kept/free flag, so the kept volume
//! moves smoothly with the grid: keeping whole entries by their centre
//! moved a thin cantilever's compliance by up to 46% across five grid
//! offsets (measured 2026-10-13).

use super::super::solve::Prepared;
use super::filter::DesignEntry;

/// Lattice points per cell side.
const LATTICE: usize = 8;

/// Each entry's density floor in [0, 1], in (solid, entry) order.
pub(super) fn keep_floors(
    prepared: &Prepared,
    faces: &[crate::model::FaceId],
    entries: &[DesignEntry],
    thickness: f64,
) -> Vec<f64> {
    let mut floors = vec![0.0; entries.len()];
    if faces.is_empty() {
        return floors;
    }
    let grid = &prepared.grid;
    let h = grid.cell_size();
    let reach = (thickness / h).ceil().max(0.0) as i64 + 1;
    let mut first = 0;
    for (_, clipped, moments, _) in &prepared.solids {
        let pieces: Vec<&crate::clip::FacetPiece> = clipped
            .pieces
            .iter()
            .filter(|p| faces.contains(&p.face))
            .collect();
        for (e, cell) in moments.cells.iter().enumerate() {
            let Some(ijk) = grid.cell_ijk(cell.cell) else {
                continue;
            };
            // The kept pieces within reach of the cell.
            let near: Vec<&&crate::clip::FacetPiece> = pieces
                .iter()
                .filter(|p| {
                    grid.cell_ijk(p.cell).is_some_and(|q| {
                        (0..3).all(|a| (q[a] as i64 - ijk[a] as i64).abs() <= reach)
                    })
                })
                .collect();
            if near.is_empty() {
                continue;
            }
            let distance = |x: [f64; 3]| {
                near.iter()
                    .map(|p| {
                        super::super::infill::polygon_distance(x, p.polygon.vertices(), p.normal)
                    })
                    .fold(f64::INFINITY, f64::min)
            };
            let lo = grid.cell_min(ijk);
            let (mut inside, mut kept) = (0usize, 0usize);
            for k in 0..LATTICE {
                for j in 0..LATTICE {
                    for i in 0..LATTICE {
                        let x = [i, j, k].map(|q| (q as f64 + 0.5) / LATTICE as f64);
                        let p = [0, 1, 2].map(|a| lo[a] + h * x[a]);
                        if !crate::moments::lattice::contains(grid, clipped, ijk, p) {
                            continue;
                        }
                        // A split cell's fragments share its lattice; the
                        // floor is the cell's.
                        inside += 1;
                        kept += usize::from(distance(p) <= thickness);
                    }
                }
            }
            if inside > 0
                && let Some(slot) = floors.get_mut(first + e)
            {
                *slot = kept as f64 / inside as f64;
            }
        }
        first += moments.cells.len();
    }
    floors
}
