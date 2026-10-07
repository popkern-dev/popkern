//! Long-boundary-segment split requests for boundary-first charts
//! (split from `boundary_first.rs` under its line ratchet).
//!
//! LONG boundary segments on curved charts are split toward the
//! lattice spacing along them (requests feed the next tessellation
//! round, shared by both faces of the edge, so the mesh stays
//! watertight). Boundary samples graded like the interior keep the
//! triangles beside a seam well shaped: without this pass the display
//! mesh's seam normals roughened up to 5x (zebra attribution, cube
//! corner cap-cap slope 4.0e-3 -> 2.0e-2; PK-651). The threshold
//! follows the segment's DOMINANT direction on directional lattices:
//! judged by the much smaller u-cell, a v-running seam split to half
//! the row spacing and Bowyer–Watson fanned long-u chords between the
//! rows (the PK-146 class: ground-spring volume 0.0502 -> 0.0456).
//!
//! Segments are measured in the chart-centre metric applied to their
//! uv span, against a lattice that is uniform in parameter space, so
//! per direction the comparison is span against parameter cell. (A 3D
//! chord reads the same on regular charts but disagrees where the
//! centre metric is not local — on NURBS walls it refined boundaries
//! the lattice does not follow and deepened their in-face folds:
//! coil spring worst 31° -> 92°, propeller 85° -> 120°.) A pole run's
//! segments carry an empty span and never request (PK-626).

use super::{BoundaryRing, SplitRequest};

/// Chart lattice description for the split rule: the chart-centre
/// metric scale `(su, sv)`, the per-direction metric cell sizes, the
/// square (single-count) cell size, and whether the lattice is
/// directional (per-direction counts differ from the square count).
pub(super) struct SegmentLattice {
    pub(super) scale: (f64, f64),
    pub(super) cells: (f64, f64),
    pub(super) square_cell: f64,
    pub(super) directional: bool,
}

/// Request a midpoint split of every ring segment longer than 1.5x the
/// lattice spacing along its dominant direction. Segments whose
/// canonical span is empty (`f_lo == f_hi`: a chart-degenerate run
/// such as a pole row, whose samples are one 3D point) never request.
pub(super) fn request_long_segment_splits(
    rings: &[BoundaryRing],
    lattice: &SegmentLattice,
    split_requests: &mut Vec<SplitRequest>,
) {
    let (su, sv) = lattice.scale;
    for ring in rings {
        let m = ring.uv.len();
        for k in 0..m {
            let (a, b) = (ring.uv[k], ring.uv[(k + 1) % m]);
            let thr = if lattice.directional {
                if (su * (b.0 - a.0)).abs() > (sv * (b.1 - a.1)).abs() {
                    lattice.cells.0
                } else {
                    lattice.cells.1
                }
            } else {
                lattice.square_cell
            };
            let length = ((su * (b.0 - a.0)).powi(2) + (sv * (b.1 - a.1)).powi(2)).sqrt();
            let (eh, f_lo, f_hi) = ring.sources[k];
            if length > 1.5 * thr && f_lo != f_hi {
                split_requests.push((eh, 0.5 * (f_lo + f_hi)));
            }
        }
    }
}
