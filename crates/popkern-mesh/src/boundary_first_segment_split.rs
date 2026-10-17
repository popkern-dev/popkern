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

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::{BoundaryRing, SplitRequest};

/// With an angular facet bound, a segment whose ends' surface normals
/// spread past this share of the bound splits at the lattice spacing
/// itself rather than 1.5x it (PK-687): beside a trim that cuts the
/// lattice obliquely, a 1.5-cell segment is the base of facets at the
/// full vertex-normal spread, which tilt past the bound (the coil
/// spring's ground-cut ends: four at 12.1-13.9 deg that no local split
/// could repair, every candidate made them worse).
const TURN_SPLIT_SHARE: f64 = 0.6;

/// Chart lattice description for the split rule: the chart-centre
/// metric scale `(su, sv)`, the per-direction metric cell sizes, the
/// square (single-count) cell size, and whether the lattice is
/// directional (per-direction counts differ from the square count).
pub(super) struct SegmentLattice<'a> {
    pub(super) scale: (f64, f64),
    pub(super) cells: (f64, f64),
    pub(super) square_cell: f64,
    pub(super) directional: bool,
    /// The chart and the angular facet bound (radians), when the preset
    /// bounds facet tilt and the chart is NURBS. Analytic charts already
    /// meet the tilt contract at 1.5 cells, and splitting a sphere
    /// corner's rims to one cell sent the boss plate's four corner
    /// patches to a fallback mesh with 1542 open edges each (filed as
    /// its own defect).
    pub(super) turn: Option<(&'a SurfaceGeometry, f64)>,
}

/// Request a midpoint split of every ring segment longer than 1.5x the
/// lattice spacing along its dominant direction (1x where its normals
/// turn, see [`TURN_SPLIT_SHARE`]). Segments whose
/// canonical span is empty (`f_lo == f_hi`: a chart-degenerate run
/// such as a pole row, whose samples are one 3D point) never request.
pub(super) fn request_long_segment_splits(
    rings: &[BoundaryRing],
    lattice: &SegmentLattice<'_>,
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
            let factor = match lattice.turn {
                Some((surface, bound)) if turns_past(surface, a, b, TURN_SPLIT_SHARE * bound) => {
                    1.0
                }
                _ => 1.5,
            };
            if length > factor * thr && f_lo != f_hi {
                split_requests.push((eh, 0.5 * (f_lo + f_hi)));
            }
        }
    }
}

/// Whether the surface normals at chart points `a` and `b` differ by
/// more than `angle` (radians). A chart-singular end never does.
fn turns_past(surface: &SurfaceGeometry, a: (f64, f64), b: (f64, f64), angle: f64) -> bool {
    let (Some(na), Some(nb)) = (
        surface.normal_at(a.0, a.1).normalized(),
        surface.normal_at(b.0, b.1).normalized(),
    ) else {
        return false;
    };
    na.dot(&nb) < popkern_core::geometry::PortableMath::pcos(angle)
}
