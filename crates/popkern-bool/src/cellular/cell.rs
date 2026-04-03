//! Cell and CellDecomposition types for volumetric boolean operations.
//!
//! The cellular approach partitions space into disjoint volumetric regions
//! (cells) bounded by the refined face fragments. Each cell knows which
//! input bodies contain it. Boolean operations reduce to selecting cells
//! by a predicate on their membership.

use popkern_core::geometry::Point3;

/// Unique identifier for a cell in the decomposition.
pub type CellId = usize;

/// A volumetric region in the cellular decomposition.
#[derive(Debug, Clone)]
pub struct Cell {
    /// Unique identifier.
    pub id: CellId,
    /// Which input bodies contain this cell.
    /// Bit 0 = body A, bit 1 = body B, etc.
    pub membership: u64,
    /// A representative interior point (for diagnostics and classification).
    pub sample_point: Point3,
}

impl Cell {
    /// Returns true if this cell is inside body A.
    pub fn in_a(&self) -> bool {
        self.membership & 1 != 0
    }

    /// Returns true if this cell is inside body B.
    pub fn in_b(&self) -> bool {
        self.membership & 2 != 0
    }
}

/// A face in the decomposition that separates two cells.
#[derive(Debug, Clone)]
pub struct CellBoundaryFace {
    /// Index of this face in the fragment list.
    pub fragment_index: usize,
    /// Cell on the positive-normal side of this face.
    pub cell_plus: CellId,
    /// Cell on the negative-normal side of this face.
    pub cell_minus: CellId,
}

/// Complete cellular decomposition of the space partitioned by two bodies.
///
/// The decomposition contains all cells (volumetric regions) and the
/// boundary faces that separate them. Boolean operations select cells
/// by their membership predicate.
#[derive(Debug, Clone)]
pub struct CellDecomposition {
    /// All cells in the decomposition.
    pub cells: Vec<Cell>,
    /// Boundary faces with references to cells on each side.
    pub boundaries: Vec<CellBoundaryFace>,
}

impl CellDecomposition {
    /// Number of cells.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Select cells matching a predicate on membership.
    ///
    /// Returns the indices of cells where `predicate(membership)` is true.
    pub fn select_cells(&self, predicate: impl Fn(u64) -> bool) -> Vec<CellId> {
        self.cells
            .iter()
            .filter(|c| predicate(c.membership))
            .map(|c| c.id)
            .collect()
    }

    /// Select cells for union: cells in A OR B.
    pub fn select_union(&self) -> Vec<CellId> {
        self.select_cells(|m| m & 1 != 0 || m & 2 != 0)
    }

    /// Select cells for intersection: cells in A AND B.
    pub fn select_intersection(&self) -> Vec<CellId> {
        self.select_cells(|m| m & 1 != 0 && m & 2 != 0)
    }

    /// Select cells for subtraction: cells in A AND NOT B.
    pub fn select_subtraction(&self) -> Vec<CellId> {
        self.select_cells(|m| m & 1 != 0 && m & 2 == 0)
    }

    /// Select cells for symmetric difference: cells in A XOR B.
    pub fn select_symmetric_difference(&self) -> Vec<CellId> {
        self.select_cells(|m| (m & 1 != 0) ^ (m & 2 != 0))
    }

    /// Get the boundary faces of a set of selected cells.
    ///
    /// A boundary face is included in the result if exactly one of its
    /// two adjacent cells is selected. Internal faces (both cells selected)
    /// and external faces (neither cell selected) are excluded.
    pub fn boundary_faces_of(&self, selected: &[CellId]) -> Vec<(usize, bool)> {
        let selected_set: std::collections::HashSet<CellId> = selected.iter().copied().collect();
        let mut result = Vec::new();

        for bf in &self.boundaries {
            let plus_in = selected_set.contains(&bf.cell_plus);
            let minus_in = selected_set.contains(&bf.cell_minus);

            if plus_in && !minus_in {
                // Face normal points outward (into non-selected region)
                result.push((bf.fragment_index, false));
            } else if minus_in && !plus_in {
                // Face normal points inward — flip it
                result.push((bf.fragment_index, true));
            }
            // Both or neither selected: internal/external face, skip
        }

        result
    }
}
