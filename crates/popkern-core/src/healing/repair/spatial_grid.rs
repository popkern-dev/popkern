//! Spatial grid for accelerating proximity queries in boundary assembly.
//!
//! Partitions 3D space into uniform cells, allowing O(1) insertion
//! and O(27) neighbor-cell lookup. Reduces pairwise matching from
//! O(n^2) to O(n) for uniformly distributed points.

use std::collections::HashMap;

use crate::geometry::Point3;

/// Spatial grid for accelerating edge proximity queries.
pub(super) struct SpatialGrid {
    /// Map from integer cell coordinates to element indices.
    cells: HashMap<(i64, i64, i64), Vec<usize>>,
    /// Width of each cell in world-space units.
    cell_size: f64,
}

impl SpatialGrid {
    /// Create a new spatial grid with the given cell size.
    pub(super) fn new(cell_size: f64) -> Self {
        Self {
            cells: HashMap::new(),
            cell_size: if cell_size > 0.0 { cell_size } else { 1.0 },
        }
    }

    /// Compute the integer cell coordinates for a point.
    fn cell_coords(&self, point: &Point3) -> (i64, i64, i64) {
        let inv = 1.0 / self.cell_size;
        (
            (point.x * inv).floor() as i64,
            (point.y * inv).floor() as i64,
            (point.z * inv).floor() as i64,
        )
    }

    /// Insert an element index at the given point's cell.
    pub(super) fn insert(&mut self, index: usize, point: &Point3) {
        let coords = self.cell_coords(point);
        self.cells.entry(coords).or_default().push(index);
    }

    /// Query all element indices in the 3x3x3 neighborhood around a point.
    pub(super) fn query_nearby(&self, point: &Point3) -> Vec<usize> {
        let (cx, cy, cz) = self.cell_coords(point);
        let mut result = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let key = (cx + dx, cy + dy, cz + dz);
                    if let Some(indices) = self.cells.get(&key) {
                        result.extend_from_slice(indices);
                    }
                }
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_and_query_nearby() {
        let mut grid = SpatialGrid::new(1.0);
        let p0 = Point3::new(0.5, 0.5, 0.5);
        let p1 = Point3::new(0.8, 0.8, 0.8);
        let p_far = Point3::new(10.0, 10.0, 10.0);

        grid.insert(0, &p0);
        grid.insert(1, &p1);
        grid.insert(2, &p_far);

        // Query near p0 should find indices 0 and 1, not 2.
        let nearby = grid.query_nearby(&p0);
        assert!(nearby.contains(&0));
        assert!(nearby.contains(&1));
        assert!(!nearby.contains(&2));

        // Query near p_far should find index 2, not 0 or 1.
        let nearby_far = grid.query_nearby(&p_far);
        assert!(nearby_far.contains(&2));
        assert!(!nearby_far.contains(&0));
    }

    #[test]
    fn empty_grid_returns_empty() {
        let grid = SpatialGrid::new(1.0);
        let nearby = grid.query_nearby(&Point3::new(0.0, 0.0, 0.0));
        assert!(nearby.is_empty());
    }
}
