//! The triangle-to-cell facet clipper: one geometry for every surface
//! integral (cross-cutting Thread 1).
//!
//! Every triangle of a solid is clipped to each grid cell it crosses.
//! The resulting facet pieces carry their cell, face, triangle, area
//! and the triangle's unit normal. The moment rasteriser, loads,
//! supports and Robin terms all integrate over these same pieces.
//!
//! Cells are half-open for ties: a piece lying flat on a cell's upper
//! plane belongs to the next cell, so a facet on a grid plane is
//! counted once. Pieces are ordered by (cell, triangle).

pub mod polygon;

use crate::grid::Grid;
use crate::model::{FaceId, SolidMesh};

pub use polygon::{MAX_VERTICES, Polygon};

use polygon::{cross, sub};

/// One triangle's part inside one cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FacetPiece {
    /// Linear index of the cell.
    pub cell: usize,
    /// Index of the source triangle in the solid.
    pub triangle: u32,
    /// Face of the source triangle.
    pub face: FaceId,
    /// The clipped polygon, in the triangle's winding order.
    pub polygon: Polygon,
    /// Area of the polygon.
    pub area: f64,
    /// Unit outward normal of the source triangle.
    pub normal: [f64; 3],
}

/// All facet pieces of one solid.
#[derive(Debug, Clone, PartialEq)]
pub struct ClippedSolid {
    /// Pieces ordered by (cell, triangle).
    pub pieces: Vec<FacetPiece>,
    /// Triangles skipped because they have zero area.
    pub degenerate_triangles: usize,
}

/// Why a solid could not be clipped.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ClipError {
    /// A triangle reaches outside the grid.
    #[error("triangle {triangle} reaches outside the grid")]
    OutsideGrid {
        /// The triangle.
        triangle: u32,
    },
    /// A triangle index is out of range or the mesh is malformed.
    #[error("triangle {triangle} references a missing vertex or face tag")]
    Malformed {
        /// The triangle.
        triangle: u32,
    },
}

/// Inclusive cell range on one axis covering [lo, hi], widened by one
/// cell each side (clamped) so points on cell planes are never missed.
fn axis_range(grid: &Grid, axis: usize, lo: f64, hi: f64) -> Option<(usize, usize)> {
    let h = grid.cell_size();
    let o = grid.origin()[axis];
    let n = grid.dims()[axis];
    let first = ((lo - o) / h).floor();
    let last = ((hi - o) / h).floor();
    if !(first.is_finite() && last.is_finite()) || first < 0.0 || last >= n as f64 {
        return None;
    }
    let first = (first as usize).saturating_sub(1);
    let last = (last as usize + 1).min(n - 1);
    Some((first, last))
}

/// Clips `polygon` to the closed box of cell `ijk`, then applies the
/// half-open tie rule.
pub fn clip_to_cell(grid: &Grid, ijk: [usize; 3], polygon: &Polygon) -> Polygon {
    let lo = grid.cell_min(ijk);
    let hi = grid.cell_min([ijk[0] + 1, ijk[1] + 1, ijk[2] + 1]);
    let mut p = *polygon;
    for axis in 0..3 {
        p = p.clip(axis, lo[axis], true);
        if p.is_empty() {
            return p;
        }
        p = p.clip(axis, hi[axis], false);
        if p.is_empty() {
            return p;
        }
    }
    if (0..3).any(|axis| p.lies_on(axis, hi[axis])) {
        return Polygon::empty();
    }
    p
}

/// Clips every triangle of `solid` to the cells of `grid`.
pub fn clip_solid(grid: &Grid, solid: &SolidMesh) -> Result<ClippedSolid, ClipError> {
    let mut pieces = Vec::new();
    let mut degenerate_triangles = 0;
    for (t, corners) in solid.triangles.iter().enumerate() {
        let triangle = t as u32;
        let malformed = ClipError::Malformed { triangle };
        let face = *solid.triangle_face.get(t).ok_or(malformed.clone())?;
        let mut v = [[0.0; 3]; 3];
        for (slot, &i) in v.iter_mut().zip(corners) {
            *slot = *solid.vertices.get(i as usize).ok_or(malformed.clone())?;
        }
        let n = cross(sub(v[1], v[0]), sub(v[2], v[0]));
        let length = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if length == 0.0 {
            degenerate_triangles += 1;
            continue;
        }
        let normal = n.map(|c| c / length);
        let source = Polygon::triangle(v[0], v[1], v[2]);
        let mut range = [(0, 0); 3];
        for (axis, r) in range.iter_mut().enumerate() {
            let lo = v.iter().map(|p| p[axis]).fold(f64::INFINITY, f64::min);
            let hi = v.iter().map(|p| p[axis]).fold(f64::NEG_INFINITY, f64::max);
            *r = axis_range(grid, axis, lo, hi).ok_or(ClipError::OutsideGrid { triangle })?;
        }
        for k in range[2].0..=range[2].1 {
            for j in range[1].0..=range[1].1 {
                for i in range[0].0..=range[0].1 {
                    let polygon = clip_to_cell(grid, [i, j, k], &source);
                    if polygon.len() < 3 {
                        continue;
                    }
                    let area = polygon.area();
                    if area > 0.0 {
                        let cell = grid
                            .cell_index([i, j, k])
                            .ok_or(ClipError::OutsideGrid { triangle })?;
                        pieces.push(FacetPiece {
                            cell,
                            triangle,
                            face,
                            polygon,
                            area,
                            normal,
                        });
                    }
                }
            }
        }
    }
    // Stable: equal cells keep triangle order.
    pieces.sort_by_key(|p| p.cell);
    Ok(ClippedSolid {
        pieces,
        degenerate_triangles,
    })
}

#[cfg(test)]
#[path = "clip_tests.rs"]
mod clip_tests;
