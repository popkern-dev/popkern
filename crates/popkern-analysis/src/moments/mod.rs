//! The column-integral moment rasteriser (cross-cutting Thread 1).
//!
//! For each active cell of a solid it computes the 27 moments
//! m_abc = ∫_frag ξ^a η^b ζ^c dV (a, b, c ≤ 2, cell-local coordinates
//! in [−1, 1]³), exactly for the tessellated polyhedron. By the
//! divergence theorem with G_a(ξ) = (ξ^(a+1) − (−1)^(a+1)) / (a+1):
//!
//! m_abc = (h/2) [ Σ_pieces ∫ G_a(ξ) η^b ζ^c n_x dA + G_a(1) · P_bc ],
//!
//! where P_bc = Σ ∫ η^b ζ^c n_x dA over the pieces of the same column
//! beyond the cell: the solid's cross-section on the cell's upper x
//! wall. Walls at ξ = −1 and the y and z walls contribute nothing.
//! Columns are swept from high x to low x with a fixed-order suffix
//! sum, so full, cut and empty cells all come out of one pass.

pub mod fragments;
pub mod lattice;
pub mod triangle;

use crate::clip::{ClippedSolid, FacetPiece};
use crate::error::{AnalysisError, Degeneracy};
use crate::grid::Grid;
use crate::model::SolidId;

use triangle::monomial_integrals;

/// Index of m_abc in a moment array.
pub const fn moment_index(a: usize, b: usize, c: usize) -> usize {
    a + 3 * b + 9 * c
}

/// How a cell meets a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellKind {
    /// No surface piece; the cell is inside the solid.
    Full,
    /// At least one surface piece crosses the cell.
    Cut,
}

/// The moments of one solid in one active cell.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellMoments {
    /// Linear cell index.
    pub cell: usize,
    /// Full or cut.
    pub kind: CellKind,
    /// m_abc at `moment_index(a, b, c)`, in physical volume units.
    pub m: [f64; 27],
    /// Fragment index within the cell (0 unless the cell is split; see
    /// [`fragments`]).
    pub fragment: u8,
}

impl CellMoments {
    /// Volume fraction η = m₀₀₀ / h³.
    pub fn volume_fraction(&self, h: f64) -> f64 {
        self.m[0] / (h * h * h)
    }
}

/// The moments of one solid over the grid.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidMoments {
    /// The solid.
    pub solid: SolidId,
    /// Active cell entries in (cell, fragment) order: one per cell, or
    /// one per fragment of a split cell.
    pub cells: Vec<CellMoments>,
    /// Σ m₀₀₀: the rasterised volume.
    pub volume: f64,
    /// Number of full cells.
    pub full: usize,
    /// Number of cut cells.
    pub cut: usize,
    /// Cells whose exact moments were inconsistent (fragment volume
    /// outside [0, h³] beyond round-off) and were replaced by the
    /// sub-cell lattice fallback.
    pub fallback_cells: usize,
    /// Split cells, piece fragments and dof copies.
    pub fragments: fragments::Fragments,
}

/// How a cell holding disconnected parts of one solid is handled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Splitting {
    /// One entry per cell: the parts weld (CC1's negative control).
    Welded,
    /// Fragments keep lattice moments plus a residual share (the
    /// containment exact moments replace; LF1's negative control).
    #[cfg(test)]
    LatticeMoments,
    /// Fragments get exact moments where a wall allows (production).
    Exact,
}

/// Relative slack on 0 ≤ m₀₀₀ ≤ h³ before a cell's exact moments are
/// treated as inconsistent.
pub const CONSISTENCY_TOL: f64 = 1.0e-9;

/// True when a cell's fragment volume is physically possible.
pub fn consistent(m: &[f64; 27], h: f64) -> bool {
    let cube = h * h * h;
    m[0].is_finite() && m[0] >= -CONSISTENCY_TOL * cube && m[0] <= (1.0 + CONSISTENCY_TOL) * cube
}

/// G_a(1) = (1 − (−1)^(a+1)) / (a+1): 2, 0, 2/3.
const G_AT_ONE: [f64; 3] = [2.0, 0.0, 2.0 / 3.0];

/// Maps a physical point to the local coordinates of the cell whose
/// minimum corner is `lo`.
fn local(p: [f64; 3], lo: [f64; 3], h: f64) -> [f64; 3] {
    [0, 1, 2].map(|a| 2.0 * (p[a] - lo[a]) / h - 1.0)
}

/// Adds one piece's contributions: `q[abc]` += ∫ G_a(ξ) η^b ζ^c n_x dA
/// and `r[bc]` += ∫ η^b ζ^c n_x dA.
fn add_piece(piece: &FacetPiece, lo: [f64; 3], h: f64, q: &mut [f64; 27], r: &mut [f64; 9]) {
    let nx = piece.normal[0];
    if nx == 0.0 {
        return;
    }
    let v = piece.polygon.vertices();
    for i in 1..v.len() - 1 {
        let tri = [v[0], v[i], v[i + 1]];
        let u = [0, 1, 2].map(|k| tri[1][k] - tri[0][k]);
        let w = [0, 1, 2].map(|k| tri[2][k] - tri[0][k]);
        let x = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        let area = 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        if area == 0.0 {
            continue;
        }
        let m = monomial_integrals(tri.map(|p| local(p, lo, h)), area);
        for c in 0..3 {
            for b in 0..3 {
                r[b + 3 * c] += nx * m[0][b][c];
                for a in 0..3 {
                    // ∫ G_a = (∫ξ^(a+1) − (−1)^(a+1) ∫1) / (a+1)
                    let sign = if a % 2 == 0 { -1.0 } else { 1.0 };
                    let g = (m[a + 1][b][c] - sign * m[0][b][c]) / (a + 1) as f64;
                    q[moment_index(a, b, c)] += nx * g;
                }
            }
        }
    }
}

/// Rasterises one solid from its clipped pieces. A cut cell whose exact
/// moments are inconsistent is recomputed on the 2^`subcell_depth`
/// lattice and counted. Declines an inside-out solid (non-positive
/// rasterised volume).
pub fn rasterise_solid(
    grid: &Grid,
    solid: SolidId,
    clipped: &ClippedSolid,
    subcell_depth: u8,
) -> Result<SolidMoments, AnalysisError> {
    rasterise_with(
        grid,
        solid,
        clipped,
        subcell_depth,
        &|_| false,
        Splitting::Exact,
    )
}

/// [`rasterise_solid`] with `force(cell)` sending chosen cut cells to
/// the fallback, so tests can exercise it on exact geometry, and an
/// explicit `split` mode (negative controls).
pub(crate) fn rasterise_with(
    grid: &Grid,
    solid: SolidId,
    clipped: &ClippedSolid,
    subcell_depth: u8,
    force: &dyn Fn(usize) -> bool,
    split: Splitting,
) -> Result<SolidMoments, AnalysisError> {
    let h = grid.cell_size();
    let [nx, ny, nz] = grid.dims();
    let pieces = &clipped.pieces;
    let mut cursor = pieces.len();
    let mut cells = Vec::new();
    let mut fallback_cells = 0;
    // Sweep columns from the last to the first so pieces (sorted by
    // cell) are consumed from the end; within a column, high x first.
    for k in (0..nz).rev() {
        for j in (0..ny).rev() {
            let mut suffix = [0.0; 9];
            for i in (0..nx).rev() {
                let ijk = [i, j, k];
                let Some(cell) = grid.cell_index(ijk) else {
                    continue;
                };
                let start = pieces[..cursor].partition_point(|p| p.cell < cell);
                let mine = &pieces[start..cursor];
                cursor = start;
                let lo = grid.cell_min(ijk);
                let mut q = [0.0; 27];
                let mut r = [0.0; 9];
                for piece in mine {
                    add_piece(piece, lo, h, &mut q, &mut r);
                }
                let mut m = [0.0; 27];
                for c in 0..3 {
                    for b in 0..3 {
                        for (a, g) in G_AT_ONE.iter().enumerate() {
                            let idx = moment_index(a, b, c);
                            m[idx] = 0.5 * h * (q[idx] + g * suffix[b + 3 * c]);
                        }
                    }
                }
                for (s, add) in suffix.iter_mut().zip(r) {
                    *s += add;
                }
                if !mine.is_empty() && (!consistent(&m, h) || force(cell)) {
                    m = lattice::lattice_moments(grid, clipped, ijk, subcell_depth);
                    fallback_cells += 1;
                }
                let kind = if !mine.is_empty() {
                    Some(CellKind::Cut)
                } else if m[0] > 0.5 * h * h * h {
                    Some(CellKind::Full)
                } else {
                    None
                };
                if let Some(kind) = kind {
                    cells.push(CellMoments {
                        cell,
                        kind,
                        m,
                        fragment: 0,
                    });
                }
            }
        }
    }
    cells.reverse();
    let (cells, fragments) = if split != Splitting::Welded {
        let exact = split == Splitting::Exact;
        fragments::fragment(grid, clipped, cells, subcell_depth, exact).map_err(|e| match e {
            fragments::FragmentError::TooMany(index) => AnalysisError::TooManyFragments {
                index,
                limit: fragments::NODE_COPIES,
                cell: grid.cell_size(),
            },
        })?
    } else {
        let whole = fragments::Fragments {
            depth: subcell_depth,
            piece_fragment: vec![0; clipped.pieces.len()],
            corner_copy: vec![[0; 8]; cells.len()],
            ..fragments::Fragments::default()
        };
        (cells, whole)
    };
    let volume: f64 = cells.iter().map(|c| c.m[0]).sum();
    if volume <= 0.0 {
        return Err(AnalysisError::DegenerateGeometry {
            solid: Some(solid),
            detail: Degeneracy::InsideOut,
        });
    }
    let full = cells.iter().filter(|c| c.kind == CellKind::Full).count();
    Ok(SolidMoments {
        solid,
        cut: cells.len() - full,
        full,
        cells,
        volume,
        fallback_cells,
        fragments,
    })
}

/// Cut cells of `moments` whose inside lattice points form more than
/// one 6-connected component: parts of the solid closer than the cell
/// size, joined by the cell's shared dofs.
pub fn slit_cells(
    grid: &Grid,
    clipped: &ClippedSolid,
    moments: &SolidMoments,
    subcell_depth: u8,
) -> Vec<usize> {
    let mut out: Vec<usize> = moments
        .cells
        .iter()
        .filter(|c| c.kind == CellKind::Cut)
        .filter_map(|c| {
            let ijk = grid.cell_ijk(c.cell)?;
            (lattice::fragment_components(grid, clipped, ijk, subcell_depth) > 1).then_some(c.cell)
        })
        .collect();
    out.dedup();
    out
}

/// Gate X5 (exact half): the rasterised volume equals the mesh's
/// enclosed volume within `rel_tol`. A non-positive volume is
/// `InsideOut`; a disagreement is `VolumeMismatch`.
pub fn volume_parity(
    moments: &SolidMoments,
    mesh_volume: f64,
    rel_tol: f64,
) -> Result<(), AnalysisError> {
    let degenerate = |detail| AnalysisError::DegenerateGeometry {
        solid: Some(moments.solid),
        detail,
    };
    if mesh_volume <= 0.0 || moments.volume <= 0.0 {
        return Err(degenerate(Degeneracy::InsideOut));
    }
    if (moments.volume - mesh_volume).abs() > rel_tol * mesh_volume {
        return Err(degenerate(Degeneracy::VolumeMismatch));
    }
    Ok(())
}

#[cfg(test)]
#[path = "lattice_tests.rs"]
mod lattice_tests;
#[cfg(test)]
#[path = "moments_tests.rs"]
mod moments_tests;
#[cfg(test)]
#[path = "sliver_tests.rs"]
mod sliver_tests;
