//! Per-solid assembly of the linear-elastic system: condensed EAS cells
//! from moments, Nitsche supports and traction loads on clipped facets,
//! reduced by the AgFEM extension applied per component.

use crate::aggregate::{Aggregation, Extension, Row};
use crate::clip::ClippedSolid;
use crate::element::operators::{CellOperators, Core};
use crate::element::vector::fragment_mean;
use crate::element::vector_boundary::{BoundaryTables, EnhancedFlux, PieceSupport};
use crate::element::{VectorTables, surface_moments};
use crate::grid::Grid;
use crate::linalg::sparse::{Csr, SparseError};
use crate::moments::SolidMoments;

#[path = "assemble_forces.rs"]
mod forces;
pub use forces::face_forces;

/// A condition on one face, in N and mm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FaceMechanics {
    /// Traction-free.
    Free,
    /// A support.
    Support(PieceSupport),
    /// A constant traction, N/mm².
    Traction([f64; 3]),
    /// A pressure, N/mm², acting along −n.
    Pressure(f64),
    /// Held (by Nitsche) to the solve's prescribed displacement field,
    /// affine per piece from the field's value and gradient at the
    /// piece's centroid: a submodel's cut faces, driven by the coarse
    /// solution.
    Prescribed,
    /// Loaded by the prescribed field's traction σ(∇u)·n at each
    /// piece's centroid, with a weak spring to its displacement against
    /// rigid motion: a force-driven submodel cap, which carries the
    /// coarse load path even where the coarse solution is too stiff (a
    /// singular corner).
    PrescribedTraction,
}

/// The coarse field a submodel's caps are driven by.
#[derive(Debug, Clone, Copy)]
pub struct Prescribed<'a> {
    /// The field.
    pub field: &'a super::DisplacementField,
    /// Stiffness of the weak spring of `PrescribedTraction` faces, N/mm³.
    pub spring: f64,
}

/// Material and Nitsche constants, MPa.
#[derive(Debug, Clone, Copy)]
pub struct Elastic {
    /// Lamé λ.
    pub lambda: f64,
    /// Lamé μ.
    pub mu: f64,
    /// Nitsche β = γ_D·(λ + 2μ)/h, N/mm³.
    pub beta: f64,
    /// Nitsche flux from the enhanced stress (production) rather than
    /// the plain Q1 stress (kept to show the latter goes indefinite
    /// near incompressibility, PH-034).
    pub enhanced_flux: bool,
    /// An infill core's Lamé parameters (`infill`).
    pub core: Option<(f64, f64)>,
}

/// The cell operators of a solve: the walls' material, and the core's
/// when there is one.
pub(super) fn cell_operators<'a>(
    vt: &'a VectorTables,
    h: f64,
    elastic: &Elastic,
) -> Result<CellOperators<'a>, crate::linalg::CholeskyError> {
    let operators = CellOperators::new(vt, h, elastic.lambda, elastic.mu)?;
    match elastic.core {
        Some(core) => operators.with_core(core),
        None => Ok(operators),
    }
}

/// Entry `entry`'s infill core.
pub(super) fn core_of<'a>(input: &'a SolidInput<'_>, entry: usize) -> &'a Core {
    static NONE: Core = Core::None;
    input.cores.and_then(|c| c.get(entry)).unwrap_or(&NONE)
}

/// Everything the assembly reads for one solid.
pub struct SolidInput<'a> {
    /// The grid.
    pub grid: &'a Grid,
    /// Clipped facet pieces.
    pub clipped: &'a ClippedSolid,
    /// Fragment moments.
    pub moments: &'a SolidMoments,
    /// Aggregation.
    pub aggregation: &'a Aggregation,
    /// Condition of each face, by `FaceId`.
    pub faces: &'a [FaceMechanics],
    /// The field `Prescribed` and `PrescribedTraction` faces are
    /// driven by.
    pub prescribed: Option<Prescribed<'a>>,
    /// A thermal strain, per entry (`thermoelastic`).
    pub thermal: Option<&'a super::thermoelastic::ThermalStrains>,
    /// Each entry's infill core (`infill`).
    pub cores: Option<&'a [Core]>,
    /// Each entry's stiffness scale (a generative design's densities
    /// through SIMP): the cell's block and its support terms scale, its
    /// applied loads do not.
    pub scales: Option<&'a [f64]>,
}

/// For a `PrescribedTraction` piece: the field's traction σ(∇u)·n, the
/// affine target (u0, G) of its spring and the spring stiffness, at the
/// piece's centroid.
#[allow(clippy::type_complexity)]
fn piece_traction(
    input: &SolidInput<'_>,
    elastic: &Elastic,
    piece: &crate::clip::FacetPiece,
) -> Option<([f64; 3], ([f64; 3], [[f64; 3]; 3]), f64)> {
    let prescribed = input.prescribed?;
    let v = piece.polygon.vertices();
    let nv = v.len().max(1) as f64;
    let c = [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv);
    let (u, g) = prescribed.field.value_and_gradient_at(c).ok()??;
    let strain = |i: usize, j: usize| 0.5 * (g[i][j] + g[j][i]);
    let trace = strain(0, 0) + strain(1, 1) + strain(2, 2);
    let n = piece.normal;
    let t = std::array::from_fn(|i| {
        (0..3)
            .map(|j| {
                let sigma = 2.0 * elastic.mu * strain(i, j)
                    + if i == j { elastic.lambda * trace } else { 0.0 };
                sigma * n[j]
            })
            .sum::<f64>()
    });
    let u0 = std::array::from_fn(|i| u[i] - (0..3).map(|j| g[i][j] * c[j]).sum::<f64>());
    Some((t, (u0, g), prescribed.spring))
}

/// The support a piece on a face with `condition` acts as: the face's
/// own, or for a prescribed face the prescribed field made affine at
/// the piece's centroid (`None` for any other condition, or where the
/// field has no value).
pub(super) fn piece_support(
    condition: FaceMechanics,
    input: &SolidInput<'_>,
    piece: &crate::clip::FacetPiece,
) -> Option<PieceSupport> {
    match condition {
        FaceMechanics::Support(support) => Some(support),
        FaceMechanics::Prescribed => {
            let v = piece.polygon.vertices();
            let nv = v.len().max(1) as f64;
            let c = [0, 1, 2].map(|a| v.iter().map(|x| x[a]).sum::<f64>() / nv);
            let (u, g) = input.prescribed?.field.value_and_gradient_at(c).ok()??;
            let u0 = std::array::from_fn(|i| u[i] - (0..3).map(|j| g[i][j] * c[j]).sum::<f64>());
            Some(PieceSupport::Fixed { u0, gradient: g })
        }
        _ => None,
    }
}

/// Why a solid's system could not be assembled.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AssemblyError {
    /// A sparse operation failed (an internal invariant).
    #[error(transparent)]
    Sparse(#[from] SparseError),
    /// A cell's enhanced block is not positive definite.
    #[error("cell {cell}: enhanced modes are not positive definite: {detail}")]
    Enhanced {
        /// Linear cell index.
        cell: usize,
        /// The factorisation's message.
        detail: String,
    },
}

/// One solid's reduced system and the parts kept for post-processing.
pub struct SolidSystem {
    /// Eᵀ A E over free dofs (3 per free node).
    pub a: Csr,
    /// Eᵀ f.
    pub b: Vec<f64>,
    /// Cells that skipped EAS (η below `ETA_EAS`).
    pub eas_disabled: usize,
}

/// E ⊗ I₃ as a sparse matrix (rows 3·node + i, columns 3·free + i).
pub fn vector_extension(e: &Extension) -> Result<Csr, SparseError> {
    let mut triplets = Vec::new();
    for (r, row) in e.rows.iter().enumerate() {
        for i in 0..3 {
            match row {
                Row::Free(d) => triplets.push((3 * r + i, 3 * *d as usize + i, 1.0)),
                Row::Extended(w) => {
                    triplets.extend(w.iter().map(|(d, c)| (3 * r + i, 3 * *d as usize + i, *c)));
                }
            }
        }
    }
    Csr::from_triplets(3 * e.nodes.len(), 3 * e.free, &triplets)
}

/// Positions of entry `entry`'s 8 nodes among the extension's active
/// nodes.
pub(super) fn cell_positions(
    grid: &Grid,
    moments: &SolidMoments,
    e: &Extension,
    entry: usize,
) -> Result<[usize; 8], SparseError> {
    let nodes = moments
        .node_keys(grid, entry)
        .ok_or(SparseError::Shape("cell outside grid"))?;
    let mut out = [0; 8];
    for (slot, n) in out.iter_mut().zip(nodes) {
        *slot = e
            .position(n)
            .ok_or(SparseError::Shape("cell node missing"))?;
    }
    Ok(out)
}

/// Support and load contributions of the pieces in one cell, by face
/// condition. Returns (support matrix, support load, applied load).
fn piece_terms(
    input: &SolidInput<'_>,
    entry: usize,
    op: &crate::element::vector::VectorCell,
    bt: &BoundaryTables,
    elastic: &Elastic,
) -> ([[f64; 24]; 24], [f64; 24], [f64; 24]) {
    let moments = input.moments;
    let enhanced = enhanced_flux(op, &moments.cells[entry], elastic);
    let cell = moments.cells[entry].cell;
    let grid = input.grid;
    let h = grid.cell_size();
    let lo = grid.cell_min(grid.cell_ijk(cell).unwrap_or_default());
    let pieces = &input.clipped.pieces;
    let start = pieces.partition_point(|p| p.cell < cell);
    let end = pieces.partition_point(|p| p.cell <= cell);
    let mut ks = [[0.0; 24]; 24];
    let mut fs = [0.0; 24];
    let mut fl = [0.0; 24];
    for (index, piece) in pieces[start..end].iter().enumerate() {
        if moments.entry_of_piece(start + index, cell) != Some(entry) {
            continue;
        }
        let condition = input
            .faces
            .get(piece.face.0 as usize)
            .copied()
            .unwrap_or(FaceMechanics::Free);
        let s = surface_moments(std::iter::once(piece), lo, h);
        match condition {
            FaceMechanics::Free => {}
            FaceMechanics::Support(_) | FaceMechanics::Prescribed => {
                let Some(support) = piece_support(condition, input, piece) else {
                    continue;
                };
                let (k, f) = bt.support(
                    support,
                    &s,
                    piece.normal,
                    lo,
                    h,
                    elastic.lambda,
                    elastic.mu,
                    elastic.beta,
                    enhanced,
                );
                for r in 0..24 {
                    fs[r] += f[r];
                    for c in 0..24 {
                        ks[r][c] += k[r][c];
                    }
                }
            }
            FaceMechanics::Traction(t) => {
                for (x, y) in fl.iter_mut().zip(bt.traction(&s, t)) {
                    *x += y;
                }
            }
            FaceMechanics::Pressure(p) => {
                let t = piece.normal.map(|n| -p * n);
                for (x, y) in fl.iter_mut().zip(bt.traction(&s, t)) {
                    *x += y;
                }
            }
            FaceMechanics::PrescribedTraction => {
                let Some((t, target, spring)) = piece_traction(input, elastic, piece) else {
                    continue;
                };
                for (x, y) in fl.iter_mut().zip(bt.traction(&s, t)) {
                    *x += y;
                }
                let (k, f) = bt.spring(&s, spring, target, lo, h);
                for r in 0..24 {
                    fs[r] += f[r];
                    for c in 0..24 {
                        ks[r][c] += k[r][c];
                    }
                }
            }
        }
    }
    (ks, fs, fl)
}

/// The cell's enhanced modes for the Nitsche flux, when it has them
/// and the flux uses them.
pub(super) fn enhanced_flux<'a>(
    op: &'a crate::element::vector::VectorCell,
    cell: &crate::moments::CellMoments,
    elastic: &Elastic,
) -> Option<EnhancedFlux<'a>> {
    if !elastic.enhanced_flux {
        return None;
    }
    op.recovery.as_ref().map(|recovery| EnhancedFlux {
        recovery,
        mean: fragment_mean(&cell.m),
    })
}

/// Entry `entry`'s unscaled block, material and supports (what a
/// stiffness scale multiplies), and its loads, the part a scale
/// multiplies and the applied part, with its Nitsche β from `betas`.
pub(super) fn unscaled_terms(
    operators: &CellOperators<'_>,
    bt: &BoundaryTables,
    input: &SolidInput<'_>,
    elastic: &Elastic,
    betas: &[f64],
    entry: usize,
) -> Result<UnscaledTerms, AssemblyError> {
    let cell = &input.moments.cells[entry];
    let operator = operators
        .of_core(
            cell,
            input.moments.is_fragment(entry),
            core_of(input, entry),
        )
        .map_err(|e| AssemblyError::Enhanced {
            cell: cell.cell,
            detail: e.to_string(),
        })?;
    let (op, _) = operator.get();
    let elastic = &Elastic {
        beta: betas.get(entry).copied().unwrap_or(elastic.beta),
        ..*elastic
    };
    let (ks, mut fs, fl) = piece_terms(input, entry, op, bt, elastic);
    if let Some(th) = super::thermoelastic::entry_load(bt, input, entry, op, elastic) {
        for (x, y) in fs.iter_mut().zip(th) {
            *x += y;
        }
    }
    let k = std::array::from_fn(|r| std::array::from_fn(|c| op.k[r][c] + ks[r][c]));
    Ok((k, fs, fl))
}

/// An entry's unscaled block, its scaled load and its applied load.
pub(super) type UnscaledTerms = ([[f64; 24]; 24], [f64; 24], [f64; 24]);

/// Entries assembled per parallel chunk (each block is 4.6 kB).
const ASSEMBLY_CHUNK: usize = 4096;

/// Entries in a chunk below which it is assembled on one thread.
const PARALLEL_MIN_CELLS: usize = 256;

/// One entry's condensed block and load, before the scatter.
struct EntryBlock {
    k: Box<[[f64; 24]; 24]>,
    load: [f64; 24],
    used: bool,
}

/// Assembles and reduces one solid's system.
pub fn assemble_solid(
    vt: &VectorTables,
    bt: &BoundaryTables,
    input: &SolidInput<'_>,
    elastic: &Elastic,
) -> Result<SolidSystem, AssemblyError> {
    let grid = input.grid;
    let h = grid.cell_size();
    let e = &input.aggregation.extension;
    let n = 3 * e.nodes.len();
    // Each cell couples its eight nodes: build that pattern once and add
    // the cell blocks into it in cell order (the order a triplet sort
    // would sum them in, so the bits match).
    let mut elements = Vec::with_capacity(input.moments.cells.len());
    for entry in 0..input.moments.cells.len() {
        elements.push(cell_positions(grid, input.moments, e, entry)?);
    }
    let mut a_nodes = Csr::block_pattern(e.nodes.len(), 3, &elements)?;
    let mut load = vec![0.0; n];
    let mut eas_disabled = 0;
    let operators = cell_operators(vt, h, elastic).map_err(|e| AssemblyError::Enhanced {
        cell: 0,
        detail: e.to_string(),
    })?;
    let betas = super::stability::entry_betas(&operators, bt, input, elastic)?;
    // Each entry's block and load, computed on threads a chunk at a
    // time and added in entry order exactly as one thread would.
    let block = |entry: usize| -> Result<EntryBlock, AssemblyError> {
        let cell = &input.moments.cells[entry];
        let operator = operators
            .of_core(
                cell,
                input.moments.is_fragment(entry),
                core_of(input, entry),
            )
            .map_err(|e| AssemblyError::Enhanced {
                cell: cell.cell,
                detail: e.to_string(),
            })?;
        let (op, used) = operator.get();
        let elastic = &Elastic {
            beta: betas.get(entry).copied().unwrap_or(elastic.beta),
            ..*elastic
        };
        let (ks, mut fs, fl) = piece_terms(input, entry, op, bt, elastic);
        if let Some(th) = super::thermoelastic::entry_load(bt, input, entry, op, elastic) {
            for (x, y) in fs.iter_mut().zip(th) {
                *x += y;
            }
        }
        // A scaled entry's material holds its supports with it (its
        // Nitsche terms are the material's); applied loads stay.
        let scale = input
            .scales
            .and_then(|s| s.get(entry))
            .copied()
            .unwrap_or(1.0);
        let mut k = Box::new([[0.0; 24]; 24]);
        for r in 0..24 {
            for c in 0..24 {
                k[r][c] = scale * (op.k[r][c] + ks[r][c]);
            }
        }
        Ok(EntryBlock {
            k,
            load: std::array::from_fn(|r| scale * fs[r] + fl[r]),
            used: *used,
        })
    };
    let entries = input.moments.cells.len();
    for start in (0..entries).step_by(ASSEMBLY_CHUNK) {
        let range = start..(start + ASSEMBLY_CHUNK).min(entries);
        let blocks = crate::parallel::map_range(range.clone(), PARALLEL_MIN_CELLS, block);
        for (entry, computed) in range.zip(blocks) {
            let computed = computed?;
            if !computed.used {
                eas_disabled += 1;
            }
            let nodes = elements[entry];
            for r in 0..24 {
                load[3 * nodes[r / 3] + r % 3] += computed.load[r];
            }
            for (a, &na) in nodes.iter().enumerate() {
                for (b, &nb) in nodes.iter().enumerate() {
                    let block = std::array::from_fn(|i| {
                        std::array::from_fn(|j| computed.k[3 * a + i][3 * b + j])
                    });
                    a_nodes.add_block3(3 * na, 3 * nb, &block)?;
                }
            }
        }
    }
    let ev = vector_extension(e)?;
    let a = ev.transpose().mul(&a_nodes.mul(&ev)?)?;
    let mut b = vec![0.0; ev.cols()];
    ev.transpose().mul_vec(&load, &mut b)?;
    Ok(SolidSystem { a, b, eas_disabled })
}
