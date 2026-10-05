//! The coupled saddle-point solve of the KP20 bi-9 fair interior
//! (PK-601): `[[G + ρI, Cᵀ], [C, −εI]] [D; λ] = [r; 0]` over the free
//! coefficients × 3 coordinates, where `G` is block-diagonal by piece
//! (dense `300 × 300` blocks, coordinates coupled by the frozen normal)
//! and `C` is the homogeneous relation rows replicated per coordinate.
//!
//! The system is quasi-definite, so it is factored by an AMD-ordered
//! supernodal LDLᵀ without pivoting. The relation block is
//! rank-deficient (redundant collocated rows), so the regularisation
//! `ε` must be large enough to keep the factorisation stable, and
//! iterative refinement against the unregularised system removes its
//! bias. Measured (design doc, PK-601): `ε = 1e-12` with no refinement
//! gave a wrong answer on n = 3 (cube 2.95 against 1.82); `ε = 1e-8`
//! with refinement reproduces the pivoted LU to a 1e-16 residual.
//!
//! The pattern, its fill-reducing ordering and the symbolic factor
//! depend only on `(n, eq (10) coefficients)`: they are memoised, with
//! each generated entry's position in the pattern, so a corner scatters
//! its values and factors numerically.

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use faer::Mat;
use faer::sparse::SymbolicSparseColMat;
use popkern_math::linalg::LdltSymbolic;

use super::kp20_bi9_fair_curvature::PIECE_DIM;
use super::kp20_bi9_relations::RelationSpace;

/// Ridge on the unknown block, in the Gram's normalised scale.
const RIDGE: f64 = 1e-10;
/// Regularisation of the relation block (see the module docs).
const RELATION_REGULARISATION: f64 = 1e-8;
/// Refinement passes against the unregularised system.
const REFINEMENT_PASSES: usize = 3;

/// Why the coupled solve declined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum KktDecline {
    /// The pattern or its symbolic factor could not be built.
    Structure,
    /// A non-finite value in the solution.
    NonFinite,
    /// The refined residual stayed above its gate (relative to `‖r‖`).
    Residual(f64),
}

/// Where each generated entry's value comes from.
#[derive(Clone, Copy)]
enum Source {
    /// Piece `p`'s Gram, local `(a, b)`.
    Gram(u32, u16, u16),
    /// A relation row weight.
    Relation(f64),
    /// A diagonal of the unknown block (ridge only).
    UnknownDiagonal,
    /// A diagonal of the relation block (−ε only).
    RelationDiagonal,
}

/// The memoised structure for one `(n, b)`.
struct Pattern {
    size: usize,
    /// Lower-triangular CSC pattern.
    symbolic: SymbolicSparseColMat<usize>,
    /// Position of each generated entry in the value array.
    positions: Vec<usize>,
    sources: Vec<Source>,
    /// Diagonal positions of the unknowns and of the relation rows.
    unknown_diagonals: Vec<usize>,
    relation_diagonals: Vec<usize>,
    /// AMD order and `L D Lᵀ` analysis of the pattern.
    ldlt: LdltSymbolic,
}

type Key = (usize, [u64; 6]);

fn cache() -> &'static Mutex<HashMap<Key, Arc<Pattern>>> {
    static CACHE: OnceLock<Mutex<HashMap<Key, Arc<Pattern>>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Every lower-triangular entry in a fixed generation order.
fn entries(n: usize, space: &RelationSpace) -> Vec<(usize, usize, Source)> {
    let count = 3 * space.count;
    let unknown = |s: usize, k: usize, a: usize| {
        let x = space.index[s][k][a / 30][(a / 3) % 10];
        (x != usize::MAX).then(|| 3 * x + a % 3)
    };
    let mut out = Vec::new();
    for s in 0..n {
        for k in 0..4 {
            let piece = (4 * s + k) as u32;
            for a in 0..PIECE_DIM {
                let Some(xa) = unknown(s, k, a) else {
                    continue;
                };
                for b in a..PIECE_DIM {
                    let Some(xb) = unknown(s, k, b) else {
                        continue;
                    };
                    out.push((
                        xa.max(xb),
                        xa.min(xb),
                        Source::Gram(piece, a as u16, b as u16),
                    ));
                }
            }
        }
    }
    for x in 0..count {
        out.push((x, x, Source::UnknownDiagonal));
    }
    for (r, row) in space.rows.iter().enumerate() {
        for c in 0..3 {
            let rr = count + 3 * r + c;
            for &(x, w) in row {
                out.push((rr, 3 * x + c, Source::Relation(w)));
            }
            out.push((rr, rr, Source::RelationDiagonal));
        }
    }
    out
}

fn build_pattern(n: usize, space: &RelationSpace) -> Option<Pattern> {
    let size = 3 * space.count + 3 * space.rows.len();
    let generated = entries(n, space);
    let mut order: Vec<usize> = (0..generated.len()).collect();
    order.sort_unstable_by_key(|&t| (generated[t].1, generated[t].0));
    let mut col_ptr = vec![0usize; size + 1];
    let mut row_idx = Vec::with_capacity(generated.len());
    let mut positions = vec![0usize; generated.len()];
    let mut last: Option<(usize, usize)> = None;
    for &t in &order {
        let (row, col, _) = generated[t];
        if last != Some((row, col)) {
            row_idx.push(row);
            col_ptr[col + 1] += 1;
            last = Some((row, col));
        }
        positions[t] = row_idx.len() - 1;
    }
    for c in 0..size {
        col_ptr[c + 1] += col_ptr[c];
    }
    let ldlt = super::kp20_bi9_ldlt::analyse(size, &col_ptr, &row_idx)?;
    let symbolic = SymbolicSparseColMat::new_checked(size, size, col_ptr, None, row_idx);
    let mut unknown_diagonals = Vec::new();
    let mut relation_diagonals = Vec::new();
    for (t, &(_, _, source)) in generated.iter().enumerate() {
        match source {
            Source::UnknownDiagonal => unknown_diagonals.push(positions[t]),
            Source::RelationDiagonal => relation_diagonals.push(positions[t]),
            _ => {}
        }
    }
    Some(Pattern {
        size,
        symbolic,
        sources: generated.iter().map(|e| e.2).collect(),
        positions,
        unknown_diagonals,
        relation_diagonals,
        ldlt,
    })
}

fn pattern(n: usize, b: [[f64; 3]; 2], space: &RelationSpace) -> Option<Arc<Pattern>> {
    let key: Key = (
        n,
        [
            b[0][0].to_bits(),
            b[0][1].to_bits(),
            b[0][2].to_bits(),
            b[1][0].to_bits(),
            b[1][1].to_bits(),
            b[1][2].to_bits(),
        ],
    );
    if let Some(p) = cache()
        .lock()
        .unwrap_or_else(|poison| poison.into_inner())
        .get(&key)
    {
        return Some(p.clone());
    }
    let built = Arc::new(build_pattern(n, space)?);
    cache()
        .lock()
        .unwrap_or_else(|poison| poison.into_inner())
        .insert(key, built.clone());
    Some(built)
}

/// `y += A x` for the symmetric matrix stored as its lower triangle.
fn symmetric_product(pattern: &Pattern, values: &[f64], x: &[f64], y: &mut [f64]) {
    let col_ptr = pattern.symbolic.col_ptrs();
    let rows = pattern.symbolic.row_indices();
    for col in 0..pattern.size {
        for idx in col_ptr[col]..col_ptr[col + 1] {
            let (row, w) = (rows[idx], values[idx]);
            y[row] += w * x[col];
            if row != col {
                y[col] += w * x[row];
            }
        }
    }
}

/// Solve for the displacement of the free coefficients (`3 · count`,
/// ordered `3 · unknown + coordinate`). `grams[p]` is piece `p`'s
/// combined `300 × 300` Gram (`p = 4 · sector + piece`); `rhs` has one
/// entry per free coefficient × coordinate. Returns the displacement and
/// the refined residual relative to `‖rhs‖`.
pub(super) fn solve_coupled(
    n: usize,
    b: [[f64; 3]; 2],
    space: &RelationSpace,
    grams: &[Mat<f64>],
    rhs: &[f64],
) -> Result<(Vec<f64>, f64), KktDecline> {
    let pattern = pattern(n, b, space).ok_or(KktDecline::Structure)?;
    let mut exact = vec![0.0f64; pattern.symbolic.row_indices().len()];
    for (t, source) in pattern.sources.iter().enumerate() {
        let value = match *source {
            Source::Gram(p, a, b) => {
                let g = &grams[p as usize];
                let (a, b) = (a as usize, b as usize);
                // A diagonal pair is generated once; an off-diagonal
                // pair once for its lower entry.
                g.read(a, b)
            }
            Source::Relation(w) => w,
            Source::UnknownDiagonal | Source::RelationDiagonal => 0.0,
        };
        exact[pattern.positions[t]] += value;
    }
    let mut regularised = exact.clone();
    for &d in &pattern.unknown_diagonals {
        regularised[d] += RIDGE;
    }
    for &d in &pattern.relation_diagonals {
        regularised[d] -= RELATION_REGULARISATION;
    }
    let factor = pattern
        .ldlt
        .factor(&regularised)
        .ok_or(KktDecline::NonFinite)?;
    let size = pattern.size;
    let mut target = vec![0.0f64; size];
    target[..rhs.len()].copy_from_slice(rhs);
    let norm = target.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
    let (mut x, mut residual) = (vec![0.0f64; size], target.clone());
    let mut relative = 1.0;
    for _ in 0..REFINEMENT_PASSES {
        let mut step = residual.clone();
        pattern.ldlt.solve(&factor, &mut step);
        for (xi, si) in x.iter_mut().zip(&step) {
            *xi += si;
        }
        let mut product = vec![0.0f64; size];
        symmetric_product(&pattern, &exact, &x, &mut product);
        residual = (0..size).map(|i| target[i] - product[i]).collect();
        relative = residual.iter().map(|v| v * v).sum::<f64>().sqrt() / norm;
    }
    if !x.iter().all(|v| v.is_finite()) {
        return Err(KktDecline::NonFinite);
    }
    if relative.is_nan() || relative > 1e-9 {
        return Err(KktDecline::Residual(relative));
    }
    x.truncate(rhs.len());
    Ok((x, relative))
}
