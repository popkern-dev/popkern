//! Symbolic eigenvector extraction over `Q[√d]`.
//!
//! Given a rank-2 mixed-sign rank-deficient quadric matrix `D` from
//! `find_rational_degenerate_fiber`, this module computes the two
//! real plane normals **exactly** in `Q[√d]`, where the radicand
//! `d = tr(D)² − 4·S₂(D)` arises from the 2×2 reduction of the
//! 4×4 characteristic polynomial:
//!
//! ```text
//! χ_D(μ) = μ²·(μ² − tr·μ + S₂)
//! ```
//!
//! (The `μ²` factor is the kernel — rank-2 ⇒ double zero eigenvalue.)
//! The non-zero eigenvalues `μ_± = (tr ± √d) / 2` live in `Q[√d]`.
//! Their eigenvectors `v_± ∈ Q[√d]⁴` come from Gauss-elim solving
//! `(D − μ_± I)·v = 0` in the field. The plane factorisation
//!
//! ```text
//! D = ½·((a + b)(a − b)ᵀ + (a − b)(a + b)ᵀ),  where
//!   a = √λ_+ · v_+,  b = √(-λ_-) · v_-
//! ```
//!
//! works fully in `Q[√d]` *iff* both `λ_+` and `−λ_-` are perfect
//! squares in `Q[√d]` (so we don't need a *nested* radical). For the
//! special case `tr = 0` (e.g., the Steinmetz pencil's λ=−1 fibre,
//! where `D = diag(−1, 0, 1, 0)`), `λ_± = ±√(d)/... ` simplifies — see
//! the test `steinmetz_planes_are_rational`. General cases that
//! require `Q[√d, √λ_+]` are out of scope here; the type system
//! catches the panic at `Mul`/`Add` between incompatible radicands.

#![allow(dead_code)]

use super::qsqrtd::QSqrtD;
use super::rational::{Mat4SymRat, Rat};
use num_traits::{One, Signed, Zero};

/// 4×4 symmetric matrix of `Q[√d]` values.
#[derive(Clone, Debug)]
pub(crate) struct Mat4SymQSqrtD(pub(crate) [[QSqrtD; 4]; 4]);

/// Lift a rational matrix into `Q[√d]` by embedding each entry as a
/// pure rational (`b = 0`). The radicand `d` is supplied so that
/// later operations (subtracting an `μ ∈ Q[√d]` from the diagonal,
/// row reductions involving multiplications by surds) stay in the
/// canonical form.
pub(crate) fn lift_to_qsqrtd(m: &Mat4SymRat) -> Mat4SymQSqrtD {
    let z = QSqrtD::zero();
    let mut out = [
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z],
    ];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = QSqrtD::from_rat(m.0[i][j].clone());
        }
    }
    Mat4SymQSqrtD(out)
}

/// `tr(D)` and `S₂(D)` (sum of 2×2 principal minors), both rational
/// for any rational matrix. `S₂ = ½·(tr² − tr(D²))` — but for our use
/// case it's cheaper to enumerate the six 2×2 principal minors.
pub(crate) fn trace_and_s2(d: &Mat4SymRat) -> (Rat, Rat) {
    let tr = &d.0[0][0] + &d.0[1][1] + &d.0[2][2] + &d.0[3][3];
    let mut s2 = Rat::zero();
    for i in 0..4 {
        for j in (i + 1)..4 {
            // 2×2 principal minor at rows/cols (i, j):
            //   d[i][i]·d[j][j] − d[i][j]·d[j][i]
            // (Symmetric matrix ⇒ d[i][j] = d[j][i].)
            s2 += &d.0[i][i] * &d.0[j][j] - &d.0[i][j] * &d.0[j][i];
        }
    }
    (tr, s2)
}

/// The two non-zero eigenvalues of a rank-2 mixed-sign rank-deficient
/// matrix `D` as `Q[√disc]` values, where `disc = tr² − 4·S₂`.
/// Returns `None` if `disc < 0` (eigenvalues are complex — `D` isn't
/// the rank-2 mixed-sign case after all).
pub(crate) fn rank2_eigenvalues(d: &Mat4SymRat) -> Option<(QSqrtD, QSqrtD, Rat)> {
    let (tr, s2) = trace_and_s2(d);
    let four = Rat::from(num_rational::BigRational::from_integer(
        num_bigint::BigInt::from(4),
    ));
    let disc = &tr * &tr - &four * &s2;
    if disc.is_negative() {
        return None;
    }
    // μ_± = (tr ± √disc) / 2.
    let two = Rat::from(num_rational::BigRational::from_integer(
        num_bigint::BigInt::from(2),
    ));
    let half_tr = &tr / &two;
    let half = Rat::one() / &two;
    let mu_plus = QSqrtD::new(half_tr.clone(), half.clone(), disc.clone());
    let mu_minus = QSqrtD::new(half_tr, -half, disc.clone());
    Some((mu_plus, mu_minus, disc))
}

/// Subtract `μ · I` from a `Mat4SymQSqrtD`. (Modifies a clone.)
pub(crate) fn subtract_mu_diag(m: &Mat4SymQSqrtD, mu: &QSqrtD) -> Mat4SymQSqrtD {
    let mut out = m.clone();
    for i in 0..4 {
        out.0[i][i] = out.0[i][i].clone() - mu.clone();
    }
    out
}

/// Find a non-zero null vector of a 4×4 matrix `m` over `Q[√d]` via
/// Gaussian elimination. Returns the first non-zero kernel basis
/// vector found, or `None` if the matrix has full rank.
///
/// Algorithm: row-reduce to RREF; identify the first free variable
/// (column without a pivot); set it to 1 and back-substitute for the
/// pivot variables. Pure field arithmetic — no tolerance.
pub(crate) fn null_vector_qsqrtd(m: &Mat4SymQSqrtD) -> Option<[QSqrtD; 4]> {
    // Copy rows.
    let mut rows: [[QSqrtD; 4]; 4] = m.0.clone();
    // Track which row contains the pivot for each column.
    let mut pivot_row_of_col: [Option<usize>; 4] = [None, None, None, None];
    let mut current_row = 0_usize;

    for col in 0..4 {
        if current_row >= 4 {
            break;
        }
        // Find a row at or below `current_row` with a non-zero entry
        // in this column.
        let pivot = (current_row..4).find(|&r| !rows[r][col].is_zero());
        let Some(pivot_row) = pivot else {
            continue; // free column
        };
        rows.swap(current_row, pivot_row);
        // Normalise the pivot row so the pivot entry = 1.
        let pivot_val = rows[current_row][col].clone();
        for slot in rows[current_row].iter_mut() {
            *slot = slot.clone() / pivot_val.clone();
        }
        // Eliminate this column in every other row.
        for r in 0..4 {
            if r != current_row && !rows[r][col].is_zero() {
                let factor = rows[r][col].clone();
                let pivot_row_clone = rows[current_row].clone();
                for (slot, p) in rows[r].iter_mut().zip(pivot_row_clone.iter()) {
                    let sub = factor.clone() * p.clone();
                    *slot = slot.clone() - sub;
                }
            }
        }
        pivot_row_of_col[col] = Some(current_row);
        current_row += 1;
    }

    // Find the first free column.
    let free_col = (0..4).find(|c| pivot_row_of_col[*c].is_none())?;
    // Build the null vector: free var = 1, others computed from RREF.
    let z = QSqrtD::zero();
    let mut v: [QSqrtD; 4] = [z.clone(), z.clone(), z.clone(), z];
    v[free_col] = QSqrtD::from_rat(Rat::one());
    for c in 0..4 {
        if c == free_col {
            continue;
        }
        if let Some(pr) = pivot_row_of_col[c] {
            // RREF row says: x_c + Σ_{free f} rows[pr][f]·x_f = 0.
            // We've set x_{free_col} = 1; other free vars = 0.
            // ⇒ x_c = −rows[pr][free_col].
            v[c] = -rows[pr][free_col].clone();
        }
    }
    Some(v)
}

/// Symbolically extract the two plane normals (as 4-vectors of
/// `Q[√d]` values) from a rank-2 mixed-sign quadric matrix `D`.
///
/// Returns `None` if:
///   - `disc = tr² − 4·S₂ < 0` (eigenvalues complex),
///   - `disc = 0` (rank-1 case — single doubled plane; caller can
///     handle separately),
///   - or the eigenvalue arithmetic exits Q[√d] (incompatible
///     radicand panics caught upstream).
pub(crate) fn extract_planes_qsqrtd(d: &Mat4SymRat) -> Option<[[QSqrtD; 4]; 2]> {
    let (mu_plus, mu_minus, disc) = rank2_eigenvalues(d)?;
    if disc.is_zero() {
        return None;
    }
    let lifted = lift_to_qsqrtd(d);
    let m_plus = subtract_mu_diag(&lifted, &mu_plus);
    let v_plus = null_vector_qsqrtd(&m_plus)?;
    let m_minus = subtract_mu_diag(&lifted, &mu_minus);
    let v_minus = null_vector_qsqrtd(&m_minus)?;
    // The plane normals are (proportional to) v_+ and v_- when both
    // λ_+ and −λ_- are perfect squares in Q[√d] — then a = √λ_+·v_+
    // collapses to v_+ up to a scalar. For the general case we'd
    // need Q[√d, √λ_+]. We return v_+ and v_- as "candidate planes"
    // — caller verifies via norm / f64 evaluation.
    Some([v_plus, v_minus])
}
