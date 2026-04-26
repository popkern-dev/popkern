//! The quadric pencil `λA + μB`, its characteristic polynomial
//! `det(λA + μB)`, and the search for a degenerate (rank ≤ 3)
//! fiber.
//!
//! Levin's algorithm rests on the fact that the *intersection* of
//! two quadrics `A` and `B` equals the intersection of either one
//! with any other quadric in the pencil (every point on both `A`
//! and `B` satisfies `λA + μB = 0` for all `(λ, μ)`). So we are
//! free to swap one of the originals for a more tractable surface
//! in the pencil — specifically, a *ruled* (rank ≤ 3) quadric:
//! cone, plane pair, doubled plane, or pair of lines.
//!
//! Because each entry of `λA + μB` is linear in `(λ, μ)`, the
//! determinant `det(λA + μB)` is a homogeneous polynomial of total
//! degree 4 in `(λ, μ)`. Its roots `(λᵢ, μᵢ)` are the rank-deficient
//! fibers of the pencil. Generically there are 4 in projective
//! `(λ:μ)`-space; for real-symmetric pencils arising from real
//! quadric surfaces, at least one real root almost always exists
//! (Wang-Goldman-Tu §3 catalogues the failure modes).

use super::quadric::Mat4Sym;

const EPS_DET: f64 = 1e-9;

/// Construct `λ·A + μ·B`.
pub(crate) fn linear_comb(a: &Mat4Sym, b: &Mat4Sym, lambda: f64, mu: f64) -> Mat4Sym {
    let mut out = [[0.0; 4]; 4];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = lambda * a.0[i][j] + mu * b.0[i][j];
        }
    }
    Mat4Sym(out)
}

/// Determinant of a 4×4 matrix via cofactor expansion along row 0.
pub(crate) fn det4(m: &[[f64; 4]; 4]) -> f64 {
    let mut det = 0.0;
    for j in 0..4 {
        let cofactor = minor3x3_det(m, 0, j);
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        det += sign * m[0][j] * cofactor;
    }
    det
}

fn minor3x3_det(m: &[[f64; 4]; 4], skip_row: usize, skip_col: usize) -> f64 {
    let mut sub = [[0.0; 3]; 3];
    let mut ri = 0;
    for (i, row) in m.iter().enumerate() {
        if i == skip_row {
            continue;
        }
        let mut cj = 0;
        for (j, &val) in row.iter().enumerate() {
            if j == skip_col {
                continue;
            }
            sub[ri][cj] = val;
            cj += 1;
        }
        ri += 1;
    }
    sub[0][0] * (sub[1][1] * sub[2][2] - sub[1][2] * sub[2][1])
        - sub[0][1] * (sub[1][0] * sub[2][2] - sub[1][2] * sub[2][0])
        + sub[0][2] * (sub[1][0] * sub[2][1] - sub[1][1] * sub[2][0])
}

/// Coefficients of `det(λA + B)` as a univariate polynomial of degree
/// 4 in `λ`, returned in **monomial-descending** order:
/// `c[0]·λ⁴ + c[1]·λ³ + c[2]·λ² + c[3]·λ + c[4]`.
///
/// Computed via Lagrange interpolation at five λ values
/// (0, ±1, ±2). The closed-form inverse-Vandermonde for these
/// specific points lets us avoid a general 5×5 solve.
pub(crate) fn char_poly_in_lambda(a: &Mat4Sym, b: &Mat4Sym) -> [f64; 5] {
    // Evaluate det(λA + B) at λ ∈ {0, 1, -1, 2, -2}.
    let lambdas = [0.0_f64, 1.0, -1.0, 2.0, -2.0];
    let mut s = [0.0_f64; 5];
    for k in 0..5 {
        let m = linear_comb(a, b, lambdas[k], 1.0);
        s[k] = det4(&m.0);
    }
    // s[0] = c4
    // s[1] = c0 + c1 + c2 + c3 + c4
    // s[2] = c0 - c1 + c2 - c3 + c4
    // s[3] = 16c0 + 8c1 + 4c2 + 2c3 + c4
    // s[4] = 16c0 - 8c1 + 4c2 - 2c3 + c4
    let c4 = s[0];
    let aa = s[1] + s[2] - 2.0 * c4; // 2c0 + 2c2
    let bb = s[3] + s[4] - 2.0 * c4; // 32c0 + 8c2
    let c0 = (bb - 4.0 * aa) / 24.0;
    let c2 = aa / 2.0 - c0;
    let cc = s[1] - s[2]; // 2c1 + 2c3
    let dd = s[3] - s[4]; // 16c1 + 4c3
    let c1 = (dd - 2.0 * cc) / 12.0;
    let c3 = cc / 2.0 - c1;
    [c0, c1, c2, c3, c4]
}

/// A degenerate quadric in the pencil, found by zeroing the
/// characteristic polynomial.
#[derive(Copy, Clone, Debug)]
pub(crate) struct DegenerateFiber {
    /// The rank-deficient quadric `λA + μB`.
    pub(crate) quadric: Mat4Sym,
    /// Pencil coordinate `λ`.
    pub(crate) lambda: f64,
    /// Pencil coordinate `μ`.
    pub(crate) mu: f64,
}

/// Search the pencil for a real degenerate fiber. Tries:
///   1. Real roots of `det(λA + B) = 0` in finite λ.
///   2. The fiber `(λ, μ) = (1, 0)` itself (i.e., `A`).
///   3. The fiber `(λ, μ) = (0, 1)` itself (i.e., `B`).
///
/// Returns `None` if no real degenerate fiber exists in the pencil.
pub(crate) fn find_degenerate_fiber(a: &Mat4Sym, b: &Mat4Sym) -> Option<DegenerateFiber> {
    let coefs = char_poly_in_lambda(a, b);
    let lambdas = super::polynomial::solve_quartic(coefs);
    let scale = matrix_norm(&a.0).max(matrix_norm(&b.0)).max(1.0);
    let det_tol = EPS_DET * scale.powi(4);
    for lam in &lambdas {
        let m = linear_comb(a, b, *lam, 1.0);
        if det4(&m.0).abs() < det_tol {
            return Some(DegenerateFiber {
                quadric: m,
                lambda: *lam,
                mu: 1.0,
            });
        }
    }
    if det4(&a.0).abs() < det_tol {
        return Some(DegenerateFiber {
            quadric: *a,
            lambda: 1.0,
            mu: 0.0,
        });
    }
    if det4(&b.0).abs() < det_tol {
        return Some(DegenerateFiber {
            quadric: *b,
            lambda: 0.0,
            mu: 1.0,
        });
    }
    None
}

fn matrix_norm(m: &[[f64; 4]; 4]) -> f64 {
    let mut n = 0.0_f64;
    for row in m.iter() {
        for &v in row {
            n = n.max(v.abs());
        }
    }
    n
}
