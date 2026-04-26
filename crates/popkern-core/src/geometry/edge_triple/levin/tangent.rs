//! Stage 2b.40: symbolic tangent-intersection detection.
//!
//! Two quadric surfaces are *tangent along a curve* iff their pencil
//! `λA + μB` contains a **rank-1** fibre — i.e., a quadric matrix
//! `D = c · n · nᵀ` for some scalar `c` and 4-vector `n`. Geometrically,
//! `pᵀ·D·p = c·(nᵀp)²` vanishes precisely on the plane `nᵀp = 0`
//! (with multiplicity 2), so the intersection of the two original
//! quadrics lies entirely in that plane and is tangent to it
//! (multiplicity-2 contact along the intersection curve).
//!
//! ## Why symbolic
//!
//! The f64 Levin path treats every degenerate fibre as rank-2 and
//! computes two plane normals from a 2×2 eigendecomposition. When the
//! true rank is 1, the two computed normals are *coincident in
//! arithmetic* but differ by floating-point noise — the algorithm
//! produces two slightly-different planes for what should be a single
//! doubled plane, and the downstream branch-pair stitching chokes.
//!
//! The symbolic test sidesteps this entirely: in `Q`, "every 2×2
//! minor is zero" is a clean exact predicate. When it holds, the
//! plane normal is recovered from any non-zero column of `D` (which
//! is `c·n_j·n` — parallel to `n`).
//!
//! ## CAD relevance
//!
//! Cylinder-with-inscribed-sphere (tangent along the equator) is the
//! canonical case. Two cones sharing apex and axis but with the same
//! half-angle — coincident, handled by `surfaces_are_coincident` —
//! are NOT a tangent case. Two cones sharing apex and axis with
//! distinct half-angles intersect at the apex only (0-D, out of
//! scope). Sphere-tangent-to-cylinder externally is a tangent-curve
//! case; sphere-kissing-sphere is *not* (0-D contact, no curve).

#![allow(dead_code)] // helper module; dispatch wiring lands in 2b.41+.

use num_traits::{One, Signed, Zero};

use super::rational::{Mat4SymRat, Rat, char_poly_rat, linear_comb_rat, rational_roots};

/// True iff every 2×2 minor of `d` vanishes — equivalently, `rank(d) ≤ 1`.
/// Symbolic, exact: no tolerance, no FP noise.
pub(crate) fn is_rank1_or_less(d: &Mat4SymRat) -> bool {
    for i in 0..4 {
        for j in (i + 1)..4 {
            for k in 0..4 {
                for l in (k + 1)..4 {
                    let minor = &d.0[i][k] * &d.0[j][l] - &d.0[i][l] * &d.0[j][k];
                    if !minor.is_zero() {
                        return false;
                    }
                }
            }
        }
    }
    true
}

/// Recover the plane normal `n` from a rank-1 quadric `D = c · n · nᵀ`.
///
/// Strategy: pick the diagonal index `j` with the largest `|d_jj|`.
/// Since `d_ij = c · n_i · n_j`, column `j` reads `c · n_j · n` — i.e.,
/// it's parallel to `n` (with `n_j ≠ 0` because we picked the largest
/// non-zero diagonal). Returning that column gives `n` up to overall
/// scaling, which is exactly what plane-normal callers need.
///
/// Returns `None` only when `D` is the zero matrix (no plane).
pub(crate) fn extract_plane_from_rank1(d: &Mat4SymRat) -> Option<[Rat; 4]> {
    let mut best_j: Option<usize> = None;
    let mut best_abs: Option<Rat> = None;
    for j in 0..4 {
        if !d.0[j][j].is_zero() {
            let abs_jj = d.0[j][j].abs();
            let pick = match &best_abs {
                None => true,
                Some(prev) => &abs_jj > prev,
            };
            if pick {
                best_j = Some(j);
                best_abs = Some(abs_jj);
            }
        }
    }
    let j = best_j?;
    Some([
        d.0[0][j].clone(),
        d.0[1][j].clone(),
        d.0[2][j].clone(),
        d.0[3][j].clone(),
    ])
}

/// Find a tangent plane in the pencil `λA + B`, if one exists.
///
/// Walks every rational root `λ*` of `det(λA + B) = 0`, builds the
/// fibre `M = λ*·A + B`, and returns the plane normal of the first
/// `M` with `rank(M) ≤ 1`. Falls back to `A` and `B` themselves
/// (covers the degenerate input case where one surface is already
/// rank-1, e.g., a `Plane` quadric).
///
/// Returns `None` if no rank-1 fibre with rational `λ` exists. (A
/// rank-1 fibre at irrational `λ` would need the `Q[√d]` extension —
/// out of scope for this stage; in practice CAD inputs round-trip to
/// rational quadrics so the rational-root search catches the
/// engineering-relevant cases.)
pub(crate) fn find_tangent_plane(a: &Mat4SymRat, b: &Mat4SymRat) -> Option<[Rat; 4]> {
    let coefs = char_poly_rat(a, b);
    let roots = rational_roots(&coefs);
    let one = Rat::one();
    for lam in roots {
        let m = linear_comb_rat(a, b, &lam, &one);
        if is_rank1_or_less(&m)
            && let Some(n) = extract_plane_from_rank1(&m)
        {
            return Some(n);
        }
    }
    if is_rank1_or_less(a)
        && let Some(n) = extract_plane_from_rank1(a)
    {
        return Some(n);
    }
    if is_rank1_or_less(b)
        && let Some(n) = extract_plane_from_rank1(b)
    {
        return Some(n);
    }
    None
}

/// Convenience: detect tangent intersection at the surface level.
/// Lifts both surfaces into the rational quadric form and runs
/// `find_tangent_plane`. Returns `None` for non-rational inputs
/// (e.g., NURBS surfaces with no quadric representation) or for pairs
/// without a tangent plane in the pencil.
pub(crate) fn surfaces_are_tangent(
    surf_a: &crate::geometry::surface::SurfaceGeometry,
    surf_b: &crate::geometry::surface::SurfaceGeometry,
) -> Option<[Rat; 4]> {
    let a = super::rational::surface_to_quadric_rat(surf_a)?;
    let b = super::rational::surface_to_quadric_rat(surf_b)?;
    find_tangent_plane(&a, &b)
}
