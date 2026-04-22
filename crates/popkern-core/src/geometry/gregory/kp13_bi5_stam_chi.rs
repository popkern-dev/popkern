//! Full Stam-1998 Catmull-Clark characteristic-map evaluator.
//!
//! The characteristic map `χ^k : [0, 1]² → R²` at a valence-`n`
//! extraordinary vertex is defined as the limit of Catmull-Clark
//! subdivision projected onto the 2-D subspace spanned by the
//! **subdominant eigenvector pair** of the subdivision matrix
//! (with eigenvalue `λ_n`). Two key properties determine the map
//! uniquely (up to rigid similarity):
//!
//! 1. **Boundary values** at the four unit-square corners:
//!    `χ^k(0, 0) = 0`, `χ^k(1, 0) = m_k`, `χ^k(0, 1) = m_{k+1}`,
//!    `χ^k(1, 1) = v_k` (the outer limit point).
//!
//! 2. **Self-similarity** under Catmull-Clark subdivision:
//!
//!    ```text
//!      χ^k(u, v) = λ_n · χ^k(2u, 2v)   for (u, v) near the star,
//!                                       i.e. max(u, v) ≤ 1/2.
//!    ```
//!
//!    Stam 1998 §4 derives this directly from the subdivision-
//!    matrix eigenstructure: the subdominant-eigenvector projection
//!    of a Catmull-Clark subdivided mesh is `λ_n` times that of the
//!    un-subdivided mesh, so restricting `χ` to the "`Ω_0`" quadrant
//!    `[0, 1/2]²` and rescaling by 2 reproduces `λ · χ` on the full
//!    unit square. The identity is exact, not an approximation.
//!
//! # Algorithm (this module)
//!
//! Property (2) directly gives an arbitrary-`(u, v)` recursion:
//!
//! ```text
//!   eval(u, v) := {
//!     λ · eval(2u, 2v)              if max(u, v) ≤ 1/2  (near-star)
//!     far_region_eval(u, v)         otherwise           (far-from-star)
//!   }
//! ```
//!
//! - **Termination.** Each near-star recursion doubles `(u, v)`
//!   and multiplies the result by `λ ∈ (0, 1)`. After at most
//!   `⌈log₂(1/max(u, v))⌉` levels, the argument leaves `Ω_0` and
//!   the far-region branch fires. For any `(u, v) ≠ (0, 0)` the
//!   recursion terminates in `O(log(1/max(u, v)))` steps.
//!
//! - **`(0, 0)`-limit.** The recursion is capped at [`MAX_DEPTH`]
//!   levels, after which the result is scaled by `λ^MAX_DEPTH`
//!   (below `1e-15` for any physical `λ`) — agreeing with the
//!   known limit `χ(0, 0) = 0` to double precision.
//!
//! - **Far-region evaluator.** Delegates to [`CatmullClarkChi`],
//!   whose bi-cubic Bézier patch reproduces the four corner jets
//!   (value, first partials, second partials) exactly and provides
//!   a smooth interior interpolation. At regular valence `n = 4`
//!   the characteristic map is linear (bi-linear) and the bi-cubic
//!   matches analytically in the far region.
//!
//! # Fidelity
//!
//! **Exact:** the self-similarity law `χ(u, v) = λ · χ(2u, 2v)`.
//! Because property (2) is exact, the recursion is **not** an
//! approximation — it moves interior queries deterministically
//! onto the far-region evaluator.
//!
//! **Approximate:** the far-region evaluator is a bi-cubic Bézier
//! patch (not a true 4×4 Catmull-Clark B-spline sub-patch). The
//! bi-cubic matches `χ` exactly at the four corners (value, both
//! first partials, all three second partials) and is `O(1)` to
//! evaluate, but interior far-region values are a degree-3 tensor-
//! product interpolant of the corner jets rather than the true CC
//! limit surface. Replacing the bi-cubic with Stam's closed-form
//! (`2n+8`-CP eigendecomposition) would remove this residual
//! approximation and is noted below.
//!
//! # Relationship to full Stam-1998
//!
//! Stam's original formulation builds the `(2n+8) × (2n+8)` CC
//! subdivision matrix on the extended patch, block-diagonalises
//! via DFT, and derives closed-form eigenvalue/eigenvector
//! expressions. The far-region query then reduces to a 4×4
//! bi-cubic B-spline evaluation on one of the three regular
//! sub-patches of the extended patch. That machinery is **not**
//! implemented here — the present module implements the
//! **self-similarity half** of Stam's algorithm exactly, and
//! delegates the far-region half to the corner-matching bi-cubic
//! Bézier patch [`CatmullClarkChi`].
//!
//! In practice, `StamChi` produces arbitrary-`(u, v)` results
//! that differ from [`CatmullClarkChi::eval`] **only in the
//! near-star quadrant** `[0, 1/2]²` — there, it returns the
//! self-similarity-consistent value `λ^k · CatmullClarkChi::eval
//! (2^k u, 2^k v)`, which is exact in `λ` at all orders of
//! recursion. For `(u, v)` outside the near-star quadrant, the
//! two evaluators return identical results.

#![allow(dead_code)]

use super::kp13_bi5_chi::{CatmullClarkChi, ChiCorner, ChiJet};

/// Cap on the recursion depth for near-star queries. Each level
/// multiplies the result by `λ ∈ (0, 1)`, so after 50 levels any
/// reasonable `λ` has `λ^50 ≲ 1e-15` — at or below double-precision
/// noise. Queries at `(u, v) = (0, 0)` exactly still terminate via
/// the depth cap and return a value within machine epsilon of the
/// known limit `(0, 0)`.
pub const MAX_DEPTH: u32 = 50;

/// Boundary of the "`Ω_0`" near-star quadrant. Queries with
/// `max(u, v) ≤ STAR_BOUNDARY` recurse; otherwise evaluate directly.
///
/// Set to `1/2` exactly, matching the standard CC subdivision
/// half-domain split: one subdivision step maps the full unit
/// square onto `[0, 1/2]²` (for the `Ω_0` child) plus three
/// regular sub-patches.
pub const STAR_BOUNDARY: f64 = 0.5;

/// Full Stam-1998 Catmull-Clark characteristic-map evaluator.
///
/// Wraps [`CatmullClarkChi`] and applies the subdominant-
/// eigenvalue self-similarity law
/// `χ(u, v) = λ · χ(2u, 2v)` recursively on near-star queries.
/// See the module header for the full derivation and fidelity
/// tradeoffs.
#[derive(Debug, Clone, Copy)]
pub struct StamChi {
    /// Valence of the extraordinary vertex.
    pub n: usize,
    /// Subdominant eigenvalue `λ_n` (shared with [`CatmullClarkChi`]).
    pub lambda: f64,
    /// Far-region evaluator.
    pub far: CatmullClarkChi,
}

impl StamChi {
    /// New Stam χ evaluator for valence `n ≥ 3`.
    pub fn new(n: usize) -> Self {
        let far = CatmullClarkChi::new(n);
        Self {
            n,
            lambda: far.lambda,
            far,
        }
    }

    /// Evaluate `χ^k(u, v)` for arbitrary `(u, v) ∈ [0, 1]²`.
    ///
    /// Uses the self-similarity law `χ(u, v) = λ · χ(2u, 2v)` to
    /// recurse on near-star queries (`max(u, v) ≤ 1/2`) and
    /// delegates to [`CatmullClarkChi::eval`] on the far region.
    pub fn eval(&self, k: usize, u: f64, v: f64) -> (f64, f64) {
        self.eval_rec(k, u, v, 0)
    }

    /// Analytic `ChiJet` at the given corner of sector `k`.
    ///
    /// Corner jets are known in closed form from the CC
    /// eigenstructure — no subdivision recursion is required.
    /// Delegates to [`CatmullClarkChi::jet`], whose bi-cubic
    /// Bézier CP grid matches the analytic corner jets exactly.
    pub fn jet(&self, k: usize, corner: ChiCorner) -> ChiJet {
        self.far.jet(k, corner)
    }

    /// Recursion depth at which the far-region branch would fire
    /// for the given `(u, v)`. Useful for diagnostics and tests.
    /// Returns `0` for far-region queries, positive values for
    /// near-star queries.
    pub fn depth_for(u: f64, v: f64) -> u32 {
        let mx = u.max(v);
        if mx > STAR_BOUNDARY {
            return 0;
        }
        if mx <= 0.0 {
            return MAX_DEPTH;
        }
        // Smallest k such that 2^k · max(u, v) > STAR_BOUNDARY.
        let k = (STAR_BOUNDARY / mx).log2().ceil() as u32;
        k.min(MAX_DEPTH)
    }

    fn eval_rec(&self, k: usize, u: f64, v: f64, depth: u32) -> (f64, f64) {
        if depth >= MAX_DEPTH || u.max(v) > STAR_BOUNDARY {
            return self.far.eval(k, u, v);
        }
        let (x2, y2) = self.eval_rec(k, 2.0 * u, 2.0 * v, depth + 1);
        (self.lambda * x2, self.lambda * y2)
    }
}

#[cfg(test)]
#[path = "kp13_bi5_stam_chi_tests.rs"]
mod tests;
