//! Shared helpers for the trig-pair analytic dispatch.

use popkern_math::PortableMath;

/// Decompose a 2×2 matrix `M = [[a, b], [c, d]]` such that
///
/// ```text
/// M · (cos α, sin α) = R(rotation) · (u_radius · cos β, v_radius · sin β)
/// ```
///
/// for **all** `α`, where `β = α − phase` is a phase-shifted
/// parameter. Returns `(u_radius, v_radius, rotation, phase)`.
///
/// This is the parameterisation-preserving form of the 2×2 SVD
/// `M = U · Σ · Vᵀ`:
/// - `R(rotation) = U` is the left singular vector frame,
/// - `R(phase)ᵀ = Vᵀ` is the right singular vector frame,
/// - `(u_radius, v_radius)` are the singular values, signed if `det M < 0`
///   (the reflection part of `U` is absorbed into a negative
///   `v_radius` so callers can still treat `R(rotation)` as a
///   proper rotation).
///
/// Consumers building a `Pcurve::elliptical(center, u_r, v_r, rot,
/// start_angle, end_angle, t_start, t_end)` from this output should
/// pass `start_angle = α_start − phase` and `end_angle = α_end − phase`
/// so that the pcurve evaluates to `M · (cos α, sin α)` at
/// parameter `α`. If a caller assumed `phase = 0` (earlier code
/// path), they got correct output only when the matrix's columns
/// were already perpendicular and `m11 ≥ m22` with positive
/// determinant — the fix surfaces the bug for all other inputs.
///
/// Robust to rank-1 / zero matrices through guarded normalisation;
/// near-degenerate cases default to `phase = rotation = 0` and
/// `u_r = v_r = trace/2`-style equal singular values.
pub(crate) fn decompose_2x2(a: f64, b: f64, c: f64, d: f64) -> (f64, f64, f64, f64) {
    // M^T M is symmetric 2×2; its eigenvalues are the squared
    // singular values of M.
    let m11 = a * a + c * c;
    let m22 = b * b + d * d;
    let m12 = a * b + c * d;
    let trace = m11 + m22;
    let det = m11 * m22 - m12 * m12;
    let disc = ((trace * trace) / 4.0 - det).max(0.0).sqrt();
    let s1_sq = trace / 2.0 + disc;
    let s2_sq = (trace / 2.0 - disc).max(0.0);
    let s1 = s1_sq.sqrt();
    let s2 = s2_sq.sqrt();

    // Right singular vector for s1: (M^T M - s1²I) v = 0. Using
    // row 1 yields v ∝ (-m12, m11 - s1²); using row 2 yields v ∝
    // (s1² - m22, m12). Both are mathematically equivalent, but
    // they suffer different FP cancellation depending on which of
    // m11/m22 is closer to s1². When m12 ≈ 0 and m11 ≈ s1², row 1
    // collapses to (0, ε) with ε purely numerical noise — that
    // would steer the eigenvector toward the WRONG axis (the s2
    // eigenvector instead of the s1). Pick whichever row produces
    // the larger pre-normalisation magnitude.
    //
    // Threshold note: the eigenvector formula coefficients have
    // magnitudes ≲ |M| (the spectral norm). When `disc` (the SVD
    // discriminant) is within FP noise of 0 — the near-degenerate
    // case where M is essentially a scaled rotation — both formulas
    // give zero-norm vectors and any rotation is valid. We default
    // to identity in that case so callers building pcurves anchored
    // at α = 0 get a consistent (cos α, sin α) parameterisation.
    let scale = m11.max(m22).max(m12.abs()).max(1.0);
    let degeneracy_threshold = scale * 1.0e-12;
    let v1_a = (-m12, m11 - s1_sq);
    let v1_b = (s1_sq - m22, m12);
    let n_a = (v1_a.0 * v1_a.0 + v1_a.1 * v1_a.1).sqrt();
    let n_b = (v1_b.0 * v1_b.0 + v1_b.1 * v1_b.1).sqrt();
    let v1 = if n_a < degeneracy_threshold && n_b < degeneracy_threshold {
        // Near-degenerate (eigenvalues coincident) — pick identity.
        (1.0, 0.0)
    } else if n_a >= n_b {
        (v1_a.0 / n_a, v1_a.1 / n_a)
    } else {
        (v1_b.0 / n_b, v1_b.1 / n_b)
    };

    // Left singular vector u1 = (1/s1) · M · v1.
    let u1 = if s1 > 1.0e-30 {
        let ux = (a * v1.0 + b * v1.1) / s1;
        let uy = (c * v1.0 + d * v1.1) / s1;
        (ux, uy)
    } else {
        (1.0, 0.0)
    };
    let rotation = u1.1.patan2(u1.0);
    let phase = v1.1.patan2(v1.0);

    // Determinant sign of M: positive ⇒ proper rotation U; negative
    // ⇒ U is a reflection. Since `Pcurve::elliptical` only takes a
    // single `rotation` angle, we encode the reflection by negating
    // `v_radius`. Pcurve then evaluates `(u_r·cos β, −|v_r|·sin β)`
    // which traces the same ellipse with the orientation needed.
    let det_m = a * d - b * c;
    let v_radius = if det_m < 0.0 { -s2 } else { s2 };

    (s1, v_radius, rotation, phase)
}
