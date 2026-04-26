//! Shared helpers for the trig-pair analytic dispatch.

/// Decompose a 2×2 matrix `A = [[a,b],[c,d]]` into
/// `R(rotation) · diag(u_radius, v_radius)` where `R` is a rotation.
/// Returns `(u_radius, v_radius, rotation)` such that the ellipse
/// traced by `A · (cos α, sin α)` equals
/// `R · (u_r cos α, v_r sin α)`.
///
/// Uses a closed-form 2×2 SVD; rotation is the angle of the left
/// singular vector matching the larger singular value. Robust to
/// rank-1 / zero matrices through guarded normalisation.
pub(super) fn decompose_2x2(a: f64, b: f64, c: f64, d: f64) -> (f64, f64, f64) {
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

    // Left singular vector u1 = (1/s1) · A · v1.
    let u1 = if s1 > 1.0e-30 {
        let ux = (a * v1.0 + b * v1.1) / s1;
        let uy = (c * v1.0 + d * v1.1) / s1;
        (ux, uy)
    } else {
        (1.0, 0.0)
    };
    let rotation = u1.1.atan2(u1.0);

    (s1, s2, rotation)
}
