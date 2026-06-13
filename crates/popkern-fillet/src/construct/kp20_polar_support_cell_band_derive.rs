//! Row-1 band derivation for the constructive support cells:
//! `D(t) = α(t)·C'(t) + β(t)·R(t)` with linear `α, β` matched to
//! the exact end legs, expanded into degree-9 Bezier
//! coefficients. Degenerate corners (chain-collinear legs:
//! T-junction stems at 180 degrees, sliver tips near 0) derive
//! one-sided from the healthy end; a transverse floor guards the
//! mid-row against parametric collapse.

use popkern_core::geometry::{Point3, Vec3};

use super::chain_ribbon_field::row_derivative_coefficients;

/// Determinant floor for the 2x2 end-leg solves (relative to the
/// basis vectors' scale).
const LEG_SOLVE_DET_TOL: f64 = 1.0e-10;

/// Relative floor on the band's transverse (β·R) component: the
/// weak end's transverse magnitude is clamped to this fraction
/// of the strong end's, with matching sign. Without it, sliver
/// tips whose legs tilt almost tangent to the chain make β(t)
/// cross zero mid-row — ∂u ∥ ∂v there, a parametrically
/// degenerate stretch whose normals are noise on both sides of
/// the trim. The clamp keeps everything in-span (only the
/// in-plane shape near the corner changes), so G1 is unaffected
/// by construction.
const TRANSVERSE_FLOOR: f64 = 0.6;

/// Sine bound under which an end leg counts as chain-collinear
/// (a degenerate corner). T-junction stem vertices have a
/// 180-degree interior angle BY CONSTRUCTION (a vertex inserted
/// mid-edge of a straight boundary), so their legs are exactly
/// anti-parallel to the chain; sliver tips approach 0 degrees.
/// Interpolating such a leg forces a tangent-dominated band
/// whose trim normals are noise — those corners derive
/// one-sided from the healthy end instead.
const DEGENERATE_CORNER_SIN_TOL: f64 = 0.1;

/// Derive the full row-1 band from `(C', R)` and the exact end
/// legs: `D(t) = α(t)·C'(t) + β(t)·R(t)` with linear `α, β`
/// matched to `D(0) = E₀`, `D(1) = E₁`, expanded into degree-9
/// Bezier coefficients, then `row1ₖ = row0ₖ + Dₖ / 9`.
pub(super) fn derive_row1(
    row0: &[Point3; 10],
    row1: &[Point3; 10],
    field_segment: &[Vec3; 9],
    clamped_count: &mut usize,
    one_sided_count: &mut usize,
) -> Option<[Point3; 10]> {
    let derivative = row_derivative_coefficients(row0);
    let leg_start = (row1[0] - row0[0]) * 9.0;
    let leg_end = (row1[9] - row0[9]) * 9.0;
    // One-sided derivation at degenerate corners: a chain-
    // collinear leg cannot define a transverse blend; copying the
    // healthy end's transverse magnitude (scaled across the field)
    // keeps the band in span AND non-tangential. Both crosses of
    // the trim stay in span(C', R), so pair G1 is unaffected.
    let degenerate_start = leg_is_chain_collinear(leg_start, derivative[0]);
    let degenerate_end = leg_is_chain_collinear(leg_end, derivative[8]);
    let field_start_mag = field_segment[0].dot(&field_segment[0]).sqrt().max(1.0e-300);
    let field_end_mag = field_segment[8].dot(&field_segment[8]).sqrt().max(1.0e-300);
    let (alpha_start, mut beta_start, alpha_end, mut beta_end) =
        match (degenerate_start, degenerate_end) {
            (true, false) => {
                let (a1, b1) = solve_leg(leg_end, derivative[8], field_segment[8])?;
                *one_sided_count += 1;
                (
                    chain_projection(leg_start, derivative[0]),
                    b1 * field_end_mag / field_start_mag,
                    a1,
                    b1,
                )
            }
            (false, true) => {
                let (a0, b0) = solve_leg(leg_start, derivative[0], field_segment[0])?;
                *one_sided_count += 1;
                (
                    a0,
                    b0,
                    chain_projection(leg_end, derivative[8]),
                    b0 * field_start_mag / field_end_mag,
                )
            }
            _ => {
                let (a0, b0) = solve_leg(leg_start, derivative[0], field_segment[0])?;
                let (a1, b1) = solve_leg(leg_end, derivative[8], field_segment[8])?;
                (a0, b0, a1, b1)
            }
        };
    // Transverse-floor clamp: each end's transverse magnitude
    // must reach a fraction of the LEG SCALE (not merely of the
    // other end — rows whose both legs are near-tangent would
    // otherwise stay collapsed), with one consistent sign taken
    // from the dominant end (or positive when both vanish).
    let leg_scale = leg_start
        .dot(&leg_start)
        .sqrt()
        .max(leg_end.dot(&leg_end).sqrt());
    let floor_magnitude = TRANSVERSE_FLOOR * leg_scale;
    let transverse_start = beta_start * field_segment[0].dot(&field_segment[0]).sqrt();
    let transverse_end = beta_end * field_segment[8].dot(&field_segment[8]).sqrt();
    let dominant_sign = if transverse_start.abs() >= transverse_end.abs() {
        if transverse_start >= 0.0 { 1.0 } else { -1.0 }
    } else if transverse_end >= 0.0 {
        1.0
    } else {
        -1.0
    };
    // Measured note (PK-085): restricting this clamp to
    // near-collinear legs collapses the worst span residual
    // (0.39 -> 0.12, clamps 90 -> 15) but reads 7.10e-3 vs
    // 6.84e-3 at the gate — the broad clamp stays on measured
    // acceptance; revisit when the binding pair changes.
    if transverse_start * dominant_sign < floor_magnitude {
        let magnitude = field_segment[0].dot(&field_segment[0]).sqrt().max(1.0e-300);
        beta_start = dominant_sign * floor_magnitude / magnitude;
        *clamped_count += 1;
    }
    if transverse_end * dominant_sign < floor_magnitude {
        let magnitude = field_segment[8].dot(&field_segment[8]).sqrt().max(1.0e-300);
        beta_end = dominant_sign * floor_magnitude / magnitude;
        *clamped_count += 1;
    }
    let mut out = *row1;
    for k in 1..9 {
        // Product of a linear blend with a degree-8 field:
        // coeffₖ = (9−k)/9·(α₀ d_k + β₀ r_k)
        //        +    k/9·(α₁ d_{k−1} + β₁ r_{k−1}).
        let head = derivative[k] * alpha_start + field_segment[k] * beta_start;
        let tail = derivative[k - 1] * alpha_end + field_segment[k - 1] * beta_end;
        let d_k = head * ((9 - k) as f64 / 9.0) + tail * (k as f64 / 9.0);
        out[k] = row0[k] + d_k * (1.0 / 9.0);
    }
    Some(out)
}

/// Whether a corner leg is (anti)parallel to the chain tangent
/// within [`DEGENERATE_CORNER_SIN_TOL`].
fn leg_is_chain_collinear(leg: Vec3, chain: Vec3) -> bool {
    let cross = leg.cross(&chain);
    let denominator = (leg.dot(&leg) * chain.dot(&chain)).sqrt();
    denominator > 1.0e-300 && cross.dot(&cross).sqrt() / denominator < DEGENERATE_CORNER_SIN_TOL
}

/// Projection coefficient of a degenerate leg onto the chain
/// tangent (its alpha; the leg has no usable transverse part).
fn chain_projection(leg: Vec3, chain: Vec3) -> f64 {
    let chain_squared = chain.dot(&chain);
    if chain_squared <= 1.0e-300 {
        0.0
    } else {
        leg.dot(&chain) / chain_squared
    }
}

/// Solve `leg = α·u + β·v` in the least-squares sense (exact
/// when the leg is coplanar with the basis, which the node
/// tangent-plane constraints guarantee).
pub(super) fn solve_leg(leg: Vec3, u: Vec3, v: Vec3) -> Option<(f64, f64)> {
    let uu = u.dot(&u);
    let uv = u.dot(&v);
    let vv = v.dot(&v);
    let det = uu * vv - uv * uv;
    if det.abs() <= LEG_SOLVE_DET_TOL * uu.max(vv).max(1.0e-300) {
        return None;
    }
    let bu = u.dot(&leg);
    let bv = v.dot(&leg);
    Some(((vv * bu - uv * bv) / det, (uu * bv - uv * bu) / det))
}
