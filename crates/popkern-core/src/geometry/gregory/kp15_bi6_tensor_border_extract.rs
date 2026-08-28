//! Extract a bi-3 [`TensorBorder`] from a
//! [`PolarBoundaryInput`].
//!
//! The K-P 2015 bi-6 polar cap solver consumes **depth-2 Hermite
//! data** along each sector's outer edge ([`TensorBorderBi3`]).
//! For our fillet-engine use case, that data isn't directly
//! provided — only the geometric inputs (spokes, face tangent
//! planes, per-spoke radii) are. This module bridges the gap with
//! a pragmatic extractor that produces a valid non-degenerate
//! tensor-border.
//!
//! # Semantics
//!
//! For sector `k`, the tensor-border rows represent BB-CPs of a
//! bi-3 patch patch-edge that the polar cap joins along. Each
//! row has 4 CPs indexed `i ∈ 0..3`.
//!
//! - **Row 0** (position): samples along spoke `k` at BB-CP
//!   positions (end-interpolating the spoke's degree-7 Bezier).
//!   For a bi-3 representation of a curve, the BB-CPs are at
//!   specific positions — we degree-reduce from the degree-7
//!   input via equally-spaced parameter samples plus end-point
//!   anchoring.
//!
//! - **Row 1** (first cross-derivative): offsets each row-0 CP
//!   by `(1/3)·r_k · n̂_k` where `r_k` is the linearly-
//!   interpolated radius along the spoke and `n̂_k` is the
//!   sector `k` face tangent plane's `v` basis vector. The
//!   factor `1/n = 1/3` reflects the bi-3 BB derivative
//!   convention: `∂f/∂v = 3·(row1 − row0)`, so row1 − row0
//!   encodes `(1/3)·∂f/∂v`.
//!
//! - **Row 2** (second cross-derivative): offsets row-0 by
//!   `(2/3)·r_k · n̂_k`. This is a **linear** approximation
//!   assuming zero second derivative — matches a flat extension
//!   of the tensor-border. Curvature-aware extraction will
//!   replace this when incident-surface curvature data is
//!   plumbed through (future work).
//!
//! # What this produces vs. what KP15 ideally wants
//!
//! KP15 expects the tensor-border to encode second-order contact
//! between the cap and the surrounding `b` surface. Our flat
//! approximation at row 2 means the cap will have approximately
//! zero cross-curvature contribution from row 2 — good enough
//! for G1 continuity but not yet full G2 match to the incident
//! surface. F_5 minimisation (I.10b) still produces a well-posed
//! optimal cap.
//!
//! The full-fidelity extraction requires sampling ∂²/∂v² of the
//! fillet strip surface at each BB-CP position — which in turn
//! requires access to the fillet strip geometry (cylinder/torus
//! curvatures for the standard cases). Deferred.

use super::kp_polar_boundary_input::{PolarBoundaryInput, TensorBorder, TensorBorderBi3};
use crate::geometry::{Point3, Vec3};

/// Extract a bi-3 [`TensorBorder`] from a polar boundary input,
/// using the simple linear-offset heuristic described in the
/// module docs.
///
/// # Errors
/// - `polar.validate()` fails (non-finite radii, length mismatch).
/// - `n < 3` (can't form a polar cap).
pub fn extract_bi3_tensor_border_from_polar_input(
    polar: &PolarBoundaryInput,
) -> Result<TensorBorder, String> {
    polar.validate()?;
    let n = polar.n();
    let mut sectors: Vec<[[Point3; 4]; 3]> = Vec::with_capacity(n);
    for k in 0..n {
        let spoke = &polar.spokes[k];
        let row0 = bi3_row_from_spoke(&spoke.boundary_cps);

        // Offset direction: face tangent plane v-axis, expressed
        // in world coordinates. Sign convention: use +v (callers
        // can flip by replacing with -v if needed — the bi-6
        // solver is sign-agnostic for scaffold-level use).
        let n_hat: Vec3 = polar.face_tangent_planes[k].v;

        // Row-1 and row-2 offsets per paper bi-3 convention:
        //   row1 = row0 + (1/3) · r · n̂
        //   row2 = row0 + (2/3) · r · n̂
        // Radius interpolated linearly from radius_at_c to
        // radius_at_v across the 4 BB-CP positions.
        let r_c = spoke.radius_at_c;
        let r_v = spoke.radius_at_v;
        // BB-CP "positions" for a bi-3 curve correspond to
        // parameter values t_i = i / 3 in de Casteljau frames.
        let radii_at_cps = [r_c, (2.0 * r_c + r_v) / 3.0, (r_c + 2.0 * r_v) / 3.0, r_v];

        let mut row1 = row0;
        let mut row2 = row0;
        for i in 0..4 {
            let r = radii_at_cps[i];
            let delta = n_hat * (r / 3.0);
            row1[i] = row0[i] + delta;
            row2[i] = row0[i] + delta + delta; // 2·delta
        }
        sectors.push([row0, row1, row2]);
    }
    Ok(TensorBorder::Bi3(TensorBorderBi3 { sectors }))
}

/// Degree-reduce a degree-7 Bezier curve (8 CPs) to degree 3
/// (4 CPs), preserving the two endpoints exactly and picking
/// the two interior CPs via parameter-value evaluation at
/// `t = 1/3` and `t = 2/3` (de Casteljau). This is a pragmatic
/// sampling approach rather than a least-squares projection —
/// the bi-6 cap solver tolerates the introduced boundary error
/// because F_5 minimisation refines the overall shape.
fn bi3_row_from_spoke(spoke_cps_bi7: &[Point3; 8]) -> [Point3; 4] {
    let p0 = spoke_cps_bi7[0];
    let p3 = spoke_cps_bi7[7];
    let p1 = eval_bi7_at(spoke_cps_bi7, 1.0 / 3.0);
    let p2 = eval_bi7_at(spoke_cps_bi7, 2.0 / 3.0);
    // Degree-3 BB-CPs reproduce a cubic Bezier; the sampled
    // points are on-curve values, not BB-CPs of the original.
    // For the scaffold's purposes, we treat the sampled points
    // as the 4 BB-CPs of a new cubic that approximates the
    // degree-7 curve. The resulting cubic won't exactly
    // reproduce the degree-7 curve, but the endpoint match +
    // two interior samples captures the rough shape.
    [p0, p1, p2, p3]
}

/// Evaluate a degree-7 Bezier curve (8 CPs) at `t` via
/// de Casteljau.
fn eval_bi7_at(cps: &[Point3; 8], t: f64) -> Point3 {
    let mut work: [Point3; 8] = *cps;
    for r in 1..=7 {
        for i in 0..(8 - r) {
            let a = work[i];
            let b = work[i + 1];
            work[i] = Point3::new(
                (1.0 - t) * a.x + t * b.x,
                (1.0 - t) * a.y + t * b.y,
                (1.0 - t) * a.z + t * b.z,
            );
        }
    }
    work[0]
}

#[cfg(test)]
#[path = "kp15_bi6_tensor_border_extract_tests.rs"]
mod tests;
