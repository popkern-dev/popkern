//! The ring metric: continuity across the seams of a Gregory ring
//! whose patch k's u = 1 side meets patch k + 1's v = 0 side reversed
//! (the polar/three-quad layout). Split out of `kp_gate_helpers.rs`
//! for the file budget (PK-498).

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::GregoryPatch;

use super::{dot, normal_curvature_along, vec_norm};

/// Worst ring-seam residuals (PK-498 honest measures).
/// - `worst_g1_normal` is the sin of the angle between the two sides'
///   normals.
/// - `worst_g2_residual` is the normal-curvature jump along, across and
///   diagonal to the seam, normalized by `1 + κ_max`.
pub struct RingMetrics {
    pub worst_g0_gap: f64,
    pub worst_g1_normal: f64,
    pub worst_g2_residual: f64,
    /// Per ring seam k (patch k's u = 1 against patch k + 1's v = 0):
    /// worst `(g1, g2)`.
    pub per_seam: Vec<(f64, f64)>,
}

/// Sample the ring on `t ∈ [0.05, 0.95]` at `n_samples` points.
/// Endpoint margins exclude the Gregory corner-singularity
/// transition zones — not a bi-9 construction concern.
pub fn measure_ring(patches: &[GregoryPatch], n_samples: usize) -> RingMetrics {
    let n = patches.len();
    let mut worst_g0 = 0.0_f64;
    let mut worst_g1 = 0.0_f64;
    let mut worst_g2 = 0.0_f64;
    let ts: Vec<f64> = (0..n_samples)
        .map(|k| 0.05 + 0.90 * (k as f64) / ((n_samples - 1) as f64))
        .collect();
    let mut per_seam = Vec::with_capacity(n);
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        let (mut seam_g1, mut seam_g2) = (0.0_f64, 0.0_f64);
        for &t in &ts {
            let t_b = 1.0 - t;

            let pa = a.point_at(1.0, t);
            let pb = b.point_at(t_b, 0.0);
            let gap =
                ((pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2)).sqrt();
            if gap > worst_g0 {
                worst_g0 = gap;
            }

            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            // G1: sin of the angle between the two sides' normals, the
            // live-seam measure. The former cross-tangent/spoke
            // coplanarity triple product collapses when the cross
            // tangents are nearly collinear (a continuing seam): a
            // PK-498 control passed a 53° crease at 1.7e-2.
            // Skip samples where either chart basis is near-collinear
            // (the ring centre): the normal there is chart noise, the
            // same conditioning guard the G2 comparison uses.
            let (n_side_a, n_side_b) =
                (xtan_a.cross(&t_spoke), b.partial_u(t_b, 0.0).cross(&xtan_b));
            let basis_sine =
                |n: Vec3, p: Vec3, q: Vec3| vec_norm(n) / (vec_norm(p) * vec_norm(q)).max(1e-300);
            let conditioned = basis_sine(n_side_a, xtan_a, t_spoke) >= 5.0e-2
                && basis_sine(n_side_b, b.partial_u(t_b, 0.0), xtan_b) >= 5.0e-2;
            let denom = vec_norm(n_side_a) * vec_norm(n_side_b);
            if conditioned && denom > 1e-300 {
                let sine = vec_norm(n_side_a.cross(&n_side_b)) / denom;
                seam_g1 = seam_g1.max(sine);
                if sine > worst_g1 {
                    worst_g1 = sine;
                }
            }

            let Some(t_unit) = t_spoke.normalized() else {
                continue;
            };
            // Fold the two charts' natural normals to a consistent
            // side before comparing II: a chart-handedness reversal
            // flips the SIGN of the measured normal curvature while
            // the geometry stays curvature-continuous (measured on
            // the PK-240 balanced rings: kb = -ka to 2-3 digits with
            // |dk| under 4e-2, while the unfolded comparison read
            // 1.2e0 — the same unfolded-normal artifact class as
            // PK-244's flags). Chart-reversal bands themselves are
            // recorded separately; this metric measures GEOMETRY.
            let na = a.partial_u(1.0, t).cross(&a.partial_v(1.0, t));
            let nb = b.partial_u(t_b, 0.0).cross(&b.partial_v(t_b, 0.0));
            // Conditioning guard: II is unreadable through a chart
            // whose basis is near-collinear (measured on the PK-240
            // balanced rings: the worst "G2" sample sat at a 0.0-deg
            // basis angle with healthy factor magnitudes — the
            // absolute det floor cannot catch it). Skip, don't clamp.
            let sin_a = vec_norm(na)
                / (vec_norm(a.partial_u(1.0, t)) * vec_norm(a.partial_v(1.0, t))).max(1e-300);
            let sin_b = vec_norm(nb)
                / (vec_norm(b.partial_u(t_b, 0.0)) * vec_norm(b.partial_v(t_b, 0.0))).max(1e-300);
            if sin_a < 5.0e-2 || sin_b < 5.0e-2 {
                continue;
            }
            let denom = vec_norm(na) * vec_norm(nb);
            if denom < 1e-300 {
                continue;
            }
            let cos = dot(na, nb) / denom;
            // Convergence guard: near a chart-degeneracy crossing the
            // natural normal spins through zero and II reads noise —
            // the G1 gate (coplanarity, orientation-free) separately
            // covers plane agreement, so a low |cos| here can only be
            // chart noise, never geometry. Skip, don't clamp.
            if cos.abs() < 0.9 {
                continue;
            }
            // Compare II along the seam AND across it (transverse and
            // diagonal in the shared tangent plane). Along the seam,
            // agreement is implied by G0+G1 (Meusnier: κ_n along a
            // shared direction is c''·n), so the along-seam term alone
            // is blind to curvature continuity. The PK-498 negative
            // control measured it at 1.0e-12 on a seam whose
            // transverse curvature was perturbed far past the gate.
            let Some(n_a) = na.normalized() else {
                continue;
            };
            let Some(t_perp) = n_a.cross(&t_unit).normalized() else {
                continue;
            };
            let directions = [
                t_unit,
                t_perp,
                (t_unit + t_perp).normalized().unwrap_or(t_unit),
            ];
            for dir in directions {
                let (Some(ka), Some(kb)) = (
                    normal_curvature_along(a, 1.0, t, dir),
                    normal_curvature_along(b, t_b, 0.0, dir),
                ) else {
                    continue;
                };
                let kb = if cos < 0.0 { -kb } else { kb };
                let kmax = ka.abs().max(kb.abs());
                let resid = (ka - kb).abs() / (1.0 + kmax);
                seam_g2 = seam_g2.max(resid);
                if resid > worst_g2 {
                    worst_g2 = resid;
                }
            }
        }
        per_seam.push((seam_g1, seam_g2));
    }
    RingMetrics {
        worst_g0_gap: worst_g0,
        worst_g1_normal: worst_g1,
        worst_g2_residual: worst_g2,
        per_seam,
    }
}
