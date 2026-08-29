//! Shared ring-quality metrics for the K-P bi-nonic promotion —
//! worst-case spoke coplanarity (G1) and spoke gap (G0) over the
//! promotion's standard 50-sample grid. Split from
//! `three_quad_kp_promote.rs` for the 500-line structural cap.

use popkern_core::geometry::gregory::GregoryPatch;

/// Worst-case G1 coplanarity across the full ring's shared spokes,
/// using 50 uniform samples on `t ∈ [0.05, 0.95]` per spoke. The
/// endpoint margin matches the commercial-gate test's convention
/// (avoids the rational-Gregory corner-singularity finite-angle
/// transition zone). Returns 0 if normals are ill-conditioned at
/// every sample (can't happen on a real ring, but guards against
/// log-scale division).
pub(crate) fn worst_ring_coplanarity(patches: &[GregoryPatch]) -> f64 {
    let n = patches.len();
    let mut worst = 0.0_f64;
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let t_b = 1.0 - t;
            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let cross = xtan_a.cross(&xtan_b);
            let numer = (cross.x * t_spoke.x + cross.y * t_spoke.y + cross.z * t_spoke.z).abs();
            let na = (xtan_a.x * xtan_a.x + xtan_a.y * xtan_a.y + xtan_a.z * xtan_a.z).sqrt();
            let nb = (xtan_b.x * xtan_b.x + xtan_b.y * xtan_b.y + xtan_b.z * xtan_b.z).sqrt();
            let ns = (t_spoke.x * t_spoke.x + t_spoke.y * t_spoke.y + t_spoke.z * t_spoke.z).sqrt();
            let denom = na * nb * ns;
            if denom > 1e-300 {
                let coplan = numer / denom;
                if coplan > worst {
                    worst = coplan;
                }
            }
        }
    }
    worst
}

/// Worst-case G0 gap across the ring's shared spokes on the same
/// 50-sample grid as [`worst_ring_coplanarity`].
pub(crate) fn worst_ring_g0_gap(patches: &[GregoryPatch]) -> f64 {
    let n = patches.len();
    let mut worst = 0.0_f64;
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let t_b = 1.0 - t;
            let pa = a.point_at(1.0, t);
            let pb = b.point_at(t_b, 0.0);
            let gap = pa.distance_to(&pb);
            if gap > worst {
                worst = gap;
            }
        }
    }
    worst
}
