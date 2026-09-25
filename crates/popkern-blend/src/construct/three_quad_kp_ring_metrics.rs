//! Shared ring-quality metrics for the K-P bi-nonic promotion —
//! worst-case spoke normal sine (G1) and spoke gap (G0) over the
//! promotion's standard 50-sample grid. Split from
//! `three_quad_kp_promote.rs` for the 500-line structural cap.

use popkern_core::geometry::gregory::GregoryPatch;

/// Worst-case G1 across the full ring's shared spokes — the sine of
/// the angle between the two patches' normals (PK-539; see
/// [`crate::construct::kp_spoke_normal_sine`] for why the normalized
/// triple product it replaces was blind) — using 50 uniform samples
/// on `t ∈ [0.05, 0.95]` per spoke. The endpoint margin matches the
/// commercial-gate test's convention (avoids the rational-Gregory
/// corner-singularity finite-angle transition zone).
pub(crate) fn worst_ring_normal_sine(patches: &[GregoryPatch]) -> f64 {
    let n = patches.len();
    let mut worst = 0.0_f64;
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        let samples: Vec<_> = (0..50)
            .map(|sample| {
                let t = 0.05 + 0.90 * (sample as f64) / 49.0;
                (
                    a.partial_u(1.0, t),
                    b.partial_v(1.0 - t, 0.0),
                    a.partial_v(1.0, t),
                )
            })
            .collect();
        worst = worst.max(crate::construct::kp_spoke_normal_sine::worst_spoke_measure(
            &samples,
        ));
    }
    worst
}

/// Worst-case G0 gap across the ring's shared spokes on the same
/// 50-sample grid as [`worst_ring_normal_sine`].
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
