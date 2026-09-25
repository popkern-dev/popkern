//! Shared ring-quality metrics for the K-P bi-nonic promotion —
//! worst-case spoke normal sine (G1) and spoke gap (G0) over the
//! promotion's standard 50-sample grid. Split from
//! `three_quad_kp_promote.rs` for the 500-line structural cap.

use popkern_core::geometry::gregory::GregoryPatch;

/// Worst-case G1 across the full ring's shared spokes — the sine of
/// the angle between the two patches' normals (PK-539; see
/// [`crate::construct::kp_spoke_normal_sine`] for why the normalized
/// triple product it replaces was blind) — over the whole spoke, end
/// zones included
/// ([`crate::construct::kp_spoke_normal_sine::whole_spoke_parameters`],
/// PK-547). The former
/// `[0.05, 0.95]` window was justified by the rational-Gregory corner
/// transition, but the promoted rings are polynomial bi-nonics, and the
/// window hid a 5%-long crease at the rim end of the served
/// extreme-ladder rings.
pub(crate) fn worst_ring_normal_sine(patches: &[GregoryPatch]) -> f64 {
    let n = patches.len();
    let parameters = crate::construct::kp_spoke_normal_sine::whole_spoke_parameters();
    let mut worst = 0.0_f64;
    for k in 0..n {
        let a = &patches[k];
        let b = &patches[(k + 1) % n];
        let samples: Vec<_> = parameters
            .iter()
            .map(|&t| {
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

/// The continuity a committed ring may claim: G1 when its
/// [`worst_ring_normal_sine`] (whole spokes, end zones included) is
/// within the commercial bar, else G0 (PK-539 measured claims,
/// extended to the ends under PK-544/PK-547).
pub(crate) fn measured_ring_continuity(patches: &[GregoryPatch]) -> crate::ContinuityLevel {
    if worst_ring_normal_sine(patches) <= crate::construct::kp_spoke_normal_sine::RING_SPOKE_G1_GATE
    {
        crate::ContinuityLevel::G1
    } else {
        crate::ContinuityLevel::G0
    }
}
