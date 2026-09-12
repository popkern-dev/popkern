//! Negative-branch / multi-wrap chart alignment for record pcurves
//! (split from `clip_windows.rs` under the 500-line structural
//! budget; measured rationale lives with the code — PK-480/PK-481).

use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

/// Whole-period translation aligning a record pcurve onto a face
/// chart that lives on a NEGATIVE u-branch — `None` when the face is
/// non-periodic, its boundary is centred at or above `−period/4`
/// (principal-range charts keep the tuned machinery), or no
/// whole-period shift is called for. See the alignment comment in
/// [`clip_curve_uv_space`] (PK-480).
pub(super) fn align_negative_branch_side(
    brep: &BRep,
    face: FaceHandle,
    pcurve: &NurbsCurve,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::SurfaceEval as _;
    let period = brep
        .faces
        .get(face)?
        .surface
        .as_ref()?
        .geometry
        .parameter_periodicity()
        .0?;
    if !(period.is_finite() && period > 0.0) {
        return None;
    }
    let ring = crate::clip::face_uv_boundary(brep, face)?;
    if pcurve.control_points.is_empty() {
        return None;
    }
    let mut samples: Vec<(f64, f64)> = Vec::new();
    for (_, pc) in &ring {
        let (a, b) = pc.parameter_range();
        if !(a.is_finite() && b.is_finite()) {
            continue;
        }
        for i in 0..=4 {
            let q = pc.point_at(a + (b - a) * i as f64 / 4.0);
            samples.push((q.x, q.y));
        }
    }
    if samples.is_empty() {
        return None;
    }
    let global_center = samples.iter().map(|s| s.0).sum::<f64>() / samples.len() as f64;
    let (u_lo, u_hi) = samples.iter().fold((f64::MAX, f64::MIN), |(lo, hi), s| {
        (lo.min(s.0), hi.max(s.0))
    });
    let multi_turn = u_hi - u_lo > period * (1.0 + 1e-9);
    // MULTI-WRAP bands (PK-481): the global mean lands mid-band
    // regardless of which wrap the record crosses (the PK-464 trap).
    // Use the boundary samples near the QUERY's v range — the wrap
    // the record actually addresses — mirroring the internal
    // normalize's multi-turn provision.
    let boundary_center = if multi_turn {
        let (vq_lo, vq_hi) = pcurve
            .control_points
            .iter()
            .map(|p| p.y)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                (lo.min(v), hi.max(v))
            });
        let pad = 0.5 * (vq_hi - vq_lo).max(period * 0.05);
        let (sum, count) = samples
            .iter()
            .filter(|s| s.1 >= vq_lo - pad && s.1 <= vq_hi + pad)
            .fold((0.0_f64, 0usize), |(s, n), q| (s + q.0, n + 1));
        if count > 0 {
            sum / count as f64
        } else {
            global_center
        }
    } else {
        global_center
    };
    // Scope: negative-branch charts (the rem_euclid fold can never
    // reach them) and multi-wrap bands (no single fold was ever
    // right). Principal-range single-period charts stay with the
    // tuned machinery.
    if boundary_center >= -period * 0.25 && !multi_turn {
        return None;
    }
    let query_mean =
        pcurve.control_points.iter().map(|p| p.x).sum::<f64>() / pcurve.control_points.len() as f64;
    let k = ((boundary_center - query_mean) / period).round();
    if k == 0.0 {
        return None;
    }
    let mut translated = pcurve.clone();
    for p in &mut translated.control_points {
        p.x += k * period;
    }
    Some(translated)
}
