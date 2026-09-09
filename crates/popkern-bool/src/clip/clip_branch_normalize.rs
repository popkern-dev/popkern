//! Period-branch normalization for the UV clip query (split from
//! `clip_uv_crossings.rs` under the 500-line file cap).

use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_seam_shift::shift_curve_u;

/// Shift the query curve by whole periods so it addresses the same
/// chart branch as the face's boundary pcurves. See the block
/// comments inline — the multi-turn/plain-translation split and the
/// single-period fold are both measured, load-bearing behaviors.
pub(super) fn normalize_query_branch(
    brep: &BRep,
    face: FaceHandle,
    boundary_nurbs: &[(popkern_core::topology::HalfEdgeHandle, NurbsCurve)],
    curve_uv: &NurbsCurve,
    surface_is_periodic: bool,
) -> NurbsCurve {
    // PERIOD-BRANCH NORMALIZATION: the query curve may live on any
    // whole-period u-branch of a periodic chart (a plane×cylinder
    // marcher parameterized one gusset-web ellipse at u ∈ [2π, 4π]),
    // while the face's pcurves live on the principal branch. UV
    // curve×curve intersection then finds no true crossings — and
    // the ±π seam dance below cannot bridge a full period (measured:
    // the arc-web×underside record clipped to the EXTENDED
    // cylinder's phantom re-crossing instead of the true weld line,
    // and the fanned arm's gusset union went non-manifold). Shift
    // the query curve by whole periods so its u-centre lands nearest
    // the boundary's u-centre; crossing t-parameters are unaffected.
    let mut normalized = curve_uv.clone();
    if surface_is_periodic && !boundary_nurbs.is_empty() {
        let period = {
            use popkern_core::geometry::SurfaceEval as _;
            brep.faces
                .get(face)
                .and_then(|f| f.surface.as_ref())
                .and_then(|s| s.geometry.parameter_periodicity().0)
                .unwrap_or(std::f64::consts::TAU)
        };
        let u_mean = |c: &NurbsCurve| {
            c.control_points.iter().map(|p| p.x).sum::<f64>() / c.control_points.len() as f64
        };
        let global_center = boundary_nurbs.iter().map(|(_, pc)| u_mean(pc)).sum::<f64>()
            / boundary_nurbs.len() as f64;
        // A MULTI-TURN chart (a helical thread band whose
        // boundary pcurves unroll several periods) has no single
        // meaningful u-centre: the global mean lands the query on
        // a mid-band turn regardless of which turn it addresses
        // (measured, PK-464: the top mouth-chamfer circle at
        // v = 4.83 shifted −2 periods toward the band's middle,
        // classified A=Outside, and the whole record dropped —
        // the ring never closed). There the centre comes from the
        // boundary control points near the QUERY's v range — the
        // turn the query actually crosses. Single-period charts
        // keep the global centre bit-identically.
        let boundary_u_extent = {
            let (lo, hi) = boundary_nurbs
                .iter()
                .flat_map(|(_, pc)| pc.control_points.iter().map(|p| p.x))
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
                    (lo.min(u), hi.max(u))
                });
            hi - lo
        };
        let boundary_center = if boundary_u_extent > period * (1.0 + 1e-9) {
            let (vq_lo, vq_hi) = normalized
                .control_points
                .iter()
                .map(|p| p.y)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
                    (lo.min(v), hi.max(v))
                });
            let pad = 0.5 * (vq_hi - vq_lo).max(period * 0.05);
            let (sum, count) = boundary_nurbs
                .iter()
                .flat_map(|(_, pc)| pc.control_points.iter())
                .filter(|p| p.y >= vq_lo - pad && p.y <= vq_hi + pad)
                .fold((0.0_f64, 0usize), |(s, n), p| (s + p.x, n + 1));
            if count > 0 {
                sum / count as f64
            } else {
                global_center
            }
        } else {
            global_center
        };
        let k = ((boundary_center - u_mean(&normalized)) / period).round();
        if k != 0.0 {
            if boundary_u_extent > period * (1.0 + 1e-9) {
                // MULTI-TURN charts: PLAIN translation.
                // `shift_curve_u` folds every control point by
                // rem_euclid(2π), which cannot land on a NEGATIVE
                // branch — the −2π shift toward the nut crest
                // band's unrolled chart folded the query to a
                // zero-span point at u = 0 and the record dropped
                // whole (PK-464, the top mouth-chamfer crest
                // arc). Multi-turn boundaries never had a working
                // fold to preserve.
                for p in &mut normalized.control_points {
                    p.x += k * period;
                }
            } else {
                // Single-period charts keep the tuned fold: for a
                // query already inside the principal range a k≠0
                // fold is a NO-OP, and tuned fixtures depend on
                // exactly that (measured: the filleted-torus
                // boolean canary stranded under an honored
                // translation).
                normalized = shift_curve_u(&normalized, k * period);
            }
            if popkern_core::trace::trace_enabled("clip") {
                eprintln!("    [uvclip-branch] {face:?} query curve shifted by {k:+.0} period(s)");
            }
        }
    }

    normalized
}
