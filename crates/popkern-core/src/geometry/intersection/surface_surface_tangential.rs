//! Crossing certificate for demoted-tangent intersection traces.
//!
//! A tangent-classified seed that marches (the grazing-crossing rescue
//! in `intersect_surfaces_with_handler`) traces either a genuine
//! crossing whose normals merely graze, or a near-contact blob where
//! two touching surfaces stay within tolerance around a corner or
//! boundary. The two are indistinguishable at the seed; they differ in
//! whether the surfaces change SIDE across the traced curve. This
//! module certifies the side change with the analytic operand's O(1)
//! implicit signed distance.

use crate::geometry::analytic_surfaces::analytic_signed_distance;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// True when probes offset to either side of the traced curve land on
/// opposite sides of the analytic operand (implicit signed distance
/// changes sign with margin). Returns false — keeping the pair on the
/// tangency-silent contract — when neither operand is analytic, when
/// the offsets cannot be constructed, or when no probe shows a
/// definite sign change (the near-contact-blob signature).
pub(super) fn tangential_crossing_certificate(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    samples: &[MarchSample],
    tolerance: &Tolerance,
) -> bool {
    if samples.len() < 3 {
        return false;
    }
    // Pick the implicit (analytic) side and the offset (free-form or
    // second analytic) side. Probing offsets on the non-implicit side
    // keeps the probe points ON one surface while the other provides
    // the side measure.
    let probe_point = samples[samples.len() / 2].point;
    let (implicit, offset_is_surface1) =
        if analytic_signed_distance(surface2, &probe_point).is_some() {
            (surface2, true)
        } else if analytic_signed_distance(surface1, &probe_point).is_some() {
            (surface1, false)
        } else {
            return false;
        };
    let offset_surface = if offset_is_surface1 {
        surface1
    } else {
        surface2
    };

    let linear = tolerance.linear();
    let total_len: f64 = samples
        .windows(2)
        .map(|w| w[0].point.distance_to(&w[1].point))
        .sum();

    for idx in [samples.len() / 2, samples.len() / 4, 3 * samples.len() / 4] {
        if idx == 0 || idx + 1 >= samples.len() {
            continue;
        }
        let s = &samples[idx];
        let (u, v) = if offset_is_surface1 {
            (s.u1, s.v1)
        } else {
            (s.u2, s.v2)
        };
        let su = offset_surface.partial_u(u, v);
        let sv = offset_surface.partial_v(u, v);
        let n_ff = su.cross(&sv);
        let t3d = samples[idx + 1].point - samples[idx - 1].point;
        let Some(across) = t3d.cross(&n_ff).normalized() else {
            continue;
        };

        // Offset distance: clear of the |S1-S2| < ε contact band,
        // whose half-width across a crossing at normal angle θ is
        // ~ε/sin θ; capped by the trace scale so probes stay local.
        let (iu, iv) = if offset_is_surface1 {
            (s.u2, s.v2)
        } else {
            (s.u1, s.v1)
        };
        let n_imp = implicit.normal_at(iu, iv);
        let norms = n_ff.length() * n_imp.length();
        if norms < tolerance.numerical_floor() {
            continue;
        }
        let sin_angle = (n_ff.cross(&n_imp).length() / norms).max(1.0e-6);
        let h = (20.0 * linear / sin_angle)
            .min(0.25 * total_len)
            .max(20.0 * linear);

        // Map the 3-D across offset into the offset surface's UV via
        // least squares on the first partials.
        let g00 = su.dot(&su);
        let g01 = su.dot(&sv);
        let g11 = sv.dot(&sv);
        let det = g00 * g11 - g01 * g01;
        if det.abs() < tolerance.numerical_floor() {
            continue;
        }
        let target = across * h;
        let b0 = su.dot(&target);
        let b1 = sv.dot(&target);
        let du = (b0 * g11 - b1 * g01) / det;
        let dv = (b1 * g00 - b0 * g01) / det;

        let p_pos = offset_surface.point_at(u + du, v + dv);
        let p_neg = offset_surface.point_at(u - du, v - dv);
        let (Some(d_pos), Some(d_neg)) = (
            analytic_signed_distance(implicit, &p_pos),
            analytic_signed_distance(implicit, &p_neg),
        ) else {
            continue;
        };
        if d_pos * d_neg < 0.0 && d_pos.abs().min(d_neg.abs()) > 2.0 * linear {
            return true;
        }
    }
    false
}
