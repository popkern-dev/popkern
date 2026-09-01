//! Boundary-clearance sampling for the clip's hugging tests — split
//! from `clip_windows.rs` (500-line cap). See PK-114 (open records),
//! PK-393/PK-394 (closed flush rims) for the measured classes.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, FaceHandle};

/// Sampled face boundary in UV plus the clearance band (2 % of the
/// trim scale).
pub(super) fn boundary_clearance_samples(
    brep: &BRep,
    face: FaceHandle,
) -> Option<(Vec<(f64, f64)>, f64)> {
    crate::clip::face_uv_boundary(brep, face).map(|boundary| {
        let mut samples = Vec::new();
        let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for (_, pcurve) in &boundary {
            let curve = pcurve.to_nurbs();
            let (b0, b1) = curve.parameter_range();
            for k in 0..=24 {
                let point = curve.point_at(b0 + (b1 - b0) * k as f64 / 24.0);
                lo_u = lo_u.min(point.x);
                hi_u = hi_u.max(point.x);
                lo_v = lo_v.min(point.y);
                hi_v = hi_v.max(point.y);
                samples.push((point.x, point.y));
            }
        }
        let trim_scale = ((hi_u - lo_u).powi(2) + (hi_v - lo_v).powi(2))
            .sqrt()
            .max(1e-9);
        (samples, trim_scale * 0.02)
    })
}

/// Signed periodic fold of a coordinate difference into ±period/2.
pub(super) fn fold_signed(d: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) if p > 0.0 => {
            let mut r = d.rem_euclid(p);
            if r > p / 2.0 {
                r -= p;
            }
            r
        }
        _ => d,
    }
}

/// Distance from `p` to the boundary POLYLINE (segments, not samples;
/// 25 samples over a 2π chart sit ~0.26 apart while the clearance is
/// ~0.13), with the query folded onto each segment's periodic branch.
pub(super) fn near_boundary(
    samples: &[(f64, f64)],
    clearance: f64,
    p: (f64, f64),
    periods: (Option<f64>, Option<f64>),
) -> bool {
    samples.windows(2).any(|seg| {
        let (au, av) = seg[0];
        let (bu, bv) = seg[1];
        let qu = au + fold_signed(p.0 - au, periods.0);
        let qv = av + fold_signed(p.1 - av, periods.1);
        let (du, dv) = (bu - au, bv - av);
        let len2 = du * du + dv * dv;
        let s = if len2 > 0.0 {
            (((qu - au) * du + (qv - av) * dv) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let (cu, cv) = (au + s * du, av + s * dv);
        ((qu - cu).powi(2) + (qv - cv).powi(2)).sqrt() < clearance
    })
}
