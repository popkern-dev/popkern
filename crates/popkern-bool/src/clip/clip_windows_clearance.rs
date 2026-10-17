//! Boundary-clearance sampling for the clip's hugging tests — split
//! from `clip_windows.rs` (500-line cap). See PK-114 (open records),
//! PK-393/PK-394 (closed flush rims) for the measured classes.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, FaceHandle};

/// One UV sample run per boundary contour (outer loop first, then
/// each hole).
pub(super) type ContourRuns = Vec<Vec<(f64, f64)>>;

/// Sampled face boundary in UV plus the clearance band (2 % of the
/// trim scale). One sample RUN per contour — the outer loop and each
/// inner loop (hole) separately, so the distance test never sees a
/// phantom segment connecting unrelated contours; a curve hugging a
/// HOLE rim is hugging the face boundary exactly like one hugging
/// the outer rim.
pub(super) fn boundary_clearance_samples(
    brep: &BRep,
    face: FaceHandle,
) -> Option<(ContourRuns, f64)> {
    crate::clip::face_uv_boundary(brep, face).map(|boundary| {
        // OUTER contour only: the hugging guards are tuned against
        // outer-boundary flush classes; adding hole rims here changed
        // near-hole verdicts on tuned socket fixtures.
        let contours = vec![boundary];
        let mut runs = Vec::with_capacity(contours.len());
        let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for contour in &contours {
            let mut samples = Vec::new();
            for (_, pcurve) in contour {
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
            runs.push(samples);
        }
        let trim_scale = ((hi_u - lo_u).powi(2) + (hi_v - lo_v).powi(2))
            .sqrt()
            .max(1e-9);
        (runs, trim_scale * 0.02)
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

/// Distance from `p` to the boundary POLYLINES (segments, not
/// samples; 25 samples over a 2π chart sit ~0.26 apart while the
/// clearance is ~0.13), with the query folded onto each segment's
/// periodic branch. Each contour is its own run — segments never
/// bridge contours.
pub(super) fn near_boundary(
    runs: &[Vec<(f64, f64)>],
    clearance: f64,
    p: (f64, f64),
    periods: (Option<f64>, Option<f64>),
) -> bool {
    runs.iter()
        .any(|samples| near_boundary_run(samples, clearance, p, periods))
}

fn near_boundary_run(
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

/// Whether `curve_uv` lies along the face's boundary contour for its
/// whole length: every one of 32 samples within the boundary band
/// ([`super::clip_uv::uv_boundary_band`]) of a boundary pcurve's
/// polyline (segments, not sample points: a band finer than the
/// sample spacing read a record exactly on a rim as off it).
pub(super) fn rides_face_boundary(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &popkern_core::geometry::nurbs_curve::NurbsCurve,
) -> bool {
    let band = super::clip_uv::uv_boundary_band(brep, curve_uv);
    let Some(ring) = crate::clip::face_uv_boundary(brep, face) else {
        return false;
    };
    let mut runs: Vec<Vec<(f64, f64)>> = Vec::new();
    for (_, pc) in &ring {
        let (a, b) = pc.parameter_range();
        if !(a.is_finite() && b.is_finite()) {
            return false;
        }
        runs.push(
            (0..=256)
                .map(|i| {
                    let p = pc.point_at(a + (b - a) * f64::from(i) / 256.0);
                    (p.x, p.y)
                })
                .collect(),
        );
    }
    let to_segment = |q: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let l2 = dx * dx + dy * dy;
        let t = if l2 > 0.0 {
            (((q.0 - a.0) * dx + (q.1 - a.1) * dy) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        ((q.0 - a.0 - t * dx).powi(2) + (q.1 - a.1 - t * dy).powi(2)).sqrt()
    };
    let (t0, t1) = curve_uv.parameter_range();
    (0..32).all(|k| {
        let p = curve_uv.point_at(t0 + (t1 - t0) * (f64::from(k) + 0.5) / 32.0);
        runs.iter()
            .flat_map(|run| run.windows(2))
            .any(|w| to_segment((p.x, p.y), w[0], w[1]) < band)
    })
}
