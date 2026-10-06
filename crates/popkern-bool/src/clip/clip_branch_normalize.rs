//! Period-branch normalization for the UV clip query (split from
//! `clip_uv_crossings.rs` under the 500-line file cap).

use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

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
            // SCORED STRIPE PICK: the control-point MEAN near the
            // query's v is off by a whole period when the two long
            // boundary helices pass that v at different unwound u
            // and their control points are sparser than the pad
            // (measured, PK-481: the chamfer cone's rim circle at
            // v = 6.0 on the bore wall's 4.6-period unwound trim —
            // the mean landed the query at u ∈ [−12.6, −6.3] while
            // the band's v = 6.0 stripe sits at ≈ [−14.6, −12.9];
            // every crossing and containment sample missed and the
            // closed record dropped whole, stranding the chamfer
            // chain). Sample each boundary pcurve finely and pick
            // the stripe (integer shift) that actually CONTAINS
            // boundary material at the query's v; the mean stays
            // as the no-material fallback.
            let mut samples: Vec<(f64, f64)> = Vec::new();
            for (_, pc) in boundary_nurbs {
                use popkern_core::geometry::curve::CurveEval as _;
                let (t0, t1) = pc.parameter_range();
                for i in 0..=16 {
                    let t = t0 + (t1 - t0) * i as f64 / 16.0;
                    let p = pc.point_at(t);
                    samples.push((p.x, p.y));
                }
            }
            let (q_lo, q_hi) = normalized
                .control_points
                .iter()
                .map(|p| p.x)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
                    (lo.min(u), hi.max(u))
                });
            let (b_lo, b_hi) = samples
                .iter()
                .map(|(u, _)| *u)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
                    (lo.min(u), hi.max(u))
                });
            let k_min = ((b_lo - q_hi) / period).floor() as i64;
            let k_max = ((b_hi - q_lo) / period).ceil() as i64;
            let mut best: Option<(usize, f64)> = None;
            for k in k_min..=k_max {
                let (w_lo, w_hi) = (q_lo + k as f64 * period, q_hi + k as f64 * period);
                let score = samples
                    .iter()
                    .filter(|(u, v)| {
                        *u >= w_lo && *u <= w_hi && *v >= vq_lo - pad && *v <= vq_hi + pad
                    })
                    .count();
                let center = 0.5 * (w_lo + w_hi);
                let better = match best {
                    None => score > 0,
                    Some((bs, _)) => score > bs,
                };
                if better {
                    best = Some((score, center));
                }
            }
            if let Some((_, center)) = best {
                center
            } else {
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
                // Single-period charts keep the tuned per-point fold
                // for SUB-period queries: a query already on the
                // principal branch stays put even when k≠0 (measured:
                // the filleted-torus boolean canary stranded under an
                // honored k-translation), and the tear a fold puts in
                // a seam-straddling query is healed by the ±π seam
                // dance in the crossing search (measured: a rigid
                // landing of a query at u ∈ [−0.45, 0.2] broke the
                // filleted-operand subtraction, loop gap 0.87). A
                // FULL-period query is the exception: every fold of it
                // tears and no ±π shift can heal it (measured, PK-681:
                // the M2 nut's mouth-chamfer circle on the crest band,
                // u ∈ [2π, 4π], folded its closing control point from
                // 2π to 0; the circle ran backwards, no crossing with
                // the band's boundary was found, the record dropped,
                // and the uncut band classified Inside and vanished).
                // It lands rigidly.
                let (q_lo, q_hi) = normalized
                    .control_points
                    .iter()
                    .map(|p| p.x)
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
                        (lo.min(u), hi.max(u))
                    });
                normalized = if q_hi - q_lo >= period * (1.0 - 1e-6) {
                    translate_onto_principal_branch(&normalized, period)
                } else {
                    super::clip_seam_shift::shift_curve_u(&normalized, k * period)
                };
            }
            if popkern_core::trace::trace_enabled("clip") {
                eprintln!("    [uvclip-branch] {face:?} query curve shifted by {k:+.0} period(s)");
            }
        }
    }

    normalized
}

/// Translate `curve` by the whole number of periods that puts the
/// most of its u-extent on the principal branch `[0, period]`: the
/// principal branch, reached rigidly. A curve already mostly on it
/// stays put (a seam-anchored circle whose lowest control point dips
/// a hair below 0 keeps its place), and a curve ending exactly on the
/// seam keeps its end at `period` instead of folding it to 0. Ties go
/// to the smaller shift.
pub(super) fn translate_onto_principal_branch(curve: &NurbsCurve, period: f64) -> NurbsCurve {
    let (lo, hi) = curve
        .control_points
        .iter()
        .map(|p| p.x)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), u| {
            (lo.min(u), hi.max(u))
        });
    if !(lo.is_finite() && hi.is_finite()) || period <= 0.0 {
        return curve.clone();
    }
    let overlap = |k: f64| (hi + k * period).min(period) - (lo + k * period).max(0.0);
    let k_lo = (-hi / period).floor();
    let k_hi = (-lo / period).ceil();
    let mut best: (f64, f64) = (0.0, overlap(0.0));
    let mut k = k_lo;
    while k <= k_hi {
        let o = overlap(k);
        if o > best.1 + 1e-12 * period || (o > best.1 - 1e-12 * period && k.abs() < best.0.abs()) {
            best = (k, o);
        }
        k += 1.0;
    }
    let mut translated = curve.clone();
    if best.0 != 0.0 {
        for p in &mut translated.control_points {
            p.x += best.0 * period;
        }
    }
    translated
}

#[cfg(test)]
#[path = "clip_branch_normalize_tests.rs"]
mod tests;
