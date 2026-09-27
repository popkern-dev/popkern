//! Per-edge two-sided seam evaluation shared by
//! [`crate::cap_seam_continuity`] and
//! [`crate::blend_continuity_evidence`] — the pcurve-true sampling
//! core (orientation resolved from endpoints, degenerate charts
//! skipped, convergence-guarded curvature). Split out of
//! `cap_seam_continuity.rs` for the 500-line structural cap; the
//! measurement semantics and their rationale live on the items.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::cap_seam_side::{SeamSide, SurfaceJetCache, seam_side};

/// The two distinct faces adjacent to an edge, resolved defensively
/// (a mid-pipeline BRep can carry stale handles; a measurement must
/// report "not a seam", never panic).
pub(crate) struct EdgeSeamFaces {
    pub(crate) half_edge: HalfEdgeHandle,
    pub(crate) twin: HalfEdgeHandle,
    pub(crate) face_a: FaceHandle,
    pub(crate) face_b: FaceHandle,
}

/// Resolve an edge's adjacency, or `None` for boundary/seam-on-one-
/// face/stale edges.
pub(crate) fn edge_seam_faces(
    brep: &BRep,
    edge: &popkern_core::topology::Edge,
) -> Option<EdgeSeamFaces> {
    let he = edge.half_edge;
    let he_data = brep.half_edges.get(he)?;
    let twin = he_data.twin;
    let twin_data = brep.half_edges.get(twin)?;
    let face_a = brep.loops.get(he_data.loop_)?.face;
    let face_b = brep.loops.get(twin_data.loop_)?.face;
    (face_a != face_b).then_some(EdgeSeamFaces {
        half_edge: he,
        twin,
        face_a,
        face_b,
    })
}

/// Both sides of a seam, ready to sample.
pub(crate) struct EdgeSeam<'a> {
    pub(crate) side_a: SeamSide<'a>,
    pub(crate) side_b: SeamSide<'a>,
    curve: &'a CurveOnEdge,
}

/// Build the two evaluable sides for an edge whose adjacency is
/// already resolved. `None` when either side lacks the pcurve or
/// surface the measurement needs — callers that care about coverage
/// count that case separately instead of silently passing it.
pub(crate) fn resolve_edge_seam<'a>(
    brep: &'a BRep,
    edge: &'a popkern_core::topology::Edge,
    adjacency: &EdgeSeamFaces,
    jets: &'a SurfaceJetCache,
) -> Option<EdgeSeam<'a>> {
    let curve = edge.curve.as_ref()?;
    let side_a = seam_side(brep, adjacency.half_edge, curve, jets)?;
    let side_b = seam_side(brep, adjacency.twin, curve, jets)?;
    Some(EdgeSeam {
        side_a,
        side_b,
        curve,
    })
}

/// Worst-case readings for one seam.
pub(crate) struct SeamWorst {
    /// Worst positional deviation of either side from the edge curve.
    pub(crate) g0: f64,
    /// Worst `|n_a × n_b|`.
    pub(crate) g1: f64,
    /// Fraction along the seam where [`Self::g1`] occurred.
    pub(crate) g1_fraction: Option<f64>,
    /// Worst normal angle in degrees.
    pub(crate) angle_deg: f64,
    /// Worst normal-curvature jump ALONG the seam direction.
    ///
    /// **Near-vacuous for a G1-clean seam.** By Meusnier's theorem
    /// the normal curvature along a shared direction is `c″ · n`:
    /// both sides share the edge curve (G0) and, at G1, the normal —
    /// so along-seam agreement is *implied* by G0+G1 and this number
    /// is a second-order echo of [`Self::g1`], not independent
    /// curvature evidence. Kept because the historical gates pinned
    /// it; honest curvature continuity is [`Self::g2_transverse`].
    pub(crate) g2: f64,
    /// Worst normal-curvature jump across the TRANSVERSE and
    /// diagonal in-tangent-plane directions — the directions that
    /// actually distinguish curvature continuity (a rolling-ball
    /// rail jumps exactly `1/r` here while reading 0 along the
    /// seam). Together with the along-seam term this pins the full
    /// second fundamental form on the shared tangent plane.
    pub(crate) g2_transverse: f64,
    /// Samples that produced a usable normal on BOTH sides. Zero
    /// means the seam contributed positional evidence only — G1/G2
    /// were not measured and must not be read as passing.
    pub(crate) samples: usize,
    /// Samples with normals on both sides whose transverse or
    /// diagonal curvature reading was skipped (a singular chart point:
    /// the first fundamental form degenerates). Counted, never
    /// silently dropped (PK-579).
    pub(crate) skipped_curvature: usize,
}

/// Dense reference-sampling resolution along a seam, for
/// nearest-point trim scoring and side-B localization.
const SEAM_REFERENCE_SAMPLES: usize = 64;

/// Worst readings over [`sample_edge_seam`]'s samples: the evidence
/// instrument's per-seam summary.
pub(crate) fn measure_edge_seam(
    seam: &EdgeSeam<'_>,
    samples_per_seam: usize,
    end_margin: f64,
) -> SeamWorst {
    let mut out = SeamWorst {
        g0: 0.0,
        g1: 0.0,
        g1_fraction: None,
        angle_deg: 0.0,
        g2: 0.0,
        g2_transverse: 0.0,
        samples: 0,
        skipped_curvature: 0,
    };
    for r in sample_edge_seam(seam, samples_per_seam, end_margin) {
        out.g0 = out.g0.max(r.g0);
        let Some(g1) = r.g1 else { continue };
        out.samples += 1;
        if g1 > out.g1 {
            out.g1 = g1;
            out.g1_fraction = Some(r.fraction);
        }
        out.angle_deg = out.angle_deg.max(r.angle_deg);
        if let Some(j) = r.jump_along {
            out.g2 = out.g2.max(j.abs());
        }
        if r.jump_perp.is_none() || r.jump_diag.is_none() {
            out.skipped_curvature += 1;
        }
        for j in [r.jump_perp, r.jump_diag].into_iter().flatten() {
            out.g2_transverse = out.g2_transverse.max(j.abs());
        }
        // The largest normal-curvature jump over EVERY tangent
        // direction: along, transverse and one diagonal determine the
        // 2×2 jump form, and its spectral norm is the honest G2
        // residual. The three sampled directions alone under-read a
        // jump that peaks on the other diagonal (pentagon5_variable
        // cap-cap: 3.1e-2 sampled vs 7.0e-2 form, certified >= 5.0e-2).
        if let (Some(a), Some(c), Some(d)) = (r.jump_along, r.jump_perp, r.jump_diag) {
            let (form, _) = crate::seam_profile::second_form_jump(a, c, d);
            out.g2_transverse = out.g2_transverse.max(form);
        }
    }
    out
}

/// Coarse reference candidates refined per nearest-point query.
const NEAREST_CANDIDATES: usize = 2;

/// Fraction in `[0, 1]` of the point nearest `p` on a curve sampled at
/// `reference` (evenly spaced fractions) whose exact distance is
/// `dist(g)`: golden-section refinement in the bracketing window of
/// each of the [`NEAREST_CANDIDATES`] nearest reference samples, the
/// best refined result kept.
///
/// One candidate is not enough on a CLOSED curve: near the closure the
/// first and last reference samples are equidistant, the tie keeps the
/// first, and a point just before the end was refined at the START,
/// reading the chord to the start as a gap (measured: 1.7e-2 on the
/// exact boss-plate rim circle once samples reached within 1% of its
/// closure, PK-579).
fn nearest_fraction(reference: &[Point3], dist: impl Fn(f64) -> f64, p: &Point3) -> f64 {
    let last = reference.len().saturating_sub(1).max(1);
    let mut order: Vec<usize> = (0..reference.len()).collect();
    order.sort_by(|&x, &y| {
        p.distance_to(&reference[x])
            .total_cmp(&p.distance_to(&reference[y]))
    });
    let window = 1.0 / last as f64;
    let mut best = (f64::INFINITY, 0.5);
    for &i in order.iter().take(NEAREST_CANDIDATES) {
        let coarse = i as f64 / last as f64;
        let (mut lo, mut hi) = (
            (coarse - window).clamp(0.0, 1.0),
            (coarse + window).clamp(0.0, 1.0),
        );
        for _ in 0..48 {
            let m1 = lo + (hi - lo) * 0.382;
            let m2 = lo + (hi - lo) * 0.618;
            if dist(m1) < dist(m2) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        let g = 0.5 * (lo + hi);
        let d = dist(g);
        if d < best.0 {
            best = (d, g);
        }
    }
    best.1
}

/// One sample of a seam: where it is, its positional and normal
/// readings, and the SIGNED normal-curvature jumps `κ_a − κ_b` along,
/// across (transverse) and diagonally (`along + transverse`, unit)
/// in the shared tangent plane.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SeamSampleReading {
    /// Fraction along the seam (curve direction).
    pub(crate) fraction: f64,
    /// Side A's point.
    pub(crate) point: Point3,
    /// Worse side's distance to the edge curve.
    pub(crate) g0: f64,
    /// `|n_a × n_b|`, when both sides give a normal.
    pub(crate) g1: Option<f64>,
    /// Normal angle in degrees (0 without normals).
    pub(crate) angle_deg: f64,
    /// Signed along-seam curvature jump.
    pub(crate) jump_along: Option<f64>,
    /// Signed transverse curvature jump.
    pub(crate) jump_perp: Option<f64>,
    /// Signed diagonal curvature jump.
    pub(crate) jump_diag: Option<f64>,
}

/// Sample one seam at `samples_per_seam` fractions inset by
/// `end_margin` at each end; degenerate charts skipped — see module
/// docs.
///
/// **Position-matched, not fraction-matched.** Side A is walked by
/// pcurve fraction; side B is located per sample by nearest-point
/// refinement to side A's 3D position. Same-fraction pairing breaks
/// on CLOSED seams (a boss-base ring fillet rail): the curve's
/// endpoints coincide, so endpoint orientation resolution is
/// degenerate and a phase-shifted pcurve gets compared against the
/// wrong point of the circle — measured: g0 = 14 on a watertight
/// housing. G0 is scored as each side's distance to the shared edge
/// curve (the trim contract validation and tessellation consume),
/// which is phase-independent.
pub(crate) fn sample_edge_seam(
    seam: &EdgeSeam<'_>,
    samples_per_seam: usize,
    end_margin: f64,
) -> Vec<SeamSampleReading> {
    let mut readings = Vec::with_capacity(samples_per_seam);
    if samples_per_seam < 2 {
        return readings;
    }
    let curve = seam.curve;
    // Dense curve reference for trim scoring.
    let reference: Vec<Point3> = (0..=SEAM_REFERENCE_SAMPLES)
        .map(|r| {
            let t = curve.t_start
                + r as f64 / SEAM_REFERENCE_SAMPLES as f64 * (curve.t_end - curve.t_start);
            (&curve.geometry as &dyn CurveEval).point_at(t)
        })
        .collect();
    let curve_at = |g: f64| -> Point3 {
        (&curve.geometry as &dyn CurveEval)
            .point_at(curve.t_start + g * (curve.t_end - curve.t_start))
    };
    // True nearest-point distance: golden-section refinement around
    // the best coarse reference samples (the coarse grid alone leaves
    // ~half a segment of false gap — measured 6.3e-3 on an exact cube
    // seam against the 1e-5 G0 gate). See [`nearest_fraction`] for why
    // more than one coarse candidate is refined.
    let dist_to_curve = |p: &Point3| -> f64 {
        let g = nearest_fraction(&reference, |g| p.distance_to(&curve_at(g)), p);
        p.distance_to(&curve_at(g))
    };
    // Coarse side-B image cache for the nearest-point search.
    let b_images: Vec<(f64, Point3)> = (0..=SEAM_REFERENCE_SAMPLES)
        .map(|r| {
            let g = r as f64 / SEAM_REFERENCE_SAMPLES as f64;
            (g, seam.side_b.eval_fraction(g).1)
        })
        .collect();

    for k in 0..samples_per_seam {
        let f = end_margin + (1.0 - 2.0 * end_margin) * k as f64 / (samples_per_seam - 1) as f64;
        let (uv_a, p_a) = seam.side_a.eval_fraction(f);
        // Locate side B at the same 3D position (coarse candidates,
        // then golden-section refinement — see [`nearest_fraction`]).
        let b_points: Vec<Point3> = b_images.iter().map(|(_, q)| *q).collect();
        let g = nearest_fraction(
            &b_points,
            |g| p_a.distance_to(&seam.side_b.eval_fraction(g).1),
            &p_a,
        );
        let (uv_b, p_b) = seam.side_b.eval_fraction(g);

        // G0: each side against the shared curve (trim contract).
        let (d_a, d_b) = (dist_to_curve(&p_a), dist_to_curve(&p_b));
        if popkern_core::trace::trace_enabled("seam-measure") {
            eprintln!(
                "  [seam-measure] f={f:.3} g={g:.3} a=({:.4},{:.4},{:.4}) d_a={d_a:.3e} b=({:.4},{:.4},{:.4}) d_b={d_b:.3e} curve t[{:.4},{:.4}] c0=({:.4},{:.4},{:.4}) c1=({:.4},{:.4},{:.4})",
                p_a.x,
                p_a.y,
                p_a.z,
                p_b.x,
                p_b.y,
                p_b.z,
                curve.t_start,
                curve.t_end,
                reference[0].x,
                reference[0].y,
                reference[0].z,
                reference[SEAM_REFERENCE_SAMPLES].x,
                reference[SEAM_REFERENCE_SAMPLES].y,
                reference[SEAM_REFERENCE_SAMPLES].z
            );
        }
        let mut reading = SeamSampleReading {
            fraction: f,
            point: p_a,
            g0: d_a.max(d_b),
            g1: None,
            angle_deg: 0.0,
            jump_along: None,
            jump_perp: None,
            jump_diag: None,
        };

        let (Some(n_a), Some(n_b)) = (seam.side_a.normal(uv_a), seam.side_b.normal(uv_b)) else {
            readings.push(reading);
            continue;
        };
        reading.g1 = Some(n_a.cross(&n_b).length());
        let cos = n_a.dot(&n_b).abs().clamp(0.0, 1.0);
        reading.angle_deg = cos.acos().to_degrees();

        // Seam direction from side A's own image (no curve-parameter
        // correspondence needed).
        let df = 1e-3;
        let p_fwd = seam.side_a.eval_fraction((f + df).min(1.0)).1;
        let p_bwd = seam.side_a.eval_fraction((f - df).max(0.0)).1;
        let Some(seam_dir) = (p_fwd - p_bwd).normalized() else {
            readings.push(reading);
            continue;
        };
        // Chart orientations are arbitrary: flip side B's normal into
        // side A's hemisphere before comparing signed curvatures, or
        // a reversed chart reads a perfect seam as a 2κ jump.
        let n_b_aligned = if n_a.dot(&n_b) < 0.0 { n_b * -1.0 } else { n_b };
        let kappa_jump = |dir: Vec3| -> Option<f64> {
            let k_a = seam.side_a.normal_curvature_signed(uv_a, dir, n_a)?;
            let k_b = seam
                .side_b
                .normal_curvature_signed(uv_b, dir, n_b_aligned)?;
            Some(k_a - k_b)
        };
        reading.jump_along = kappa_jump(seam_dir);
        // Transverse + diagonal directions carry the actual
        // curvature-continuity content — see the SeamWorst docs.
        let n_avg = ((n_a + n_b_aligned) * 0.5).normalized();
        if let Some(n_avg) = n_avg
            && let Some(t_perp) = n_avg.cross(&seam_dir).normalized()
        {
            let jump_perp = kappa_jump(t_perp);
            let jump_diag = (seam_dir + t_perp).normalized().and_then(&kappa_jump);
            reading.jump_perp = jump_perp;
            reading.jump_diag = jump_diag;
            if popkern_core::trace::trace_enabled("seam-curvature") {
                eprintln!(
                    "  [seam-curvature] f={f:.3} perp={} diag={}",
                    jump_perp.map_or("skip".into(), |j| format!("{:.3e}", j.abs())),
                    jump_diag.map_or("skip".into(), |j| format!("{:.3e}", j.abs())),
                );
            }
        }
        readings.push(reading);
    }
    readings
}
