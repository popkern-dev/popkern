//! Per-edge two-sided seam evaluation shared by
//! [`crate::cap_seam_continuity`] and
//! [`crate::blend_continuity_evidence`] — the pcurve-true sampling
//! core (orientation resolved from endpoints, degenerate charts
//! skipped, convergence-guarded curvature). Split out of
//! `cap_seam_continuity.rs` for the 500-line structural cap; the
//! measurement semantics and their rationale live on the items.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::cap_seam_continuity::SEAM_END_MARGIN;

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
) -> Option<EdgeSeam<'a>> {
    let curve = edge.curve.as_ref()?;
    let side_a = seam_side(brep, adjacency.half_edge, curve)?;
    let side_b = seam_side(brep, adjacency.twin, curve)?;
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
}

/// Dense reference-sampling resolution along a seam, for
/// nearest-point trim scoring and side-B localization.
const SEAM_REFERENCE_SAMPLES: usize = 64;

/// Sample one seam at `samples_per_seam` fractions (inset by
/// [`SEAM_END_MARGIN`], degenerate charts skipped — see module docs).
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
pub(crate) fn measure_edge_seam(seam: &EdgeSeam<'_>, samples_per_seam: usize) -> SeamWorst {
    let mut out = SeamWorst {
        g0: 0.0,
        g1: 0.0,
        g1_fraction: None,
        angle_deg: 0.0,
        g2: 0.0,
        g2_transverse: 0.0,
        samples: 0,
    };
    if samples_per_seam < 2 {
        return out;
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
    // True nearest-point distance: coarse nearest reference sample,
    // then golden-section refinement (the coarse grid alone leaves
    // ~half a segment of false gap — measured 6.3e-3 on an exact
    // cube seam against the 1e-5 G0 gate).
    let dist_to_curve = |p: &Point3| -> f64 {
        let coarse = reference
            .iter()
            .enumerate()
            .min_by(|x, y| {
                p.distance_to(x.1)
                    .partial_cmp(&p.distance_to(y.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i as f64 / SEAM_REFERENCE_SAMPLES as f64)
            .unwrap_or(0.5);
        let window = 1.0 / SEAM_REFERENCE_SAMPLES as f64;
        let (mut lo, mut hi) = (
            (coarse - window).clamp(0.0, 1.0),
            (coarse + window).clamp(0.0, 1.0),
        );
        for _ in 0..48 {
            let m1 = lo + (hi - lo) * 0.382;
            let m2 = lo + (hi - lo) * 0.618;
            if p.distance_to(&curve_at(m1)) < p.distance_to(&curve_at(m2)) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        p.distance_to(&curve_at(0.5 * (lo + hi)))
    };
    // Coarse side-B image cache for the nearest-point search.
    let b_images: Vec<(f64, Point3)> = (0..=SEAM_REFERENCE_SAMPLES)
        .map(|r| {
            let g = r as f64 / SEAM_REFERENCE_SAMPLES as f64;
            (g, seam.side_b.eval_fraction(g).1)
        })
        .collect();

    for k in 0..samples_per_seam {
        let f = SEAM_END_MARGIN
            + (1.0 - 2.0 * SEAM_END_MARGIN) * k as f64 / (samples_per_seam - 1) as f64;
        let (uv_a, p_a) = seam.side_a.eval_fraction(f);
        // Locate side B at the same 3D position: coarse nearest,
        // then golden-section refinement in the bracketing window.
        let coarse = b_images
            .iter()
            .min_by(|x, y| {
                p_a.distance_to(&x.1)
                    .partial_cmp(&p_a.distance_to(&y.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(g, _)| *g)
            .unwrap_or(f);
        let window = 1.0 / SEAM_REFERENCE_SAMPLES as f64;
        let (mut lo, mut hi) = (
            (coarse - window).clamp(0.0, 1.0),
            (coarse + window).clamp(0.0, 1.0),
        );
        for _ in 0..48 {
            let m1 = lo + (hi - lo) * 0.382;
            let m2 = lo + (hi - lo) * 0.618;
            if p_a.distance_to(&seam.side_b.eval_fraction(m1).1)
                < p_a.distance_to(&seam.side_b.eval_fraction(m2).1)
            {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        let g = 0.5 * (lo + hi);
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
        out.g0 = out.g0.max(d_a.max(d_b));

        let (Some(n_a), Some(n_b)) = (seam.side_a.normal(uv_a), seam.side_b.normal(uv_b)) else {
            continue;
        };
        out.samples += 1;
        let g1 = n_a.cross(&n_b).length();
        if g1 > out.g1 {
            out.g1 = g1;
            out.g1_fraction = Some(f);
        }
        let cos = n_a.dot(&n_b).abs().clamp(0.0, 1.0);
        out.angle_deg = out.angle_deg.max(cos.acos().to_degrees());

        // Seam direction from side A's own image (no curve-parameter
        // correspondence needed).
        let df = 1e-3;
        let p_fwd = seam.side_a.eval_fraction((f + df).min(1.0)).1;
        let p_bwd = seam.side_a.eval_fraction((f - df).max(0.0)).1;
        let Some(seam_dir) = (p_fwd - p_bwd).normalized() else {
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
            Some((k_a - k_b).abs())
        };
        if let Some(jump) = kappa_jump(seam_dir) {
            out.g2 = out.g2.max(jump);
        }
        // Transverse + diagonal directions carry the actual
        // curvature-continuity content — see the SeamWorst docs.
        let n_avg = ((n_a + n_b_aligned) * 0.5).normalized();
        if let Some(n_avg) = n_avg
            && let Some(t_perp) = n_avg.cross(&seam_dir).normalized()
        {
            if let Some(jump) = kappa_jump(t_perp) {
                out.g2_transverse = out.g2_transverse.max(jump);
            }
            if let Some(diag) = (seam_dir + t_perp).normalized()
                && let Some(jump) = kappa_jump(diag)
            {
                out.g2_transverse = out.g2_transverse.max(jump);
            }
        }
    }
    out
}

pub(crate) struct SeamSide<'a> {
    surface: &'a SurfaceGeometry,
    pcurve: &'a Pcurve,
    /// Pcurve runs opposite to the edge curve (twin half-edges trace
    /// the loop direction, not the curve direction) — resolved once
    /// per side by matching pcurve-endpoint surface images to the
    /// curve endpoints.
    reversed: bool,
}

fn seam_side<'a>(brep: &'a BRep, he: HalfEdgeHandle, curve: &CurveOnEdge) -> Option<SeamSide<'a>> {
    let he_data = brep.half_edges.get(he)?;
    let pcurve = he_data.pcurve.as_ref()?;
    let face = brep.loops.get(he_data.loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let curve_start = (&curve.geometry as &dyn CurveEval).point_at(curve.t_start);
    let uv0 = pcurve.point_at(pcurve.t_start);
    let uv1 = pcurve.point_at(pcurve.t_end);
    let s0 = surface.point_at(uv0.x, uv0.y);
    let s1 = surface.point_at(uv1.x, uv1.y);
    let reversed = s1.distance_to(&curve_start) < s0.distance_to(&curve_start);
    Some(SeamSide {
        surface,
        pcurve,
        reversed,
    })
}

impl SeamSide<'_> {
    /// Evaluate at fraction `f` of the seam in CURVE direction.
    pub(crate) fn eval_fraction(&self, f: f64) -> ((f64, f64), Point3) {
        let g = if self.reversed { 1.0 - f } else { f };
        let uv = self
            .pcurve
            .point_at(self.pcurve.t_start + g * (self.pcurve.t_end - self.pcurve.t_start));
        ((uv.x, uv.y), self.surface.point_at(uv.x, uv.y))
    }

    /// Finite-difference unit normal, or `None` where the chart is
    /// too degenerate to trust (see module docs).
    pub(crate) fn normal(&self, uv: (f64, f64)) -> Option<Vec3> {
        let h = 1e-6;
        let du = self.surface.point_at(uv.0 + h, uv.1) - self.surface.point_at(uv.0 - h, uv.1);
        let dv = self.surface.point_at(uv.0, uv.1 + h) - self.surface.point_at(uv.0, uv.1 - h);
        let (lu, lv) = (du.length(), dv.length());
        if lu.min(lv) < 1e-3 * lu.max(lv) {
            return None;
        }
        du.cross(&dv).normalized()
    }

    /// Normal curvature along `dir`, signed against the caller's
    /// `normal` (hemisphere-aligned across a seam so signed
    /// curvatures compare — a chart-orientation flip must not read a
    /// perfect seam as a `2κ` jump), with a two-step convergence
    /// guard: near a degenerate polar chart the FD truncation error
    /// scales quadratically with the step (measured 1.7 → 0.16 →
    /// 1.9e-3 across steps 1e-3/1e-4/1e-5 on the same pentagon fan
    /// seam), so a single-step estimate reads chart noise as a
    /// curvature crease. Samples whose two-step estimates disagree
    /// are skipped; agreeing ones are Richardson-extrapolated.
    fn normal_curvature_signed(&self, uv: (f64, f64), dir: Vec3, normal: Vec3) -> Option<f64> {
        let k_h = self.normal_curvature_at_step(uv, dir, normal, 1e-4)?;
        let k_h2 = self.normal_curvature_at_step(uv, dir, normal, 5e-5)?;
        if (k_h - k_h2).abs() > 0.1 * (1.0 + k_h2.abs()) {
            return None;
        }
        Some((4.0 * k_h2 - k_h) / 3.0)
    }

    fn normal_curvature_at_step(
        &self,
        uv: (f64, f64),
        dir: Vec3,
        normal: Vec3,
        step: f64,
    ) -> Option<f64> {
        let h = 1e-6;
        let du = self.surface.point_at(uv.0 + h, uv.1) - self.surface.point_at(uv.0 - h, uv.1);
        let dv = self.surface.point_at(uv.0, uv.1 + h) - self.surface.point_at(uv.0, uv.1 - h);
        let su = du * (0.5 / h);
        let sv = dv * (0.5 / h);
        let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let (r1, r2) = (dir.dot(&su), dir.dot(&sv));
        let det = e * g - f_ * f_;
        if det.abs() < 1e-30 {
            return None;
        }
        let a = (r1 * g - r2 * f_) / det;
        let b = (r2 * e - r1 * f_) / det;
        let p0 = self.surface.point_at(uv.0, uv.1);
        let p1 = self.surface.point_at(uv.0 + a * step, uv.1 + b * step);
        let p2 = self.surface.point_at(uv.0 - a * step, uv.1 - b * step);
        let second = (p1 - p0) + (p2 - p0);
        let first_len = (p1 - p2).length() * 0.5;
        if first_len < 1e-30 {
            return None;
        }
        Some(second.dot(&normal) / (first_len * first_len))
    }
}
