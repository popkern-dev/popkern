//! Continuity across a cap's own internal seams.
//!
//! A multi-sector corner cap (PnS bi-4, KP20 polar, exact-boundary
//! Coons wedge, Gregory fan) is built from `N` faces meeting at a
//! pole. Two different questions can be asked of such a cap:
//!
//! - **cap ↔ support**: does each sector meet the fillet strip /
//!   face it borders? This is what the exact-boundary route's
//!   acceptance gates.
//! - **cap ↔ cap**: do adjacent *sectors of the same cap* meet each
//!   other smoothly? This is the seam a user sees as a crease
//!   running out of the corner.
//!
//! Only the first was ever measured in the pipeline, which let a
//! route report `achieved_continuity = G1` while its internal fan was
//! G0-only. This module measures the second, so the claim and the
//! geometry can be compared.
//!
//! # Method
//!
//! For every edge shared by two faces of the given set, both sides
//! are evaluated **through their own pcurves** — not by projecting a
//! point back onto each surface, which silently picks the wrong
//! branch on a folded or periodic chart. Normals and normal
//! curvatures come from finite differences in UV, so the measurement
//! is independent of surface kind: analytic, NURBS, and Gregory caps
//! are all measured the same way.
//!
//! # Degenerate charts
//!
//! Cap sectors meet at a pole where one partial derivative vanishes
//! and a finite-difference normal is pure direction noise. Samples
//! whose chart is that ill-conditioned are **skipped**, not clamped:
//! a single such sample once read 2.9e-1 "G1" on a fan whose two
//! sides agreed positionally to 3e-15. Sampling also stays inside
//! `[0.05, 0.95]` of the seam for the same reason. A caller that
//! wants to know whether a bad number is a pole artifact should use
//! [`seam_normal_profile`] and look at where along the seam it sits.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Fraction of the seam skipped at each end. Adjacent cap sectors
/// meet at a degenerate pole; samples taken there describe the chart,
/// not the surface.
const SEAM_END_MARGIN: f64 = 0.05;

/// Worst-case continuity across the internal seams of a cap.
///
/// All magnitudes carry the model's length units except
/// [`Self::worst_normal_angle_deg`], which is an angle in degrees —
/// `0°` is a smooth join, `90°` a hard crease. Callers compare
/// against their own explicit gate.
#[derive(Debug, Clone, PartialEq)]
pub struct CapSeamContinuity {
    /// Number of cap-cap seams found and measured.
    pub seam_count: usize,
    /// Number of individual samples that produced a usable normal.
    pub sample_count: usize,
    /// Worst positional deviation of either side from the shared
    /// edge curve.
    pub worst_g0_gap: f64,
    /// Worst `|n_a × n_b|` across all samples.
    pub worst_g1_normal: f64,
    /// Worst angle between the two sides' normals, in degrees.
    pub worst_normal_angle_deg: f64,
    /// Worst normal-curvature jump along the seam direction.
    pub worst_g2_residual: f64,
    /// The seam edge carrying [`Self::worst_g1_normal`].
    pub worst_g1_edge: Option<EdgeHandle>,
    /// Fraction along that edge where the worst normal jump occurred.
    /// Useful for telling an interior crease from a pole artifact.
    pub worst_g1_fraction: Option<f64>,
}

impl CapSeamContinuity {
    fn empty() -> Self {
        Self {
            seam_count: 0,
            sample_count: 0,
            worst_g0_gap: 0.0,
            worst_g1_normal: 0.0,
            worst_normal_angle_deg: 0.0,
            worst_g2_residual: 0.0,
            worst_g1_edge: None,
            worst_g1_fraction: None,
        }
    }

    /// Whether any seam was actually measurable.
    ///
    /// A cap that commits a single face, or whose sectors were merged,
    /// has no internal seam — that is a legitimate outcome and must
    /// not be read as "perfectly continuous". Callers asserting a
    /// continuity gate should check this first.
    pub fn measured(&self) -> bool {
        self.seam_count > 0 && self.sample_count > 0
    }
}

/// Measure continuity across every edge shared by two faces in
/// `faces`.
///
/// `samples_per_seam` controls sampling density along each seam;
/// 100 is a reasonable default for gates.
pub fn measure_cap_seam_continuity(
    brep: &BRep,
    faces: &[FaceHandle],
    samples_per_seam: usize,
) -> CapSeamContinuity {
    let set: std::collections::BTreeSet<FaceHandle> = faces.iter().copied().collect();
    let mut out = CapSeamContinuity::empty();
    if samples_per_seam < 2 {
        return out;
    }
    for (edge_handle, edge) in brep.edges.iter() {
        let Some((side_a, side_b, curve)) = seam_sides(brep, edge_handle, edge, &set) else {
            continue;
        };
        out.seam_count += 1;
        for k in 0..samples_per_seam {
            let f = SEAM_END_MARGIN
                + (1.0 - 2.0 * SEAM_END_MARGIN) * k as f64 / (samples_per_seam - 1) as f64;
            let t = curve.t_start + f * (curve.t_end - curve.t_start);
            let p_ref = (&curve.geometry as &dyn CurveEval).point_at(t);
            let (uv_a, p_a) = side_a.eval_fraction(f);
            let (uv_b, p_b) = side_b.eval_fraction(f);
            out.worst_g0_gap = out
                .worst_g0_gap
                .max(p_a.distance_to(&p_ref).max(p_b.distance_to(&p_ref)));

            let (Some(n_a), Some(n_b)) = (side_a.normal(uv_a), side_b.normal(uv_b)) else {
                continue;
            };
            out.sample_count += 1;
            let g1 = n_a.cross(&n_b).length();
            if g1 > out.worst_g1_normal {
                out.worst_g1_normal = g1;
                out.worst_g1_edge = Some(edge_handle);
                out.worst_g1_fraction = Some(f);
            }
            let cos = n_a.dot(&n_b).abs().clamp(0.0, 1.0);
            out.worst_normal_angle_deg = out.worst_normal_angle_deg.max(cos.acos().to_degrees());

            let dt = 1e-4 * (curve.t_end - curve.t_start).abs().max(1e-12);
            let p_fwd = (&curve.geometry as &dyn CurveEval).point_at(t + dt);
            let Some(seam_dir) = (p_fwd - p_ref).normalized() else {
                continue;
            };
            if let (Some(k_a), Some(k_b)) = (
                side_a.normal_curvature(uv_a, seam_dir),
                side_b.normal_curvature(uv_b, seam_dir),
            ) {
                out.worst_g2_residual = out.worst_g2_residual.max((k_a - k_b).abs());
            }
        }
    }
    out
}

/// Angle between the two sides' normals, in degrees, at explicit
/// fractions along one seam.
///
/// This is the tool for attributing a bad worst-case number. A
/// misalignment that climbs steadily toward one end is a
/// degenerate-pole artifact; one that is bad mid-seam is a real
/// crease in the surface. Returns `None` if the two faces share no
/// measurable seam.
pub fn seam_normal_profile(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    fractions: &[f64],
) -> Option<Vec<(f64, f64)>> {
    let set: std::collections::BTreeSet<FaceHandle> = [face_a, face_b].into_iter().collect();
    for (edge_handle, edge) in brep.edges.iter() {
        let Some((side_a, side_b, _curve)) = seam_sides(brep, edge_handle, edge, &set) else {
            continue;
        };
        let mut out = Vec::with_capacity(fractions.len());
        for &f in fractions {
            let (uv_a, _) = side_a.eval_fraction(f);
            let (uv_b, _) = side_b.eval_fraction(f);
            let (Some(n_a), Some(n_b)) = (side_a.normal(uv_a), side_b.normal(uv_b)) else {
                continue;
            };
            let cos = n_a.dot(&n_b).abs().clamp(0.0, 1.0);
            out.push((f, cos.acos().to_degrees()));
        }
        return Some(out);
    }
    None
}

/// Resolve both sides of an edge, if it is a seam between two
/// distinct faces of `set` and both sides carry the pcurve and
/// surface the measurement needs.
#[allow(clippy::type_complexity)]
fn seam_sides<'a>(
    brep: &'a BRep,
    _edge_handle: EdgeHandle,
    edge: &'a popkern_core::topology::Edge,
    set: &std::collections::BTreeSet<FaceHandle>,
) -> Option<(SeamSide<'a>, SeamSide<'a>, &'a CurveOnEdge)> {
    let he = edge.half_edge;
    let he_data = brep.half_edges.get(he)?;
    let twin = he_data.twin;
    // Index defensively: a mid-pipeline BRep can carry a stale `loop_`
    // on a twin half-edge. A measurement must report "nothing to see
    // here", never abort on a SlotMap panic.
    let twin_data = brep.half_edges.get(twin)?;
    let f_a = brep.loops.get(he_data.loop_)?.face;
    let f_b = brep.loops.get(twin_data.loop_)?.face;
    if f_a == f_b || !set.contains(&f_a) || !set.contains(&f_b) {
        return None;
    }
    let curve = edge.curve.as_ref()?;
    let side_a = seam_side(brep, he, curve)?;
    let side_b = seam_side(brep, twin, curve)?;
    Some((side_a, side_b, curve))
}

struct SeamSide<'a> {
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
    fn eval_fraction(&self, f: f64) -> ((f64, f64), Point3) {
        let g = if self.reversed { 1.0 - f } else { f };
        let uv = self
            .pcurve
            .point_at(self.pcurve.t_start + g * (self.pcurve.t_end - self.pcurve.t_start));
        ((uv.x, uv.y), self.surface.point_at(uv.x, uv.y))
    }

    /// Finite-difference unit normal, or `None` where the chart is
    /// too degenerate to trust (see module docs).
    fn normal(&self, uv: (f64, f64)) -> Option<Vec3> {
        let h = 1e-6;
        let du = self.surface.point_at(uv.0 + h, uv.1) - self.surface.point_at(uv.0 - h, uv.1);
        let dv = self.surface.point_at(uv.0, uv.1 + h) - self.surface.point_at(uv.0, uv.1 - h);
        let (lu, lv) = (du.length(), dv.length());
        if lu.min(lv) < 1e-3 * lu.max(lv) {
            return None;
        }
        du.cross(&dv).normalized()
    }

    /// Normal curvature along `dir` by second differences of the
    /// surface's deviation from the tangent plane.
    fn normal_curvature(&self, uv: (f64, f64), dir: Vec3) -> Option<f64> {
        let normal = self.normal(uv)?;
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
        let step = 1e-4;
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
