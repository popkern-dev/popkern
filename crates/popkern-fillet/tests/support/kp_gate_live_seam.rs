//! Topology-honest cap-ring metrics: measure continuity across the
//! LIVE seam edges between cap faces instead of assuming a patch
//! ordering.
//!
//! `measure_ring` (kp_gate_helpers) encodes the legacy polar layout
//! — patch `k` side `u = 1` meets patch `k + 1` side `v = 0`
//! reversed, faces emitted in ring order. The setback vertex blend
//! commits a 24-patch fan on hexagon6_variable whose true pairwise
//! seams are machine-zero, yet the legacy convention reads 5.2e-2
//! (frame artifact, PK-087 stage 13). Here every edge shared by two
//! cap faces is measured through the faces' own pcurves: G0 as the
//! cross-face surface gap, G1 as normal misalignment, G2 as the
//! normal-curvature jump along the seam direction (both by finite
//! differences in UV, so any surface kind works).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Worst-case continuity over all cap-cap seams, plus cross-seam
/// normal-alignment statistics (|cos| of the two sides' normals per
/// sample) for scorecard-style consumers.
pub struct LiveSeamMetrics {
    pub seam_count: usize,
    pub worst_g0_gap: f64,
    pub worst_g1_normal: f64,
    pub worst_g2_residual: f64,
    pub worst_g0_edge: Option<EdgeHandle>,
    pub cos_samples: usize,
    pub cos_min: f64,
    pub cos_mean: f64,
    pub cos_max: f64,
}

/// Measure every edge whose two faces are both in `faces`.
/// Cross-seam normal alignment sampled at explicit fractions along
/// the shared curve, for localising a bad seam.
///
/// `measure_live_cap_seams` reports only the worst value over its
/// sampling window (`f in [0.05, 0.95]`). When that worst value is
/// bad, the next question is always *where*: a misalignment that
/// climbs steadily toward one end is a degenerate-pole artifact,
/// while one that is bad mid-seam is a real crease. This returns the
/// `(f, |cos|)` profile so the two can be told apart.
#[allow(dead_code)]
pub fn seam_cos_profile(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    fractions: &[f64],
) -> Option<Vec<(f64, f64)>> {
    for (_eh, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin_he) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin_he.loop_)) else {
            continue;
        };
        let (fa, fb) = (la.face, lb.face);
        if !((fa == face_a && fb == face_b) || (fa == face_b && fb == face_a)) {
            continue;
        }
        let curve = edge.curve.as_ref()?;
        let side_a = seam_side(brep, edge.half_edge, curve)?;
        let side_b = seam_side(brep, he.twin, curve)?;
        let mut out = Vec::new();
        for &f in fractions {
            let (uv_a, _) = side_a.eval_fraction(f);
            let (uv_b, _) = side_b.eval_fraction(f);
            let (Some(n_a), Some(n_b)) = (side_a.normal(uv_a), side_b.normal(uv_b)) else {
                continue;
            };
            out.push((f, n_a.dot(&n_b).abs()));
        }
        return Some(out);
    }
    None
}

pub fn measure_live_cap_seams(
    brep: &BRep,
    faces: &[FaceHandle],
    n_samples: usize,
) -> LiveSeamMetrics {
    let cap: std::collections::BTreeSet<FaceHandle> = faces.iter().copied().collect();
    let mut metrics = LiveSeamMetrics {
        seam_count: 0,
        worst_g0_gap: 0.0,
        worst_g1_normal: 0.0,
        worst_g2_residual: 0.0,
        worst_g0_edge: None,
        cos_samples: 0,
        cos_min: 1.0,
        cos_mean: 1.0,
        cos_max: 1.0,
    };
    let mut cos_sum = 0.0_f64;
    for (eh, edge) in brep.edges.iter() {
        let he = edge.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        // Index defensively: a mid-pipeline BRep can carry a stale
        // `loop_` on a twin half-edge (the PK-074 phantom-seam class).
        // A measurement helper must report that as "no measurable
        // seam here", never abort the whole gate with a SlotMap panic.
        let (Some(loop_a), Some(twin_he)) =
            (brep.loops.get(he_data.loop_), brep.half_edges.get(twin))
        else {
            continue;
        };
        let Some(loop_b) = brep.loops.get(twin_he.loop_) else {
            continue;
        };
        let (f_a, f_b) = (loop_a.face, loop_b.face);
        if f_a == f_b || !cap.contains(&f_a) || !cap.contains(&f_b) {
            continue;
        }
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let (Some(side_a), Some(side_b)) =
            (seam_side(brep, he, curve), seam_side(brep, twin, curve))
        else {
            continue;
        };
        metrics.seam_count += 1;
        for k in 0..n_samples {
            let f = 0.05 + 0.90 * k as f64 / (n_samples - 1) as f64;
            let t = curve.t_start + f * (curve.t_end - curve.t_start);
            let p_ref = (&curve.geometry as &dyn CurveEval).point_at(t);
            let (uv_a, p_a) = side_a.eval_fraction(f);
            let (uv_b, p_b) = side_b.eval_fraction(f);
            // G0 is trim consistency per side against the SHARED
            // curve — the contract validation and tessellation
            // consume. Cross-side same-fraction equality would
            // additionally demand affinely-aligned parameterizations,
            // which the kernel does not require.
            let g0 = p_a.distance_to(&p_ref).max(p_b.distance_to(&p_ref));
            if std::env::var_os("POPKERN_LIVE_SEAM_SAMPLES").is_some() && g0 > 1e-3 {
                eprintln!(
                    "[seam-sample] edge {eh:?} f={f:.3} d_a={:.3e} d_b={:.3e} uv_a=({:.4},{:.4}) uv_b=({:.4},{:.4})",
                    p_a.distance_to(&p_ref),
                    p_b.distance_to(&p_ref),
                    uv_a.0,
                    uv_a.1,
                    uv_b.0,
                    uv_b.1
                );
            }
            if g0 > metrics.worst_g0_gap {
                metrics.worst_g0_gap = g0;
                metrics.worst_g0_edge = Some(eh);
            }
            let (Some(n_a), Some(n_b)) = (side_a.normal(uv_a), side_b.normal(uv_b)) else {
                continue;
            };
            let g1 = n_a.cross(&n_b).length();
            metrics.worst_g1_normal = metrics.worst_g1_normal.max(g1);
            let cos = n_a.dot(&n_b).abs();
            metrics.cos_samples += 1;
            metrics.cos_min = metrics.cos_min.min(cos);
            metrics.cos_max = if metrics.cos_samples == 1 {
                cos
            } else {
                metrics.cos_max.max(cos)
            };
            cos_sum += cos;
            // Seam direction from the reference curve.
            let dt = 1e-4 * (curve.t_end - curve.t_start).abs().max(1e-12);
            let p_fwd = (&curve.geometry as &dyn CurveEval).point_at(t + dt);
            let Some(seam_dir) = (p_fwd - p_ref).normalized() else {
                continue;
            };
            let (Some(k_a), Some(k_b)) = (
                side_a.normal_curvature(uv_a, seam_dir),
                side_b.normal_curvature(uv_b, seam_dir),
            ) else {
                continue;
            };
            let kmax = k_a.abs().max(k_b.abs());
            let g2 = (k_a - k_b).abs() / (1.0 + kmax);
            metrics.worst_g2_residual = metrics.worst_g2_residual.max(g2);
        }
        if std::env::var_os("POPKERN_LIVE_SEAM_TRACE").is_some() {
            let mid =
                (&curve.geometry as &dyn CurveEval).point_at(0.5 * (curve.t_start + curve.t_end));
            let mut seam_g0 = 0.0_f64;
            let mut seam_g1 = 0.0_f64;
            let mut seam_g1_at = 0.0_f64;
            for k in 0..n_samples {
                let f = 0.05 + 0.90 * k as f64 / (n_samples - 1) as f64;
                let t = curve.t_start + f * (curve.t_end - curve.t_start);
                let p_ref = (&curve.geometry as &dyn CurveEval).point_at(t);
                let (uv_a, p_a) = side_a.eval_fraction(f);
                let (uv_b, p_b) = side_b.eval_fraction(f);
                seam_g0 = seam_g0
                    .max(p_a.distance_to(&p_ref))
                    .max(p_b.distance_to(&p_ref));
                if let (Some(n_a), Some(n_b)) = (side_a.normal(uv_a), side_b.normal(uv_b)) {
                    let g1 = n_a.cross(&n_b).length();
                    if g1 > seam_g1 {
                        seam_g1 = g1;
                        seam_g1_at = f;
                    }
                }
            }
            eprintln!(
                "[live-seam] edge {eh:?} faces ({f_a:?},{f_b:?}) g0={seam_g0:.3e} \
                 g1={seam_g1:.3e}@f={seam_g1_at:.2} mid=({:.4},{:.4},{:.4})",
                mid.x, mid.y, mid.z
            );
        }
    }
    if metrics.cos_samples > 0 {
        metrics.cos_mean = cos_sum / metrics.cos_samples as f64;
    }
    metrics
}

struct SeamSide<'a> {
    surface: &'a SurfaceGeometry,
    pcurve: &'a popkern_core::geometry::pcurve::Pcurve,
    /// Pcurve runs opposite to the edge curve (twin half-edges trace
    /// the loop direction, not the curve direction) — resolved once
    /// per side by matching pcurve-endpoint surface images to the
    /// curve endpoints.
    reversed: bool,
}

fn seam_side<'a>(
    brep: &'a BRep,
    he: HalfEdgeHandle,
    curve: &popkern_core::geometry::curve::CurveOnEdge,
) -> Option<SeamSide<'a>> {
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
    /// too degenerate to trust: near a collapsed patch corner one
    /// partial vanishes and the FD normal is direction noise — a
    /// single such end sample read 2.9e-1 "G1" on pentagon5_blunt's
    /// fan while both sides' positions agreed to 3e-15. Skipping the
    /// sample (like an eval failure) measures the seam at its
    /// well-conditioned samples instead of reporting chart noise.
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
        // Solve dir ≈ a·su + b·sv in the tangent plane (2x2 normal
        // equations).
        let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let (r1, r2) = (dir.dot(&su), dir.dot(&sv));
        let det = e * g - f_ * f_;
        if det.abs() < 1e-30 {
            return None;
        }
        let a = (r1 * g - r2 * f_) / det;
        let b = (r2 * e - r1 * f_) / det;
        // Second difference of the surface along (a, b) in UV.
        let step = 1e-4;
        let p0 = self.surface.point_at(uv.0, uv.1);
        let p1 = self.surface.point_at(uv.0 + a * step, uv.1 + b * step);
        let p2 = self.surface.point_at(uv.0 - a * step, uv.1 - b * step);
        let second = (p1 - p0) + (p2 - p0);
        let first_len = (p1 - p2).length() * 0.5;
        if first_len < 1e-30 {
            return None;
        }
        // κ_n = (S'' · n) / |S'|².
        Some(second.dot(&normal) / (first_len * first_len))
    }
}
