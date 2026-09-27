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

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::cap_seam_measure::{edge_seam_faces, measure_edge_seam, resolve_edge_seam};

/// Fraction of the seam skipped at each end. Adjacent cap sectors
/// meet at a degenerate pole; samples taken there describe the chart,
/// not the surface.
pub(crate) const SEAM_END_MARGIN: f64 = 0.05;

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
    let jets = crate::cap_seam_side::SurfaceJetCache::default();
    let set: std::collections::BTreeSet<FaceHandle> = faces.iter().copied().collect();
    let mut out = CapSeamContinuity::empty();
    if samples_per_seam < 2 {
        return out;
    }
    for (edge_handle, edge) in brep.edges.iter() {
        let Some(adjacency) = edge_seam_faces(brep, edge) else {
            continue;
        };
        if !set.contains(&adjacency.face_a) || !set.contains(&adjacency.face_b) {
            continue;
        }
        let Some(seam) = resolve_edge_seam(brep, edge, &adjacency, &jets) else {
            continue;
        };
        out.seam_count += 1;
        let worst = measure_edge_seam(&seam, samples_per_seam, SEAM_END_MARGIN);
        out.sample_count += worst.samples;
        out.worst_g0_gap = out.worst_g0_gap.max(worst.g0);
        if worst.g1 > out.worst_g1_normal {
            out.worst_g1_normal = worst.g1;
            out.worst_g1_edge = Some(edge_handle);
            out.worst_g1_fraction = worst.g1_fraction;
        }
        out.worst_normal_angle_deg = out.worst_normal_angle_deg.max(worst.angle_deg);
        out.worst_g2_residual = out.worst_g2_residual.max(worst.g2);
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
    let jets = crate::cap_seam_side::SurfaceJetCache::default();
    let set: std::collections::BTreeSet<FaceHandle> = [face_a, face_b].into_iter().collect();
    for (_edge_handle, edge) in brep.edges.iter() {
        let Some(adjacency) = edge_seam_faces(brep, edge) else {
            continue;
        };
        if !set.contains(&adjacency.face_a) || !set.contains(&adjacency.face_b) {
            continue;
        }
        let Some(seam) = resolve_edge_seam(brep, edge, &adjacency, &jets) else {
            continue;
        };
        let mut out = Vec::with_capacity(fractions.len());
        for &f in fractions {
            let (uv_a, _) = seam.side_a.eval_fraction(f);
            let (uv_b, _) = seam.side_b.eval_fraction(f);
            let (Some(n_a), Some(n_b)) = (seam.side_a.normal(uv_a), seam.side_b.normal(uv_b))
            else {
                continue;
            };
            let cos = n_a.dot(&n_b).abs().clamp(0.0, 1.0);
            out.push((f, cos.acos().to_degrees()));
        }
        return Some(out);
    }
    None
}
