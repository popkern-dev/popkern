//! Measured continuity evidence for a committed blend (PK-188).
//!
//! Continuity gates were fixture-scoped (apex fixtures at n ∈ {5,6,8})
//! while `achieved_continuity` elsewhere was *inferred from surface
//! kind*, not measured. This module closes that gap: given the faces a
//! blend operation created, it measures the G0/G1/G2 residual across
//! every live seam those faces participate in — blend↔support and
//! blend↔blend — and classifies an achieved [`ContinuityLevel`]
//! against explicit gates. The record distinguishes requested from
//! achieved so a silent downgrade is visible on production-shaped
//! geometry, not only on the fixtures a construction was tuned for.
//!
//! Evaluation is pcurve-true two-sided sampling (shared with
//! [`crate::cap_seam_continuity`]): orientation resolved from
//! endpoints, degenerate charts skipped, finite-difference normals and
//! normal curvatures — independent of surface kind.
//!
//! # Honesty rules
//!
//! - A seam that cannot be measured (missing pcurve/surface) is
//!   counted in [`BlendContinuityEvidence::unmeasured_seams`] and caps
//!   `achieved` at `G0` — an unmeasured seam is never a passing one.
//! - A seam where every sample's chart was too degenerate to produce a
//!   normal contributes positional evidence only and likewise caps
//!   `achieved` at `G0`.
//! - No measurable seams at all reads `achieved = G0` with an empty
//!   seam list; callers must treat that as "no evidence", not
//!   "perfect".
//!
//! # Intended-sharp seams
//!
//! Not every live seam of a blend is supposed to be smooth: an
//! open-ended fillet strip terminates against a transverse face in a
//! crease *by design* (the end cross-section), and gating G1 there
//! would make every open-ended fillet read `G0`. Seam intent is
//! classified by the measured normal angle against
//! [`ContinuityGates::smooth_intent_max_angle_deg`]: within the bound
//! the seam is gated for G1/G2; beyond it the seam is recorded as
//! intended-sharp and gated for G0 only. This is a measurement-side
//! approximation of construction intent (the honest source would be
//! rail provenance threaded through commit — a recorded follow-up):
//! a seam that degrades *past* the bound stops being gated, so
//! corpus baselines must pin the gated/sharp seam census alongside
//! the residuals — a rail collapsing into "sharp" then trips the
//! census, not silence.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use serde::{Deserialize, Serialize};

use crate::cap_seam_measure::{edge_seam_faces, measure_edge_seam, resolve_edge_seam};
use crate::types::ContinuityLevel;

/// Explicit residual gates for classifying measured continuity.
///
/// Tolerances are always explicit parameters in popkern; these are the
/// blend-seam classification bars, carried inside the evidence record
/// so a reader knows what the classification meant.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ContinuityGates {
    /// Positional gap gate (model units): both sides of a seam must
    /// lie on the shared edge curve within this distance.
    pub g0_gap: f64,
    /// Tangent gate: worst `|n_a × n_b|` across seam samples.
    pub g1_normal: f64,
    /// Curvature gate: worst normal-curvature jump along the seam.
    pub g2_curvature: f64,
    /// Intent bound (degrees): a seam whose worst normal angle stays
    /// within this is treated as intended-smooth and gated for G1/G2;
    /// beyond it the seam is recorded as intended-sharp (a designed
    /// crease, e.g. an open fillet end) and gated for G0 only. See
    /// the module docs for why this approximation is safe only when
    /// the gated/sharp census is pinned alongside the residuals.
    pub smooth_intent_max_angle_deg: f64,
}

impl ContinuityGates {
    /// The commercial blend-seam bars: `g0 = 1e-5`, `g1 = 5e-2`,
    /// `g2 = 5e-2` — the same bounds the setback-fan commit
    /// postcondition and the kp-polar commercial gates enforce on
    /// apex fixtures, applied corpus-wide — with a 30° smooth-intent
    /// bound: measured blend rails sit at ≤ 4° while designed creases
    /// start at ~45° (a fillet-strip end against a pyramid lateral
    /// measures 44.7° — an earlier 45° bound misclassified exactly
    /// that crease as a smooth-intent seam and read the part as a
    /// 0.70 G1 failure).
    pub fn commercial() -> Self {
        ContinuityGates {
            g0_gap: 1.0e-5,
            g1_normal: 5.0e-2,
            g2_curvature: 5.0e-2,
            smooth_intent_max_angle_deg: 30.0,
        }
    }
}

/// Which face family sits on the far side of a blend seam.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BlendSeamKind {
    /// Seam between a blend face and a surviving support face — the
    /// boundary the blend was asked to interpolate.
    BlendToSupport,
    /// Seam between two blend faces (strip↔strip, strip↔cap,
    /// cap↔cap) — the operation's internal wiring.
    BlendToBlend,
}

/// Worst residuals measured across one live seam.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SeamContinuityRecord {
    /// The seam edge (provenance for drill-down; valid on the BRep
    /// the measurement ran against).
    pub edge: EdgeHandle,
    /// Far-side family.
    pub kind: BlendSeamKind,
    /// Samples that produced a usable normal on both sides. Zero
    /// means G1/G2 were not measured on this seam.
    pub samples: usize,
    /// Whether this seam is gated for G1/G2 (intended-smooth) or
    /// recorded as an intended-sharp designed crease (G0-gated only).
    pub gated: bool,
    /// The edge's recorded entity tolerance, when set. A construction
    /// that accepts a bounded boundary inexactness records it on the
    /// edge (the same contract the serving curve-surface validation
    /// honors); the G0 classification arm widens to it, and the
    /// record keeps it visible so a widened claim is never silent.
    pub recorded_tolerance: Option<f64>,
    /// Worst positional deviation of either side from the edge curve.
    pub worst_g0_gap: f64,
    /// Worst `|n_a × n_b|` across samples.
    pub worst_g1_normal: f64,
    /// Worst angle between the sides' normals, in degrees (the
    /// intent-classification signal).
    pub worst_normal_angle_deg: f64,
    /// Worst normal-curvature jump across the along-seam, transverse,
    /// and diagonal in-tangent-plane directions (units 1/length).
    ///
    /// The transverse/diagonal terms carry the honest curvature-
    /// continuity content: along-seam agreement is implied by G0+G1
    /// (Meusnier), which is why a rolling-ball rail reads ~0 along
    /// the seam while genuinely jumping `1/r` transversely.
    pub worst_g2_residual: f64,
}

/// The measured continuity record for one blend operation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlendContinuityEvidence {
    /// What the caller asked for.
    pub requested: ContinuityLevel,
    /// What the seams measure, classified against [`Self::gates`] —
    /// independent of `requested`, so both silent downgrades
    /// (`achieved < requested`) and free upgrades are visible.
    pub achieved: ContinuityLevel,
    /// The gates the classification used.
    pub gates: ContinuityGates,
    /// Number of blend faces the measurement covered.
    pub blend_faces: usize,
    /// Per-seam worst residuals.
    pub seams: Vec<SeamContinuityRecord>,
    /// Live seams adjacent to a blend face that could NOT be measured
    /// (missing pcurve or surface on a side). Non-zero caps
    /// [`Self::achieved`] at `G0`.
    pub unmeasured_seams: usize,
    /// Measured seams classified as intended-sharp designed creases
    /// (G0-gated only). Baselines should pin this census — see the
    /// module docs.
    pub sharp_seams: usize,
    /// Seams whose G0 classification used a recorded entity tolerance
    /// wider than [`ContinuityGates::g0_gap`] — visible, never silent.
    pub tolerance_widened_seams: usize,
    /// Worst positional gap across all measured seams (sharp
    /// included — watertightness has no intent exemption).
    pub worst_g0_gap: f64,
    /// Worst tangent residual across GATED seams.
    pub worst_g1_normal: f64,
    /// Worst curvature jump across GATED seams.
    pub worst_g2_residual: f64,
}

impl BlendContinuityEvidence {
    /// Whether any seam was actually measured. `false` means the
    /// record carries no continuity claim at all.
    pub fn measured(&self) -> bool {
        !self.seams.is_empty()
    }

    /// Whether the blend delivered less than was requested.
    pub fn downgraded(&self) -> bool {
        self.achieved < self.requested
    }
}

/// Measure the continuity evidence for the faces a blend created.
///
/// Walks every edge with `blend_faces` on at least one side (and two
/// distinct faces), samples both sides through their own pcurves, and
/// classifies the worst residuals against `gates`.
/// `samples_per_seam` is the sampling density (33 is a reasonable
/// gate default; the apex fixtures use 100).
pub fn measure_blend_continuity_evidence(
    brep: &BRep,
    blend_faces: &[FaceHandle],
    requested: ContinuityLevel,
    gates: &ContinuityGates,
    samples_per_seam: usize,
) -> BlendContinuityEvidence {
    let set: std::collections::BTreeSet<FaceHandle> = blend_faces.iter().copied().collect();
    let mut seams = Vec::new();
    let mut unmeasured = 0usize;
    for (edge_handle, edge) in brep.edges.iter() {
        let Some(adjacency) = edge_seam_faces(brep, edge) else {
            continue;
        };
        let a_in = set.contains(&adjacency.face_a);
        let b_in = set.contains(&adjacency.face_b);
        if !a_in && !b_in {
            continue;
        }
        let kind = if a_in && b_in {
            BlendSeamKind::BlendToBlend
        } else {
            BlendSeamKind::BlendToSupport
        };
        let Some(seam) = resolve_edge_seam(brep, edge, &adjacency) else {
            unmeasured += 1;
            continue;
        };
        let worst = measure_edge_seam(&seam, samples_per_seam);
        seams.push(SeamContinuityRecord {
            edge: edge_handle,
            kind,
            samples: worst.samples,
            gated: worst.angle_deg <= gates.smooth_intent_max_angle_deg,
            recorded_tolerance: edge.tolerance,
            worst_g0_gap: worst.g0,
            worst_g1_normal: worst.g1,
            worst_normal_angle_deg: worst.angle_deg,
            worst_g2_residual: worst.g2.max(worst.g2_transverse),
        });
    }

    let sharp_seams = seams.iter().filter(|s| !s.gated).count();
    let tolerance_widened_seams = seams
        .iter()
        .filter(|s| {
            s.recorded_tolerance
                .map(|t| t > gates.g0_gap)
                .unwrap_or(false)
                && s.worst_g0_gap > gates.g0_gap
        })
        .count();
    let worst_g0_gap = seams.iter().map(|s| s.worst_g0_gap).fold(0.0, f64::max);
    let gated = || seams.iter().filter(|s| s.gated);
    let worst_g1_normal = gated().map(|s| s.worst_g1_normal).fold(0.0, f64::max);
    let worst_g2_residual = gated().map(|s| s.worst_g2_residual).fold(0.0, f64::max);
    let achieved = classify(&seams, unmeasured, gates);

    BlendContinuityEvidence {
        requested,
        achieved,
        gates: *gates,
        blend_faces: set.len(),
        seams,
        unmeasured_seams: unmeasured,
        sharp_seams,
        tolerance_widened_seams,
        worst_g0_gap,
        worst_g1_normal,
        worst_g2_residual,
    }
}

/// Classify the highest [`ContinuityLevel`] the measured seams pass:
/// G0 over every seam, G1/G2 over the gated (intended-smooth) ones.
/// See the module honesty rules for why unmeasured or normal-less
/// seams floor the result at `G0`, and why "no gated seams" is no
/// claim rather than a pass.
fn classify(
    seams: &[SeamContinuityRecord],
    unmeasured: usize,
    gates: &ContinuityGates,
) -> ContinuityLevel {
    if seams.is_empty() || unmeasured > 0 || seams.iter().any(|s| s.samples == 0) {
        return ContinuityLevel::G0;
    }
    // The G0 arm honors a recorded entity tolerance (the same
    // contract as the serving curve-surface validation): an accepted,
    // recorded boundary inexactness is not a positional failure. The
    // widening is counted in `tolerance_widened_seams`.
    if seams
        .iter()
        .any(|s| s.worst_g0_gap > gates.g0_gap.max(s.recorded_tolerance.unwrap_or(0.0)))
    {
        return ContinuityLevel::G0;
    }
    let gated: Vec<_> = seams.iter().filter(|s| s.gated).collect();
    if gated.is_empty() {
        return ContinuityLevel::G0;
    }
    let g1_ok = gated.iter().all(|s| s.worst_g1_normal <= gates.g1_normal);
    let g2_ok = gated
        .iter()
        .all(|s| s.worst_g2_residual <= gates.g2_curvature);
    if g1_ok && g2_ok {
        ContinuityLevel::G2
    } else if g1_ok {
        ContinuityLevel::G1
    } else {
        ContinuityLevel::G0
    }
}

#[cfg(test)]
#[path = "blend_continuity_evidence_tests.rs"]
mod tests;
