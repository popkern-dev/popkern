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
//! would make every open-ended fillet read `G0`.
//!
//! - **Blend-to-support seams** are classified structurally (PK-592): a
//!   seam that runs tangent somewhere, or creases by a constant sub-bound
//!   angle on a curved blend (a degraded rail), is gated; a strip end
//!   (crease varying with the profile), a planar fragment, or a crease
//!   past the bound is designed-sharp. See `blend_continuity_intent`.
//! - **Blend-to-blend seams** keep the angle bound, except inside a corner
//!   cap: gated at any angle between two blend-enclosed faces (PK-609).
//!
//! Both are measurement-side readings of construction intent (the
//! honest source would be rail provenance threaded through commit). A
//! seam misread as sharp stops being gated, so corpus baselines pin the
//! gated/sharp seam census alongside the residuals: a rail collapsing
//! into "sharp" then trips the census, not silence.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use serde::{Deserialize, Serialize};

use popkern_core::geometry::certify::SeamCertificate;

pub use crate::blend_continuity_certification::ContinuityCertification;
use crate::blend_continuity_certification::{certify_gated_seams, evidence_certify_options};
use crate::cap_seam_measure::{edge_seam_faces, measure_edge_seam, resolve_edge_seam};
use crate::types::ContinuityLevel;

#[path = "blend_continuity_intent.rs"]
mod intent;

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
    /// blend-to-blend intent-classification signal).
    pub worst_normal_angle_deg: f64,
    /// Smallest `|n_a × n_b|` across samples: whether the sides run
    /// tangent anywhere along the seam (the blend-to-support intent
    /// signal, PK-592). Infinite when no sample gave both normals.
    #[serde(default = "intent::no_tangency")]
    pub best_g1_normal: f64,
    /// Worst normal-curvature jump over every in-tangent-plane
    /// direction: the spectral norm of the second-fundamental-form jump,
    /// assembled from the along-seam, transverse and diagonal readings
    /// (units 1/length; PK-579).
    ///
    /// The transverse/diagonal terms carry the honest curvature-
    /// continuity content: along-seam agreement is implied by G0+G1
    /// (Meusnier), which is why a rolling-ball rail reads ~0 along
    /// the seam while genuinely jumping `1/r` transversely.
    pub worst_g2_residual: f64,
    /// Certified bounds over the whole seam (PK-499), for gated seams;
    /// `None` for intended-sharp seams, which are G0-gated only. When
    /// certified, its upper bounds replace the sampled maxima in the
    /// classification.
    #[serde(default)]
    pub certificate: Option<SeamCertificate>,
    /// Whether this gated seam's certification was deferred by the
    /// caller ([`ContinuityCertification::Deferred`]): it carries no
    /// certificate yet, is classified from its sampled maxima, and is
    /// counted in [`BlendContinuityEvidence::deferred_seams`] until
    /// [`BlendContinuityEvidence::certify_next_deferred`] runs on it.
    #[serde(default)]
    pub certification_deferred: bool,
    /// Samples whose transverse/diagonal curvature could not be read
    /// (a singular chart point). The sampled G2 residual covers every
    /// other sample; a nonzero count says where it is blind (PK-579).
    #[serde(default)]
    pub skipped_curvature_samples: usize,
}

impl SeamContinuityRecord {
    /// The certificate, when its upper bounds hold over the whole seam
    /// ([`SeamCertificate::bounds_cover_seam`]).
    pub fn certified(&self) -> Option<&SeamCertificate> {
        self.certificate.as_ref().filter(|c| c.bounds_cover_seam())
    }

    /// Per quantity (gap, sine, curvature): the certified decision
    /// against its gate, `Some(value)` when decided — the upper bound
    /// when it passes (`upper ≤ gate`), the certified lower bound when
    /// it fails (`lower > gate`) — and `None` when the certificate
    /// does not cover the seam or its bounds straddle the gate.
    fn certified_decisions(&self, gates: &ContinuityGates) -> [Option<f64>; 3] {
        let Some(c) = self.certified() else {
            return [None; 3];
        };
        let g0_gate = gates.g0_gap.max(self.recorded_tolerance.unwrap_or(0.0));
        let decide = |upper: f64, lower: f64, gate: f64| {
            if upper <= gate {
                Some(upper)
            } else if lower > gate {
                Some(lower)
            } else {
                None
            }
        };
        [
            decide(c.g0_upper, c.lower[0], g0_gate),
            decide(c.g1_upper, c.lower[1], gates.g1_normal),
            decide(c.g2_upper, c.lower[2], gates.g2_curvature),
        ]
    }

    /// Whether the certificate decides every gate for this seam.
    pub fn certificate_decides(&self, gates: &ContinuityGates) -> bool {
        self.certified_decisions(gates).iter().all(Option::is_some)
    }

    /// `(g0, g1, g2)` the classification uses: each certified decision
    /// where the certificate makes one, the sampled maximum otherwise
    /// (the recorded fallback). The gap always keeps the sampled
    /// distance to the edge curve too (the trim contract).
    pub fn classified_residuals(&self, gates: &ContinuityGates) -> (f64, f64, f64) {
        let [g0, g1, g2] = self.certified_decisions(gates);
        (
            self.worst_g0_gap.max(g0.unwrap_or(0.0)),
            g1.unwrap_or(self.worst_g1_normal),
            g2.unwrap_or(self.worst_g2_residual),
        )
    }
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
    /// Gated seams whose certificate decides every gate (bounds over
    /// the whole edge, each quantity certified passing or failing).
    #[serde(default)]
    pub certified_seams: usize,
    /// Gated seams classified (at least partly) from sampled maxima
    /// because their certificate does not decide every gate (the
    /// recorded fallback boundary). Deferred seams are not counted
    /// here: they have not been certified yet.
    #[serde(default)]
    pub uncertified_seams: usize,
    /// Gated seams whose certification the caller deferred (see
    /// [`ContinuityCertification::Deferred`]); classified from sampled
    /// maxima until certified.
    #[serde(default)]
    pub deferred_seams: usize,
    /// The level the certified bounds alone support: `Some` only when
    /// every gated seam certified (none deferred, none undecided).
    #[serde(default)]
    pub certified_achieved: Option<ContinuityLevel>,
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

/// Seam end inset of the evidence instrument: 1e-3 of the seam. The
/// sampled readings reach into the seam's end bands, where vertex-blend
/// corner residuals live (PK-579: the former 5% inset hid a certified
/// g2 >= 0.29 in the first 3% of pentagon5_blunt's rim seams). The
/// inset only steps off the seam ENDPOINT itself: a seam endpoint is a
/// patch corner, and a rational corner cell (dual-jet, dual-twist) is
/// 0/0 exactly there. Its ε-regularized exact form read 1e21 at f = 0
/// on the flagship. The regularized neighbourhood is ~1e-7 in parameter.
const EVIDENCE_END_MARGIN: f64 = 1.0e-3;

/// Measure the continuity evidence for the faces a blend created.
///
/// Walks every edge with `blend_faces` on at least one side (and two
/// distinct faces), samples both sides through their own pcurves, and
/// classifies the worst residuals against `gates`.
/// `samples_per_seam` is the sampling density (33 is a reasonable
/// gate default; the apex fixtures use 100). Every gated seam is
/// certified inline ([`ContinuityCertification::Inline`]).
pub fn measure_blend_continuity_evidence(
    brep: &BRep,
    blend_faces: &[FaceHandle],
    requested: ContinuityLevel,
    gates: &ContinuityGates,
    samples_per_seam: usize,
) -> BlendContinuityEvidence {
    measure_blend_continuity_evidence_with(
        brep,
        blend_faces,
        requested,
        gates,
        samples_per_seam,
        ContinuityCertification::Inline,
    )
}

/// [`measure_blend_continuity_evidence`] with an explicit
/// certification policy. Under [`ContinuityCertification::Deferred`]
/// the sampled measurement is complete and every gated seam is marked
/// deferred; certify them later against the same `brep` with
/// [`BlendContinuityEvidence::certify_next_deferred`].
pub fn measure_blend_continuity_evidence_with(
    brep: &BRep,
    blend_faces: &[FaceHandle],
    requested: ContinuityLevel,
    gates: &ContinuityGates,
    samples_per_seam: usize,
    certification: ContinuityCertification,
) -> BlendContinuityEvidence {
    let jets = crate::cap_seam_side::SurfaceJetCache::default();
    let set: std::collections::BTreeSet<FaceHandle> = blend_faces.iter().copied().collect();
    let seam_intent = intent::SeamIntent::new(brep, &set, gates);
    let certify = evidence_certify_options(brep, gates);
    let inline = certification == ContinuityCertification::Inline;
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
        let Some(seam) = resolve_edge_seam(brep, edge, &adjacency, &jets) else {
            unmeasured += 1;
            continue;
        };
        let worst = measure_edge_seam(&seam, samples_per_seam, EVIDENCE_END_MARGIN);
        let blend_face = [adjacency.face_b, adjacency.face_a][usize::from(a_in)];
        let gated = seam_intent.gated(kind, &worst, &adjacency, blend_face);
        seams.push(SeamContinuityRecord {
            edge: edge_handle,
            kind,
            samples: worst.samples,
            gated,
            recorded_tolerance: edge.tolerance,
            worst_g0_gap: worst.g0,
            worst_g1_normal: worst.g1,
            worst_normal_angle_deg: worst.angle_deg,
            best_g1_normal: worst.best_g1,
            worst_g2_residual: worst.g2.max(worst.g2_transverse),
            certificate: None,
            certification_deferred: gated && !inline,
            skipped_curvature_samples: worst.skipped_curvature,
        });
    }
    certify_gated_seams(brep, &mut seams, &certify, inline);
    let mut evidence = BlendContinuityEvidence {
        requested,
        achieved: ContinuityLevel::G0,
        gates: *gates,
        blend_faces: set.len(),
        seams,
        unmeasured_seams: unmeasured,
        sharp_seams: 0,
        tolerance_widened_seams: 0,
        worst_g0_gap: 0.0,
        worst_g1_normal: 0.0,
        worst_g2_residual: 0.0,
        certified_seams: 0,
        uncertified_seams: 0,
        deferred_seams: 0,
        certified_achieved: None,
    };
    evidence.summarize();
    evidence
}

impl BlendContinuityEvidence {
    /// Recompute every summary field from the per-seam records — the
    /// one classification both the inline measurement and deferred
    /// certification go through.
    pub(crate) fn summarize(&mut self) {
        let (seams, gates) = (&self.seams, &self.gates);
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
        let achieved = classify(seams, self.unmeasured_seams, gates);
        let certified_seams = gated().filter(|s| s.certificate_decides(gates)).count();
        let deferred_seams = gated().filter(|s| s.certification_deferred).count();
        let uncertified_seams = gated().count() - certified_seams - deferred_seams;
        let certified_achieved =
            (uncertified_seams == 0 && deferred_seams == 0).then_some(achieved);
        self.sharp_seams = sharp_seams;
        self.tolerance_widened_seams = tolerance_widened_seams;
        self.worst_g0_gap = worst_g0_gap;
        self.worst_g1_normal = worst_g1_normal;
        self.worst_g2_residual = worst_g2_residual;
        self.achieved = achieved;
        self.certified_seams = certified_seams;
        self.uncertified_seams = uncertified_seams;
        self.deferred_seams = deferred_seams;
        self.certified_achieved = certified_achieved;
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
    if seams.iter().any(|s| {
        s.classified_residuals(gates).0 > gates.g0_gap.max(s.recorded_tolerance.unwrap_or(0.0))
    }) {
        return ContinuityLevel::G0;
    }
    let gated: Vec<_> = seams.iter().filter(|s| s.gated).collect();
    if gated.is_empty() {
        return ContinuityLevel::G0;
    }
    let g1_ok = gated
        .iter()
        .all(|s| s.classified_residuals(gates).1 <= gates.g1_normal);
    let g2_ok = gated
        .iter()
        .all(|s| s.classified_residuals(gates).2 <= gates.g2_curvature);
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

#[cfg(test)]
#[path = "blend_continuity_negative_controls_tests.rs"]
mod negative_controls;
