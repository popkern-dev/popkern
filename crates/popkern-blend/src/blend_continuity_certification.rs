//! Certified bounds for blend continuity evidence (PK-499 stage 4).
//!
//! Every gated (intended-smooth) seam is certified with
//! [`popkern_core::geometry::certify`]: rigorous upper bounds on the
//! position gap, normal-angle sine and second-fundamental-form jump
//! over the whole edge, replacing the sampled maxima in the
//! classification wherever the certificate holds. The sampled maxima
//! stay on the record as diagnostics.
//!
//! # Fallback boundary
//!
//! A certificate stopped by budget with nothing uncertified still
//! holds rigorous bounds over the whole seam and is used. A seam whose
//! certificate does not cover it (a pole or blended Gregory corner the
//! enclosures cannot bound, an unsupported surface kind, a failed
//! correspondence) is classified from its sampled
//! maxima, counted in `uncertified_seams`, and keeps its certificate
//! (with the reason) on the record — never a silent pass.
//! `certified_achieved` is set only when every gated seam certified.
//!
//! # Targets and budget
//!
//! Evidence decides gates, it does not pin residuals: refinement runs
//! in decision mode against the commercial gates (a quantity certified
//! past its gate, or an interval certified within it, stops refining;
//! a small-radius rail's certified 1/r jump no longer drives ~200
//! intervals of curvature tightening), and
//! targets are a tenth of the G1/G2 gates and half the G0 gate (a
//! gap refined to 1e-6 on a radius-10 ring cost twice the intervals
//! of the decision it serves), and the sides count as
//! one seam within the G0 gate (the branch and bound settles an exact
//! seam in a handful of intervals), under a small per-seam budget so
//! part construction stays cheap.

use popkern_core::geometry::certify::{CertifyTargets, SeamCertificate, SeamCertifyOptions};
use popkern_core::topology::edge_seam_certificate::{certify_edge_seam, seam_certify_options};
use popkern_core::topology::{BRep, EdgeHandle};

use serde::{Deserialize, Serialize};

use crate::blend_continuity_evidence::{BlendContinuityEvidence, ContinuityGates};

/// When a blend's gated seams are certified.
///
/// Construction never depends on the certificate: it only upgrades a
/// sampled classification to a rigorous one (or overturns it). A host
/// that must show the part first (popkern live) defers
/// certification and runs it afterwards, seam by seam, against the
/// same BRep; batch callers (the CLI, tests, the part API) certify
/// inline.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ContinuityCertification {
    /// Certify every gated seam while measuring (the default).
    #[default]
    Inline,
    /// Record sampled evidence now and mark gated seams deferred; see
    /// [`BlendContinuityEvidence::certify_next_deferred`].
    Deferred,
}

/// One deferred seam certified by
/// [`BlendContinuityEvidence::certify_next_deferred`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeferredSeamCertified {
    /// The seam edge (valid on the BRep the evidence was measured on).
    pub edge: EdgeHandle,
    /// Deferred seams still left on the record.
    pub remaining: usize,
}

impl BlendContinuityEvidence {
    /// Certify the first deferred seam against `brep` — the BRep this
    /// evidence was measured on (its edge handles index that BRep) —
    /// and re-classify the record. `None` when nothing is deferred.
    ///
    /// Running this until it returns `None` yields exactly the record
    /// [`crate::measure_blend_continuity_evidence`] builds inline:
    /// same options, same certificates, same classification.
    pub fn certify_next_deferred(&mut self, brep: &BRep) -> Option<DeferredSeamCertified> {
        let options = evidence_certify_options(brep, &self.gates);
        let seam = self.seams.iter_mut().find(|s| s.certification_deferred)?;
        seam.certificate = Some(certify_seam_for_evidence(brep, seam.edge, &options));
        seam.certification_deferred = false;
        let edge = seam.edge;
        self.summarize();
        Some(DeferredSeamCertified {
            edge,
            remaining: self.deferred_seams,
        })
    }
}

/// Per-seam interval budget for evidence certification. The corpus
/// median certifies in 4 intervals (PK-499 E1); a projected rail (the
/// PK-389 rim) needs ~1/64 intervals before its per-branch Krawczyk
/// contracts, ~150 evaluations. Seams needing more fall back to
/// sampled classification, recorded.
pub const EVIDENCE_CERTIFY_BUDGET: usize = 256;

/// Certification options for evidence on `brep` against `gates`.
pub fn evidence_certify_options(brep: &BRep, gates: &ContinuityGates) -> SeamCertifyOptions {
    let base = seam_certify_options(brep);
    SeamCertifyOptions {
        targets: CertifyTargets {
            gap: 0.5 * gates.g0_gap,
            normal_sine: 0.1 * gates.g1_normal,
            curvature: 0.1 * gates.g2_curvature,
            relative: 1e-2,
        },
        // Sides within the G0 gate are one seam for the evidence.
        correspondence_tolerance: gates.g0_gap,
        budget: EVIDENCE_CERTIFY_BUDGET,
        // Evidence uses only covering bounds; a structural failure
        // falls back to sampled at once.
        stop_at_structural_failure: true,
        decision_gates: Some([gates.g0_gap, gates.g1_normal, gates.g2_curvature]),
        ..base
    }
}

/// Certify one seam for the evidence record.
pub(crate) fn certify_seam_for_evidence(
    brep: &BRep,
    edge: EdgeHandle,
    options: &SeamCertifyOptions,
) -> SeamCertificate {
    certify_edge_seam(brep, edge, options)
}
