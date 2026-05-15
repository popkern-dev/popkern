//! Transactional pcurve reconstruction from authoritative 3D edge curves.

use crate::geometry::trim_domain::{
    TrimDomainCandidateKind, TrimDomainDeclineReason, classify_trim_domain_candidate,
};
use crate::topology::brep_curve_trim_reconstruction_domain::{
    DomainExpansionCommit, expand_repaired_face_domain,
};
use crate::topology::brep_curve_trim_reconstruction_helpers::{
    half_edge_endpoints, project_edge_curve_to_surface_pcurve, sampled_edge_surface_disagreement,
};
use crate::topology::{BRep, HalfEdgeHandle};

/// Options for curve-backed trim reconstruction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepCurveTrimReconstructionOptions {
    /// Maximum allowed endpoint, domain, and sampled surface disagreement.
    pub tolerance: f64,
    /// Number of samples used for projection and validation.
    pub samples: usize,
    /// Whether existing pcurves may be replaced after validation.
    pub replace_existing: bool,
    /// Whether existing pcurves proven stale by validation may be replaced.
    pub replace_stale_existing: bool,
    /// Whether bounded face surfaces may be re-authored over a wider numeric
    /// UV domain when stale trim reconstruction otherwise falls outside.
    pub expand_bounded_surface_domains: bool,
    /// Whether any declined target rolls back the whole transaction.
    pub require_all_reconstructed: bool,
    /// Whether repaired output must pass full BRep validation.
    pub validate_output: bool,
}

/// Successful curve-backed trim reconstruction result.
#[derive(Debug, Clone)]
pub struct BrepCurveTrimReconstructionResult {
    /// Repaired clone of the input BRep.
    pub brep: BRep,
    /// Structured reconstruction report.
    pub report: BrepCurveTrimReconstructionReport,
}

/// Aggregate curve-backed trim reconstruction report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepCurveTrimReconstructionReport {
    /// Number of surfaced half-edges inspected.
    pub inspected_half_edges: usize,
    /// Number of half-edges selected for reconstruction.
    pub targeted_half_edges: usize,
    /// Number of pcurves reconstructed and committed.
    pub reconstructed_half_edges: usize,
    /// Number of targeted half-edges declined.
    pub declined_half_edges: usize,
    /// Maximum sampled edge-curve/surface disagreement among accepted pcurves.
    pub max_sampled_disagreement: f64,
    /// Per-half-edge decisions.
    pub decisions: Vec<BrepCurveTrimReconstructionDecision>,
}

/// Per-half-edge reconstruction decision.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepCurveTrimReconstructionDecision {
    /// Half-edge considered.
    pub half_edge: HalfEdgeHandle,
    /// Decision outcome.
    pub outcome: BrepCurveTrimReconstructionOutcome,
    /// Maximum sampled edge-curve/surface disagreement.
    pub sampled_disagreement: f64,
}

/// Reconstruction decision outcome.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BrepCurveTrimReconstructionOutcome {
    /// Existing pcurve was present and replacement was disabled.
    RetainedExisting,
    /// Pcurve was reconstructed and accepted by trim-domain validation.
    Reconstructed {
        /// Accepted trim-domain class.
        kind: TrimDomainCandidateKind,
    },
    /// Reconstruction was declined.
    Declined(BrepCurveTrimReconstructionDeclineReason),
}

/// Machine-readable decline reason.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BrepCurveTrimReconstructionDeclineReason {
    /// Required topology was missing or inconsistent.
    MissingTopology,
    /// Half-edge's owning face has no surface.
    MissingSurface,
    /// Edge has no authoritative 3D curve.
    MissingEdgeCurve,
    /// Projection of the 3D curve to the surface UV domain failed.
    ProjectionFailed,
    /// Candidate pcurve failed finite-domain or endpoint validation.
    TrimDomainDeclined(TrimDomainDeclineReason),
    /// Surface-domain expansion could not make all face trims valid.
    DomainExpansionFailed,
    /// Surface evaluated through the pcurve disagreed with the edge curve.
    SampledDisagreement,
}

/// Transactional reconstruction failure.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepCurveTrimReconstructionFailure {
    /// Failure reason.
    pub reason: BrepCurveTrimReconstructionFailureReason,
    /// Report accumulated before rollback.
    pub report: BrepCurveTrimReconstructionReport,
}

/// Reason the reconstruction transaction failed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepCurveTrimReconstructionFailureReason {
    /// At least one targeted half-edge could not be reconstructed.
    UnresolvedTrimReconstruction,
    /// Repaired clone failed validation.
    ValidationFailed(String),
}

/// Reconstruct missing or replaceable pcurves from existing 3D edge curves.
pub fn reconstruct_curve_backed_trims_transactionally(
    brep: &BRep,
    options: BrepCurveTrimReconstructionOptions,
) -> Result<BrepCurveTrimReconstructionResult, Box<BrepCurveTrimReconstructionFailure>> {
    let mut repaired = brep.clone();
    let mut report = BrepCurveTrimReconstructionReport {
        inspected_half_edges: 0,
        targeted_half_edges: 0,
        reconstructed_half_edges: 0,
        declined_half_edges: 0,
        max_sampled_disagreement: 0.0,
        decisions: Vec::new(),
    };
    let targets: Vec<_> = brep.half_edges.keys().collect();
    for half_edge in targets {
        let decision = reconstruct_half_edge_trim(brep, &mut repaired, half_edge, options);
        if !matches!(
            decision.outcome,
            BrepCurveTrimReconstructionOutcome::RetainedExisting
        ) {
            report.targeted_half_edges += 1;
        }
        report.inspected_half_edges += 1;
        match decision.outcome {
            BrepCurveTrimReconstructionOutcome::Reconstructed { .. } => {
                report.reconstructed_half_edges += 1;
                report.max_sampled_disagreement = report
                    .max_sampled_disagreement
                    .max(decision.sampled_disagreement);
            }
            BrepCurveTrimReconstructionOutcome::Declined(_) => report.declined_half_edges += 1,
            BrepCurveTrimReconstructionOutcome::RetainedExisting => {}
        }
        report.decisions.push(decision);
    }
    if options.require_all_reconstructed && report.declined_half_edges > 0 {
        return Err(Box::new(BrepCurveTrimReconstructionFailure {
            reason: BrepCurveTrimReconstructionFailureReason::UnresolvedTrimReconstruction,
            report,
        }));
    }
    if options.validate_output
        && let Err(err) = repaired.validate()
    {
        return Err(Box::new(BrepCurveTrimReconstructionFailure {
            reason: BrepCurveTrimReconstructionFailureReason::ValidationFailed(err.to_string()),
            report,
        }));
    }
    Ok(BrepCurveTrimReconstructionResult {
        brep: repaired,
        report,
    })
}

fn reconstruct_half_edge_trim(
    original: &BRep,
    repaired: &mut BRep,
    half_edge: HalfEdgeHandle,
    options: BrepCurveTrimReconstructionOptions,
) -> BrepCurveTrimReconstructionDecision {
    if original
        .half_edges
        .get(half_edge)
        .is_some_and(|he| he.pcurve.is_some())
        && !options.replace_existing
    {
        match existing_trim_status(original, half_edge, options) {
            ExistingTrimStatus::Clean {
                sampled_disagreement,
            } => {
                return decision(
                    half_edge,
                    BrepCurveTrimReconstructionOutcome::RetainedExisting,
                    sampled_disagreement,
                );
            }
            ExistingTrimStatus::NotApplicable => {
                return decision(
                    half_edge,
                    BrepCurveTrimReconstructionOutcome::RetainedExisting,
                    0.0,
                );
            }
            ExistingTrimStatus::Stale if !options.replace_stale_existing => {
                return decision(
                    half_edge,
                    BrepCurveTrimReconstructionOutcome::RetainedExisting,
                    0.0,
                );
            }
            ExistingTrimStatus::Stale => {}
        }
    }
    let candidate = match build_curve_backed_pcurve(original, half_edge, options) {
        Ok(candidate) => candidate,
        Err(BrepCurveTrimReconstructionDeclineReason::TrimDomainDeclined(
            TrimDomainDeclineReason::OutsideDomain,
        )) if options.expand_bounded_surface_domains => {
            let Some(expansion) =
                expand_repaired_face_domain(original, repaired, half_edge, options)
            else {
                return decision(
                    half_edge,
                    BrepCurveTrimReconstructionOutcome::Declined(
                        BrepCurveTrimReconstructionDeclineReason::DomainExpansionFailed,
                    ),
                    0.0,
                );
            };
            match build_curve_backed_pcurve(repaired, half_edge, options) {
                Ok(candidate) => candidate,
                Err(reason) => {
                    rollback_domain_expansion(repaired, expansion);
                    return decision(
                        half_edge,
                        BrepCurveTrimReconstructionOutcome::Declined(reason),
                        0.0,
                    );
                }
            }
        }
        Err(reason) => {
            return decision(
                half_edge,
                BrepCurveTrimReconstructionOutcome::Declined(reason),
                0.0,
            );
        }
    };
    if let Some(he) = repaired.half_edges.get_mut(half_edge) {
        he.pcurve = Some(candidate.pcurve);
    }
    decision(
        half_edge,
        BrepCurveTrimReconstructionOutcome::Reconstructed {
            kind: candidate.kind,
        },
        candidate.sampled_disagreement,
    )
}

fn rollback_domain_expansion(repaired: &mut BRep, expansion: DomainExpansionCommit) {
    if let Some(face) = repaired.faces.get_mut(expansion.face) {
        face.surface = expansion.previous_surface;
    }
}

pub(super) enum ExistingTrimStatus {
    Clean { sampled_disagreement: f64 },
    Stale,
    NotApplicable,
}

pub(super) fn existing_trim_status(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    options: BrepCurveTrimReconstructionOptions,
) -> ExistingTrimStatus {
    let Some(he) = brep.half_edges.get(half_edge) else {
        return ExistingTrimStatus::NotApplicable;
    };
    let Some(pcurve) = he.pcurve.as_ref() else {
        return ExistingTrimStatus::Stale;
    };
    let Some(face) = brep
        .loops
        .get(he.loop_)
        .and_then(|loop_| brep.faces.get(loop_.face))
    else {
        return ExistingTrimStatus::NotApplicable;
    };
    let Some(surface) = face.surface.as_ref() else {
        return ExistingTrimStatus::NotApplicable;
    };
    let Some(edge_curve) = brep.edges.get(he.edge).and_then(|edge| edge.curve.as_ref()) else {
        return ExistingTrimStatus::NotApplicable;
    };
    let Some(anchors) = half_edge_endpoints(brep, half_edge) else {
        return ExistingTrimStatus::NotApplicable;
    };
    if let Err(decline) = classify_trim_domain_candidate(
        &surface.geometry,
        pcurve,
        Some(anchors),
        options.tolerance,
        options.samples,
    ) {
        return if decline.reason == TrimDomainDeclineReason::NonFiniteDomain {
            ExistingTrimStatus::NotApplicable
        } else {
            ExistingTrimStatus::Stale
        };
    }
    let Some(sampled_disagreement) = sampled_edge_surface_disagreement(
        &edge_curve.geometry,
        (edge_curve.t_start, edge_curve.t_end),
        half_edge,
        brep,
        pcurve,
        options.samples,
    ) else {
        return ExistingTrimStatus::NotApplicable;
    };
    if sampled_disagreement > options.tolerance {
        ExistingTrimStatus::Stale
    } else {
        ExistingTrimStatus::Clean {
            sampled_disagreement,
        }
    }
}

struct AcceptedPcurve {
    pcurve: crate::geometry::pcurve::Pcurve,
    kind: TrimDomainCandidateKind,
    sampled_disagreement: f64,
}

fn build_curve_backed_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    options: BrepCurveTrimReconstructionOptions,
) -> Result<AcceptedPcurve, BrepCurveTrimReconstructionDeclineReason> {
    let he = brep
        .half_edges
        .get(half_edge)
        .ok_or(BrepCurveTrimReconstructionDeclineReason::MissingTopology)?;
    let face = brep
        .loops
        .get(he.loop_)
        .and_then(|loop_| brep.faces.get(loop_.face))
        .ok_or(BrepCurveTrimReconstructionDeclineReason::MissingTopology)?;
    let surface = face
        .surface
        .as_ref()
        .ok_or(BrepCurveTrimReconstructionDeclineReason::MissingSurface)?;
    let edge_curve = brep
        .edges
        .get(he.edge)
        .and_then(|edge| edge.curve.as_ref())
        .ok_or(BrepCurveTrimReconstructionDeclineReason::MissingEdgeCurve)?;
    let pcurve = project_edge_curve_to_surface_pcurve(brep, half_edge, options.samples)
        .ok_or(BrepCurveTrimReconstructionDeclineReason::ProjectionFailed)?;
    let anchors = half_edge_endpoints(brep, half_edge)
        .ok_or(BrepCurveTrimReconstructionDeclineReason::MissingTopology)?;
    let candidate = classify_trim_domain_candidate(
        &surface.geometry,
        &pcurve,
        Some(anchors),
        options.tolerance,
        options.samples,
    )
    .map_err(|decline| {
        BrepCurveTrimReconstructionDeclineReason::TrimDomainDeclined(decline.reason)
    })?;
    let disagreement = sampled_edge_surface_disagreement(
        &edge_curve.geometry,
        (edge_curve.t_start, edge_curve.t_end),
        half_edge,
        brep,
        &pcurve,
        options.samples,
    )
    .ok_or(BrepCurveTrimReconstructionDeclineReason::MissingTopology)?;
    if disagreement > options.tolerance {
        return Err(BrepCurveTrimReconstructionDeclineReason::SampledDisagreement);
    }
    Ok(AcceptedPcurve {
        pcurve,
        kind: candidate.kind,
        sampled_disagreement: disagreement,
    })
}

fn decision(
    half_edge: HalfEdgeHandle,
    outcome: BrepCurveTrimReconstructionOutcome,
    sampled_disagreement: f64,
) -> BrepCurveTrimReconstructionDecision {
    BrepCurveTrimReconstructionDecision {
        half_edge,
        outcome,
        sampled_disagreement,
    }
}

#[cfg(test)]
#[path = "brep_curve_trim_reconstruction_tests.rs"]
mod tests;
