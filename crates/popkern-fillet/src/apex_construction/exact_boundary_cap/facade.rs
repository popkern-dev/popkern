//! Semantic exact-boundary cap candidate facade.

use std::time::Instant;

use popkern_core::topology::{BRep, FaceHandle};

use crate::apex_construction::error::ApexBuildError;
use crate::apex_construction::gregory_cap::GregoryCapSector;
use crate::apex_construction::plan::ApexPlan;
use crate::types::{
    ApexConstructionDiagnostic, ExactBoundaryCandidatePerformanceReport, ExactBoundaryCapReport,
    ExactBoundaryRouteDecisionStep,
};

#[path = "facade_multipatch_candidate_continuity.rs"]
mod candidate_continuity;
#[path = "facade_multipatch_commit.rs"]
mod multipatch_commit;
#[path = "facade_multipatch_route.rs"]
mod multipatch_route;
#[path = "facade_report.rs"]
mod report;

use super::boundary::{ExactBoundaryCurveBudget, prepare_exact_boundary_curves};
use super::cap_quality_postcondition::worst_mid_seam_misalignment;
use super::commit::{commit_exact_boundary_sector_surfaces, post_commit_validation_diagnostic};
use super::corner::{ExactBoundaryInnerSpokeBudget, prepare_exact_boundary_inner_spokes};
use super::cross_spoke_g1::{cross_spoke_g1_enabled, enforce_cross_spoke_g1};
use super::facts::attach_shared_repair_fact_ids;
use super::pcurves::{ExactBoundaryPcurveBudget, prepare_exact_boundary_pcurves};
use super::performance::{
    finalize_performance_report, record_pcurve_report, record_support_family_counts,
    record_validation_report, support_sample_count, surface_control_point_count,
};
use super::pns_candidate::{
    pns_cap_substitution_enabled, pns_substitution_note, substitute_pns_sector_surfaces,
};
use super::ribbon::{ExactBoundaryRibbonBudget, prepare_exact_boundary_ribbons};
use super::sector::{ExactBoundaryCapSectorInput, extract_exact_boundary_cap_sectors};
use super::solve::{
    ExactBoundarySolveBudget, refine_support_normal_strips_with_pcurves,
    solve_exact_boundary_surfaces,
};
use super::topology::materialize_exact_boundary_sector_topology;
use super::validation::{
    ExactBoundaryValidationBudget, ExactBoundaryValidationInputs,
    collect_constrained_strip_band_diagnostics, validate_exact_boundary_candidate,
};
use multipatch_commit::{
    ExactBoundaryMultipatchCommitOutcome, commit_multipatch_routes_transactionally,
};
use multipatch_route::{
    attempt_multipatches_after_one_patch_decline, collect_live_one_patch_cap_support_failures,
};
pub(super) use report::build_exact_boundary_report;

/// Status of one exact-boundary cap candidate attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryCapCandidateStatus {
    /// Candidate passed validation and committed exact-boundary sector faces.
    Accepted,
    /// Candidate declined and restored the input BRep snapshot.
    Declined,
}

/// Reusable exact-boundary candidate result for construction, diagnostics, and UX.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryCapCandidateResult {
    /// Candidate status after validation/commit or rollback.
    pub status: ExactBoundaryCapCandidateStatus,
    /// Sector ownership records. On accept these reference committed sector faces.
    pub sectors: Vec<ExactBoundaryCapSectorInput>,
    /// Faces committed by the accepted candidate.
    pub committed_faces: Vec<FaceHandle>,
    /// Candidate-only performance counters and budget evidence.
    pub performance: ExactBoundaryCandidatePerformanceReport,
    /// Versioned report for UI and programmatic consumers.
    pub report: ExactBoundaryCapReport,
    /// Structured diagnostics emitted during decline.
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}

impl ExactBoundaryCapCandidateResult {
    /// Number of exact curves prepared without approximation.
    pub(crate) fn prepared_curve_count(&self) -> usize {
        self.performance.prepared_curve_count
    }

    /// Number of generated inner spokes prepared without topology mutation.
    pub(crate) fn inner_spoke_count(&self) -> usize {
        self.performance.inner_spoke_count
    }

    /// Number of support derivative ribbons prepared from support pcurve evidence.
    pub(crate) fn ribbon_count(&self) -> usize {
        self.performance.support_ribbon_count
    }

    /// Number of in-memory sector surfaces solved before BRep commit integration.
    pub(crate) fn solved_surface_count(&self) -> usize {
        self.performance.solved_surface_count
    }

    /// Number of cap/support pcurve pairs prepared before BRep commit integration.
    pub(crate) fn pcurve_pair_count(&self) -> usize {
        self.performance.pcurve_pair_count
    }
}

/// Run the exact-boundary cap candidate transactionally.
///
/// The function mutates `brep` only when all exact-boundary postconditions pass.
/// Any diagnostic or incomplete pcurve inventory restores the input snapshot and
/// returns a declined candidate result for the caller's containment path.
pub(crate) fn run_exact_boundary_cap_candidate(
    brep: &mut BRep,
    cap_face: FaceHandle,
    plan: &ApexPlan,
    sectors: &[GregoryCapSector],
) -> Result<ExactBoundaryCapCandidateResult, ApexBuildError> {
    let started = Instant::now();
    let snapshot = brep.clone();
    let baseline_validation_errors = snapshot.validate_collect_fast();
    let mut candidate_sectors = extract_exact_boundary_cap_sectors(brep, cap_face, plan, sectors)
        .map_err(|diagnostic| diagnostic.into_build_error())?;
    let mut performance = ExactBoundaryCandidatePerformanceReport {
        sector_count: candidate_sectors.len(),
        ..Default::default()
    };
    let mut diagnostics = Vec::new();
    // PK-201: diagnostics above are FAILURE EVIDENCE — `accepted`
    // below requires the vector to be empty, so anything pushed there
    // vetoes the route. `informational` is the channel for everything
    // a route wants to SAY without changing what it DOES. Without the
    // split, explaining yourself costs you the corner, which is how a
    // silent folded cap came to beat a fallback that measures 80x
    // better on the same rim.
    let mut informational: Vec<crate::types::ApexConstructionDiagnostic> = Vec::new();
    let boundary_curves = match prepare_exact_boundary_curves(
        brep,
        &candidate_sectors,
        ExactBoundaryCurveBudget::default(),
    ) {
        Ok(prepared) => prepared,
        Err(diagnostic) => {
            diagnostics.push(diagnostic.into_apex_diagnostic());
            *brep = snapshot;
            // Informational records join the report AFTER the accept decision,
            // so they can never change which route serves a corner.
            let mut diagnostics = diagnostics;
            diagnostics.extend(informational);
            attach_shared_repair_fact_ids(brep, &mut diagnostics);
            finalize_performance_report(&mut performance, started, &diagnostics);
            let report = build_exact_boundary_report(
                brep,
                ExactBoundaryCapCandidateStatus::Declined,
                &candidate_sectors,
                &[],
                &performance,
                &diagnostics,
            );
            return Ok(ExactBoundaryCapCandidateResult {
                status: ExactBoundaryCapCandidateStatus::Declined,
                sectors: candidate_sectors,
                committed_faces: Vec::new(),
                performance,
                diagnostics,
                report,
            });
        }
    };
    performance.prepared_curve_count = boundary_curves.curves.len();

    let inner_spokes = match prepare_exact_boundary_inner_spokes(
        &candidate_sectors,
        &boundary_curves,
        ExactBoundaryInnerSpokeBudget::default(),
    ) {
        Ok(prepared) => Some(prepared),
        Err(diagnostic) => {
            diagnostics.push(diagnostic.into_apex_diagnostic());
            None
        }
    };
    performance.inner_spoke_count = inner_spokes
        .as_ref()
        .map(|prepared| prepared.sectors.len() * 2)
        .unwrap_or(0);

    let topology = match inner_spokes.as_ref() {
        Some(inner_spokes) => match materialize_exact_boundary_sector_topology(
            brep,
            cap_face,
            &candidate_sectors,
            inner_spokes,
        ) {
            Ok(topology) => {
                candidate_sectors = topology.sectors.clone();
                Some(topology)
            }
            Err(diagnostic) => {
                diagnostics.push(diagnostic.into_apex_diagnostic());
                None
            }
        },
        None => None,
    };

    let ribbons = match topology.as_ref() {
        Some(_) => match prepare_exact_boundary_ribbons(
            brep,
            &candidate_sectors,
            &boundary_curves,
            ExactBoundaryRibbonBudget::default(),
        ) {
            Ok(prepared) => Some(prepared),
            Err(diagnostic) => {
                diagnostics.push(diagnostic.into_apex_diagnostic());
                None
            }
        },
        None => None,
    };
    performance.support_ribbon_count = ribbons
        .as_ref()
        .map(|prepared| prepared.ribbons.len())
        .unwrap_or(0);
    performance.support_sample_count = ribbons.as_ref().map(support_sample_count).unwrap_or(0);
    if let Some(ribbons) = ribbons.as_ref() {
        record_support_family_counts(&mut performance, ribbons);
    }

    let mut solved_surfaces = match (inner_spokes.as_ref(), ribbons.as_ref()) {
        (Some(inner_spokes), Some(ribbons)) => {
            match solve_exact_boundary_surfaces(
                inner_spokes,
                ribbons,
                ExactBoundarySolveBudget::default(),
            ) {
                Ok(solved) => Some(solved),
                Err(diagnostic) => {
                    diagnostics.push(diagnostic.into_apex_diagnostic());
                    None
                }
            }
        }
        _ => None,
    };
    // PK-192: substitute PnS bi-4 surfaces over the same rim before
    // validation, so they face the identical acceptance gates the
    // Coons wedges did. The fan is G0 across its own internal seams
    // (measured 89.9 deg worst at n=3, 89.6 at n=4) while the PnS cap
    // is G1 there by construction (0.0 deg). On failure the Coons
    // surfaces stay and behaviour is byte-identical to before.
    //
    // Deliberately silent: `accepted` at the end of this function is
    // `committed && diagnostics.is_empty() && ...`, so ANY recorded
    // diagnostic — even a purely informational one — vetoes the whole
    // exact-boundary route and drops the corner to its fallback. That
    // coupling is tracked as PK-201; until it is fixed, reporting here
    // would silently change which route serves every corner.
    if pns_cap_substitution_enabled()
        && let Some(solved) = solved_surfaces.as_mut()
    {
        let outcome = substitute_pns_sector_surfaces(
            brep,
            &candidate_sectors,
            solved,
            brep.tolerance.linear(),
        );
        informational.push(pns_substitution_note(candidate_sectors.len(), &outcome));
    }
    // PK-201 construction half: make adjacent sectors G1 across their
    // shared spokes. The wedges were only ever constrained against
    // their OUTER boundaries, so nothing held the inner spokes and
    // every measured cap was creased there. Default on, measured —
    // see `cross_spoke_g1`.
    if cross_spoke_g1_enabled()
        && let Some(solved) = solved_surfaces.as_mut()
    {
        let outcome = enforce_cross_spoke_g1(solved, brep.tolerance.linear());
        informational.push(cross_spoke_note(candidate_sectors.len(), &outcome));
    }
    performance.solved_surface_count = solved_surfaces
        .as_ref()
        .map(|solved| solved.validated_surface_count())
        .unwrap_or(0);
    performance.surface_control_point_count = solved_surfaces
        .as_ref()
        .map(surface_control_point_count)
        .unwrap_or(0);

    let mut committed = false;
    let mut committed_faces = Vec::new();
    let mut decision_chain = Vec::new();
    if let (Some(topology), Some(inner_spokes), Some(ribbons), Some(solved_surfaces)) = (
        topology.as_ref(),
        inner_spokes.as_ref(),
        ribbons.as_ref(),
        solved_surfaces.as_mut(),
    ) {
        match prepare_exact_boundary_pcurves(
            inner_spokes,
            ribbons,
            solved_surfaces,
            ExactBoundaryPcurveBudget::default(),
        ) {
            Ok(pcurves) => {
                refine_support_normal_strips_with_pcurves(solved_surfaces, &pcurves, ribbons);
                record_pcurve_report(&mut performance, &pcurves);
                let _topology_evidence = (topology.pole_vertex, topology.spoke_edges.len());
                match validate_exact_boundary_candidate(ExactBoundaryValidationInputs {
                    brep,
                    topology,
                    sectors: &candidate_sectors,
                    spokes: inner_spokes,
                    ribbons,
                    solved: solved_surfaces,
                    pcurves: &pcurves,
                    budget: ExactBoundaryValidationBudget::default(),
                }) {
                    Ok(validation_report) => {
                        record_validation_report(&mut performance, validation_report);
                        committed_faces = candidate_sectors
                            .iter()
                            .map(|sector| sector.cap_face)
                            .collect();
                        if let Err(diagnostic) = commit_exact_boundary_sector_surfaces(
                            brep,
                            cap_face,
                            &candidate_sectors,
                            solved_surfaces,
                            &pcurves,
                        ) {
                            diagnostics.push(diagnostic.into_apex_diagnostic());
                        } else if let Some(diagnostic) = post_commit_validation_diagnostic(
                            &baseline_validation_errors,
                            &brep.validate_collect_fast(),
                            brep.tolerance.linear(),
                        ) {
                            diagnostics.push(diagnostic.into_apex_diagnostic());
                        } else {
                            // PK-201 measurement (behaviour-neutral):
                            // report how continuous the committed cap
                            // actually is across its own seams. Making
                            // this a *gate* at the crate's 5e-2 G1 bar
                            // was tried and rejected 20+ exact-boundary
                            // acceptances across the direct-apex matrix
                            // — see the ticket. Recording the number is
                            // what lets that scope be sized honestly.
                            if std::env::var_os("POPKERN_CAP_SEAM_REPORT").is_some()
                                && let Some(worst) =
                                    worst_mid_seam_misalignment(brep, &committed_faces)
                            {
                                eprintln!(
                                    "[cap-seam] faces={} worst_mid_seam_normal={worst:.3e}",
                                    committed_faces.len()
                                );
                            }
                            committed = true;
                        }
                    }
                    Err(diagnostic) => {
                        let mut repair_diagnostics = collect_constrained_strip_band_diagnostics(
                            solved_surfaces,
                            &pcurves,
                            ribbons,
                            ExactBoundaryValidationBudget::default(),
                        );
                        if repair_diagnostics.is_empty() {
                            repair_diagnostics.push(diagnostic.clone());
                        }
                        let live_diagnostics = collect_live_one_patch_cap_support_failures(
                            brep,
                            cap_face,
                            &candidate_sectors,
                            solved_surfaces,
                            &pcurves,
                        );
                        for live_diagnostic in live_diagnostics {
                            if !repair_diagnostics
                                .iter()
                                .any(|existing| existing.sector_name == live_diagnostic.sector_name)
                            {
                                repair_diagnostics.push(live_diagnostic);
                            }
                        }
                        let multipatch_attempts = attempt_multipatches_after_one_patch_decline(
                            &repair_diagnostics,
                            &candidate_sectors,
                            &boundary_curves,
                            ribbons,
                            &pcurves,
                            inner_spokes,
                            &mut diagnostics,
                        );
                        let multipatch_candidates = multipatch_attempts
                            .iter()
                            .filter_map(|attempt| {
                                attempt
                                    .candidate
                                    .as_ref()
                                    .map(|candidate| (attempt, candidate))
                            })
                            .collect::<Vec<_>>();
                        if !multipatch_attempts.is_empty()
                            && multipatch_candidates.len() == multipatch_attempts.len()
                        {
                            decision_chain.push(ExactBoundaryRouteDecisionStep::declined(
                                "exact_boundary.one_patch",
                                Some(diagnostic.code.to_string()),
                            ));
                            match commit_multipatch_routes_transactionally(
                                brep,
                                cap_face,
                                &candidate_sectors,
                                solved_surfaces,
                                &pcurves,
                                &baseline_validation_errors,
                                &multipatch_candidates,
                            ) {
                                ExactBoundaryMultipatchCommitOutcome::Accepted { faces } => {
                                    committed_faces = faces;
                                    committed = true;
                                    decision_chain.push(ExactBoundaryRouteDecisionStep::accepted(
                                        "exact_boundary.multipatch",
                                    ));
                                }
                                ExactBoundaryMultipatchCommitOutcome::Declined {
                                    step,
                                    diagnostic,
                                } => {
                                    decision_chain.push(step);
                                    diagnostics.push(diagnostic);
                                }
                            }
                        } else {
                            decision_chain.push(ExactBoundaryRouteDecisionStep::declined(
                                "exact_boundary.one_patch",
                                Some(diagnostic.code.to_string()),
                            ));
                            diagnostics.push(diagnostic.clone().into_apex_diagnostic());
                            for attempt in multipatch_attempts {
                                decision_chain.push(attempt.step);
                            }
                        }
                    }
                }
            }
            Err(diagnostic) => diagnostics.push(diagnostic.into_apex_diagnostic()),
        }
    }

    let accepted = committed && diagnostics.is_empty() && performance.pcurve_pair_count > 0;
    if !accepted {
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace exact-boundary/facade] rollback: committed={} diagnostics={} pcurve_pairs={} first_diag={:?}",
                committed,
                diagnostics.len(),
                performance.pcurve_pair_count,
                diagnostics.first().map(|d| format!("{d:?}")),
            );
        }
        *brep = snapshot;
        committed_faces.clear();
    }
    attach_shared_repair_fact_ids(brep, &mut diagnostics);
    finalize_performance_report(&mut performance, started, &diagnostics);
    let status = if accepted {
        ExactBoundaryCapCandidateStatus::Accepted
    } else {
        ExactBoundaryCapCandidateStatus::Declined
    };
    let report = build_exact_boundary_report(
        brep,
        status,
        &candidate_sectors,
        &committed_faces,
        &performance,
        &diagnostics,
    );
    let report = if decision_chain.is_empty() {
        report
    } else {
        report.with_decision_chain(decision_chain)
    };

    Ok(ExactBoundaryCapCandidateResult {
        status,
        sectors: candidate_sectors,
        committed_faces,
        performance,
        report,
        diagnostics,
    })
}

/// Informational record of the cross-spoke correction.
fn cross_spoke_note(
    n_sectors: usize,
    outcome: &super::cross_spoke_g1::CrossSpokeOutcome,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "exact_boundary_cap.cross_spoke_g1".to_string(),
            code: "exact_boundary_cap.cross_spoke_g1".to_string(),
            detail_code: Some(format!(
                "spokes={};corrected_control_points={};dimension_mismatched_neighbours={}",
                outcome.spokes, outcome.corrected, outcome.dimension_mismatched
            )),
            sector_name: format!("n{n_sectors}_corner_cap"),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: None,
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
    )
}
