//! Candidate-only multi-patch exact-boundary cap construction.
#![allow(dead_code)]
//!
//! This module starts as a dormant implementation surface. Runtime routing
//! remains unchanged until the orchestration and transaction slices are built.

use std::time::Instant;

use popkern_core::geometry::{NurbsCurve, Point3};

use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSide,
};

#[path = "multipatch_budget.rs"]
pub(crate) mod budget;

#[path = "multipatch_commit.rs"]
mod commit;

#[path = "multipatch_inner.rs"]
mod inner;

#[path = "multipatch_inventory.rs"]
mod inventory;

#[path = "multipatch_report.rs"]
pub(crate) mod report;

#[path = "multipatch_strip.rs"]
pub(crate) mod strip;

#[path = "multipatch_validation.rs"]
mod validation;

use budget::ExactBoundaryMultiPatchBudget;
use inner::{
    ExactBoundaryInnerPolePatchPrototype, ExactBoundaryStripInnerInterfaceCurve,
    construct_inner_pole_patch_prototype, construct_strip_inner_interface_curve,
    measure_inner_interface_g1,
};
use inventory::{
    ExactBoundaryMultiPatchInventory, ExactBoundaryMultiPatchInventoryReport,
    ExactBoundaryOuterStripPatchInventory,
};
use report::{
    ExactBoundaryMultiPatchCandidateOutcome, ExactBoundaryMultiPatchCandidateReport,
    ExactBoundaryMultiPatchDiagnostic, ExactBoundaryMultiPatchPerformanceReport,
};
use strip::{
    ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryStripBoundaryDirection,
    scaffold_outer_strip_boundary, solve_strip_support_normal_row,
};
use validation::{
    ExactBoundaryInnerPoleValidationReport, ExactBoundaryOuterStripValidationReport,
    validate_inner_pole_patch_prototype, validate_outer_strip_candidate,
};

/// Inputs for candidate-only multi-patch orchestration.
pub(crate) struct ExactBoundaryMultiPatchCandidateRequest<'a> {
    /// Source sector being split into candidate patches.
    pub sector: &'a ExactBoundaryCapSectorInput,
    /// Failing outer sides selected by the route classifier.
    pub selected_sides: &'a [ExactBoundaryOuterSide],
    /// Exact source boundary for the first selected side.
    pub boundary: &'a NurbsCurve,
    /// Candidate parameter direction along the copied source boundary.
    pub boundary_direction: ExactBoundaryStripBoundaryDirection,
    /// Deterministic support-ribbon samples for the selected side.
    pub support_samples: &'a [ExactBoundaryRibbonSample],
    /// Candidate inner pole point.
    pub pole: Point3,
    /// Local explicit tolerances and iteration budgets.
    pub budget: ExactBoundaryMultiPatchBudget,
}

/// Candidate-only orchestration result with no live BRep mutation.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryMultiPatchCandidateAttempt {
    /// Ordered candidate report for routing, diagnostics, UI, and API clients.
    pub report: ExactBoundaryMultiPatchCandidateReport,
    /// Commit-eligible candidate payload, present only after every gate passes.
    pub candidate: Option<ExactBoundaryMultiPatchCandidate>,
}

impl ExactBoundaryMultiPatchCandidateAttempt {
    fn declined(report: ExactBoundaryMultiPatchCandidateReport) -> Self {
        debug_assert_eq!(
            report.outcome,
            ExactBoundaryMultiPatchCandidateOutcome::Declined
        );
        Self {
            report,
            candidate: None,
        }
    }

    fn commit_eligible(
        report: ExactBoundaryMultiPatchCandidateReport,
        candidate: ExactBoundaryMultiPatchCandidate,
    ) -> Self {
        debug_assert_eq!(
            report.outcome,
            ExactBoundaryMultiPatchCandidateOutcome::CommitEligible
        );
        Self {
            report,
            candidate: Some(candidate),
        }
    }
}

/// Fully validated in-memory candidate passed to later commit preflight slices.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundaryMultiPatchCandidate {
    /// Candidate ownership inventory.
    pub inventory: ExactBoundaryMultiPatchInventory,
    /// Validated inventory count report.
    pub inventory_report: ExactBoundaryMultiPatchInventoryReport,
    /// Source-backed outer strip ownership.
    pub strip_inventory: ExactBoundaryOuterStripPatchInventory,
    /// Solved outer strip candidate.
    pub strip_candidate: ExactBoundaryOuterStripBoundaryScaffold,
    /// Strip validation evidence.
    pub strip_validation: ExactBoundaryOuterStripValidationReport,
    /// Candidate strip/inner interface curve.
    pub interface: ExactBoundaryStripInnerInterfaceCurve,
    /// G1-validated inner pole patch candidate.
    pub inner_patch: ExactBoundaryInnerPolePatchPrototype,
    /// Inner-patch validation evidence.
    pub inner_validation: ExactBoundaryInnerPoleValidationReport,
}

/// Build a candidate-only multi-patch sector split behind the dormant facade.
pub(crate) fn prepare_multipatch_candidate(
    request: ExactBoundaryMultiPatchCandidateRequest<'_>,
) -> ExactBoundaryMultiPatchCandidateAttempt {
    let started_at = Instant::now();
    let mut performance = ExactBoundaryMultiPatchPerformanceReport {
        sector_count: 1,
        ..ExactBoundaryMultiPatchPerformanceReport::default()
    };

    if request.selected_sides.len() > 1 {
        let mut report =
            ExactBoundaryMultiPatchCandidateReport::declined_for_two_outer_sides(request.sector);
        report.performance.sector_count = performance.sector_count;
        report.performance.elapsed_micros = started_at.elapsed().as_micros();
        return ExactBoundaryMultiPatchCandidateAttempt::declined(report);
    }

    let Some(side) = request.selected_sides.first().copied() else {
        let mut report =
            ExactBoundaryMultiPatchCandidateReport::declined_for_two_outer_sides(request.sector);
        report.performance.sector_count = performance.sector_count;
        report.performance.elapsed_micros = started_at.elapsed().as_micros();
        return ExactBoundaryMultiPatchCandidateAttempt::declined(report);
    };

    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(request.sector, side);
    performance.outer_strip_patch_count = inventory.outer_strips.len();
    performance.interface_count = inventory.interfaces.len();
    performance.inner_patch_count = usize::from(inventory.inner_patch.is_some());

    let inventory_report = match inventory.validate_against_sector(request.sector) {
        Ok(report) => report,
        Err(diagnostic) => {
            let mut report =
                ExactBoundaryMultiPatchCandidateReport::declined_for_inventory_mismatch(diagnostic);
            report.performance = finish_performance(performance, started_at);
            return ExactBoundaryMultiPatchCandidateAttempt::declined(report);
        }
    };
    let strip = &inventory.outer_strips[0];
    let strip_inventory = strip.clone();
    let interface_inventory = &inventory.interfaces[0];

    let scaffold = match scaffold_outer_strip_boundary(
        request.sector,
        strip,
        request.boundary,
        request.boundary_direction,
        request.budget.outer_boundary_residual,
    ) {
        Ok(scaffold) => scaffold,
        Err(diagnostic) => {
            return declined_with_context(
                request.sector,
                Some(inventory_report),
                diagnostic,
                performance,
                started_at,
            );
        }
    };
    let strip_candidate = match solve_strip_support_normal_row(
        request.sector,
        strip,
        &scaffold,
        request.boundary,
        request.support_samples,
        request.budget,
    ) {
        Ok(candidate) => candidate,
        Err(diagnostic) => {
            return declined_with_context(
                request.sector,
                Some(inventory_report),
                diagnostic,
                performance,
                started_at,
            );
        }
    };
    if let Some(report) = strip_candidate.support_normal_report {
        performance.support_strip_sample_count += report.support_strip_sample_count;
        performance.bounded_iteration_count += report.bounded_iteration_count;
    }

    let strip_validation = match validate_outer_strip_candidate(
        request.sector,
        strip,
        None,
        &strip_candidate,
        request.boundary,
        request.support_samples,
        request.budget,
    ) {
        Ok(report) => report,
        Err(diagnostic) => {
            return declined_with_context(
                request.sector,
                Some(inventory_report),
                diagnostic,
                performance,
                started_at,
            );
        }
    };
    performance.support_strip_sample_count += strip_validation.support_strip_sample_count;

    let interface = match construct_strip_inner_interface_curve(
        request.sector,
        strip,
        interface_inventory,
        &strip_candidate,
        request.budget.interface_g0,
    ) {
        Ok(interface) => interface,
        Err(diagnostic) => {
            return declined_with_context(
                request.sector,
                Some(inventory_report),
                diagnostic,
                performance,
                started_at,
            );
        }
    };
    performance.interface_sample_count += interface.sample_count;

    let prototype = match construct_inner_pole_patch_prototype(
        request.sector,
        strip,
        interface_inventory,
        &interface,
        request.pole,
        request.budget.pole_g0,
    ) {
        Ok(prototype) => prototype,
        Err(diagnostic) => {
            return declined_with_context(
                request.sector,
                Some(inventory_report),
                diagnostic,
                performance,
                started_at,
            );
        }
    };
    let inner_patch = match measure_inner_interface_g1(
        request.sector,
        strip,
        interface_inventory,
        &strip_candidate,
        &interface,
        &prototype,
        request.budget.interface_g1_deficit,
    ) {
        Ok(prototype) => prototype,
        Err(diagnostic) => {
            return declined_with_context(
                request.sector,
                Some(inventory_report),
                diagnostic,
                performance,
                started_at,
            );
        }
    };
    if let Some(report) = inner_patch.interface_g1_report {
        performance.interface_sample_count += report.sample_count;
    }

    let inner_validation = match validate_inner_pole_patch_prototype(
        request.sector,
        strip,
        interface_inventory,
        &interface,
        &inner_patch,
        request.budget,
    ) {
        Ok(report) => report,
        Err(diagnostic) => {
            return declined_with_context(
                request.sector,
                Some(inventory_report),
                diagnostic,
                performance,
                started_at,
            );
        }
    };

    let report = ExactBoundaryMultiPatchCandidateReport::commit_eligible(
        request.sector,
        inventory_report.clone(),
        finish_performance(performance, started_at),
    );
    ExactBoundaryMultiPatchCandidateAttempt::commit_eligible(
        report,
        ExactBoundaryMultiPatchCandidate {
            inventory,
            inventory_report,
            strip_inventory,
            strip_candidate,
            strip_validation,
            interface,
            inner_patch,
            inner_validation,
        },
    )
}

fn declined_with_context(
    sector: &ExactBoundaryCapSectorInput,
    inventory: Option<ExactBoundaryMultiPatchInventoryReport>,
    diagnostic: ExactBoundaryMultiPatchDiagnostic,
    performance: ExactBoundaryMultiPatchPerformanceReport,
    started_at: Instant,
) -> ExactBoundaryMultiPatchCandidateAttempt {
    ExactBoundaryMultiPatchCandidateAttempt::declined(
        ExactBoundaryMultiPatchCandidateReport::declined_with_context(
            sector,
            inventory,
            diagnostic,
            finish_performance(performance, started_at),
        ),
    )
}

fn finish_performance(
    mut performance: ExactBoundaryMultiPatchPerformanceReport,
    started_at: Instant,
) -> ExactBoundaryMultiPatchPerformanceReport {
    performance.elapsed_micros = started_at.elapsed().as_micros();
    performance
}

#[cfg(test)]
#[path = "multipatch_tests.rs"]
mod tests;
