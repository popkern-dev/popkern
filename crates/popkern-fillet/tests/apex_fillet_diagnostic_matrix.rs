//! Apex/fillet diagnostic classification matrix.
//!
//! This test is intentionally lightweight: it locks the structured manifest
//! and counters that live fixture rows must satisfy as diagnostics mature.

#[path = "support/apex_fillet_diagnostic_matrix_manifest.rs"]
mod apex_fillet_diagnostic_matrix_manifest;

use apex_fillet_diagnostic_matrix_manifest::{
    DiagnosticAxis, DiagnosticMatrixTerminal, FixtureTier, UnknownDiagnosticFacet,
    corner_apex_rows, general_kernel_rows, seed_matrix_rows,
};
use popkern_core::euler::corner_apex_condition_ids;

#[test]
fn seed_matrix_has_no_unknown_candidate_families() {
    let summary = MatrixSummary::from_rows(seed_matrix_rows());
    assert_eq!(
        summary.unknown_family_rows,
        0,
        "unknown candidate family rows must be partitioned before mutator broadening: {}",
        summary.diagnostic_label()
    );
    assert_eq!(
        summary.unpartitioned_unknowns,
        0,
        "all unknowns must have a closure category and task-backed path: {}",
        summary.diagnostic_label()
    );
}

#[test]
fn general_kernel_rows_are_ready_for_broad_coverage() {
    for row in general_kernel_rows() {
        assert_eq!(
            row.fixture_tier,
            FixtureTier::SyntheticAxis,
            "{} general-kernel row must be a synthetic axis fixture: {}",
            row.case_id,
            row.diagnostic_label()
        );
        assert_complete_six_facets(*row);
    }
}

#[test]
fn matrix_covers_every_declared_diagnostic_axis() {
    let corner_rows = corner_apex_rows();
    let covered = seed_matrix_rows()
        .iter()
        .chain(general_kernel_rows().iter())
        .flat_map(|row| row.covered_axes.iter().copied())
        .chain(
            corner_rows
                .iter()
                .flat_map(|row| row.covered_axes.iter().copied()),
        )
        .collect::<std::collections::BTreeSet<_>>();
    for axis in DiagnosticAxis::all() {
        assert!(
            covered.contains(axis),
            "matrix is missing declared diagnostic axis {}",
            axis.diagnostic_label()
        );
    }
}

#[test]
fn corner_apex_rows_cover_every_declared_condition() {
    let rows = corner_apex_rows();
    let covered = rows
        .iter()
        .map(|row| row.condition_id)
        .collect::<std::collections::BTreeSet<_>>();
    for condition in corner_apex_condition_ids() {
        let label = condition.diagnostic_label();
        assert!(
            covered.contains(label.as_str()),
            "corner/apex matrix is missing condition {label}"
        );
    }
}

#[test]
fn corner_apex_rows_are_implementation_ready() {
    for row in corner_apex_rows() {
        let label = row.diagnostic_label();
        assert_eq!(
            row.fixture_tier,
            FixtureTier::SyntheticAxis,
            "{} corner/apex row must stay diagnostic-only synthetic coverage: {label}",
            row.condition_id
        );
        assert!(
            row.condition_id.starts_with("corner_"),
            "{} condition id is not scoped: {label}",
            row.condition_id
        );
        assert!(
            !row.covered_axes.is_empty(),
            "{} lacks axes: {label}",
            row.condition_id
        );
        assert_not_generic_blocker(row.expected_reason, row.condition_id);
        assert_not_generic_blocker(row.feasibility, row.condition_id);
        assert_ne!(
            row.topology_resource, "none",
            "{} lacks topology facet: {label}",
            row.condition_id
        );
        assert_ne!(
            row.geometry_feasibility, "none",
            "{} lacks geometry facet: {label}",
            row.condition_id
        );
        assert_ne!(
            row.validity_repair, "none",
            "{} lacks validity facet: {label}",
            row.condition_id
        );
        assert_ne!(
            row.operation_dependency, "none",
            "{} lacks dependency facet: {label}",
            row.condition_id
        );
        assert_ne!(
            row.semantic_history, "none",
            "{} lacks semantic history: {label}",
            row.condition_id
        );
        assert_ne!(
            row.downstream_readiness, "none",
            "{} lacks downstream facet: {label}",
            row.condition_id
        );
        assert_ne!(
            row.tolerance_source, "none",
            "{} lacks tolerance source: {label}",
            row.condition_id
        );
        if row.terminal == DiagnosticMatrixTerminal::NamedMissingPrimitive {
            assert!(
                row.operation_dependency.contains("missing-"),
                "{} missing primitive row does not name a missing dependency: {label}",
                row.condition_id
            );
        }
        if row.terminal == DiagnosticMatrixTerminal::AmbiguousDependency {
            assert!(
                row.operation_dependency.contains("ambiguous"),
                "{} ambiguous row does not name ambiguous dependency: {label}",
                row.condition_id
            );
            assert!(
                row.closure_task.is_some(),
                "{} ambiguous row lacks closure task: {label}",
                row.condition_id
            );
        }
    }
}

#[test]
fn corner_apex_rows_include_executable_invalid_missing_and_ambiguous_cases() {
    let rows = corner_apex_rows();
    let terminals = rows
        .iter()
        .map(|row| row.terminal)
        .collect::<std::collections::BTreeSet<_>>();
    for required in [
        DiagnosticMatrixTerminal::ExecutableCandidate,
        DiagnosticMatrixTerminal::NamedMissingPrimitive,
        DiagnosticMatrixTerminal::AmbiguousDependency,
        DiagnosticMatrixTerminal::InvalidTopologyReason,
    ] {
        assert!(
            terminals.contains(&required),
            "corner/apex rows do not cover terminal {}",
            required.diagnostic_label()
        );
    }
    assert!(
        rows.iter()
            .any(|row| row.downstream_readiness.contains("risk")),
        "corner/apex rows need at least one downstream-risk case"
    );
    assert!(
        rows.iter()
            .any(|row| row.semantic_history.contains("generated")),
        "corner/apex rows need generated semantic-history coverage"
    );
}

#[test]
fn matrix_rows_do_not_use_generic_blockers_as_terminal_reasons() {
    for row in seed_matrix_rows()
        .iter()
        .chain(general_kernel_rows().iter())
    {
        assert_not_generic_blocker(row.expected_reason, row.case_id);
        assert_not_generic_blocker(row.feasibility, row.case_id);
        if row.terminal == DiagnosticMatrixTerminal::NamedMissingPrimitive {
            assert!(
                row.operation_dependency
                    .is_some_and(|dependency| dependency.contains("missing-")),
                "{} named missing primitive row lacks a named missing dependency: {}",
                row.case_id,
                row.diagnostic_label()
            );
        }
        if row.unknown_facet.is_some() {
            assert!(
                row.closure_task.is_some(),
                "{} unknown facet lacks a closure task: {}",
                row.case_id,
                row.diagnostic_label()
            );
        }
    }
}

#[test]
fn seed_matrix_panic_output_names_family_and_feasibility() {
    for row in seed_matrix_rows() {
        let label = row.diagnostic_label();
        assert!(
            label.contains("candidate_family="),
            "{} lacks candidate family: {label}",
            row.case_id
        );
        assert!(
            label.contains("feasibility="),
            "{} lacks feasibility: {label}",
            row.case_id
        );
        assert_complete_six_facets(*row);
        assert!(
            label.contains("terminal="),
            "{} lacks terminal classification: {label}",
            row.case_id
        );
    }
}

fn assert_not_generic_blocker(label: &str, case_id: &str) {
    for generic in [
        "missing_operator",
        "geometry_failed",
        "invalid_topology",
        "unsupported_case",
    ] {
        assert_ne!(
            label, generic,
            "{case_id} uses generic blocker label {generic}"
        );
    }
}

fn assert_complete_six_facets(row: apex_fillet_diagnostic_matrix_manifest::DiagnosticMatrixRow) {
    let label = row.diagnostic_label();
    if row.geometry_feasibility_required {
        assert!(
            row.geometry_feasibility.is_some(),
            "{} requires geometry feasibility but has no geometry facet: {label}",
            row.case_id
        );
    }
    if row.validity_repair_required {
        assert!(
            row.validity_repair.is_some(),
            "{} requires validity/repair readiness but has no validity facet: {label}",
            row.case_id
        );
    }
    if row.operation_dependency_required {
        assert!(
            row.operation_dependency.is_some(),
            "{} requires operation dependency readiness but has no dependency facet: {label}",
            row.case_id
        );
    }
    if row.downstream_readiness_required {
        assert!(
            row.downstream_readiness.is_some(),
            "{} requires downstream readiness but has no downstream facet: {label}",
            row.case_id
        );
    }
}

#[test]
fn matrix_categories_have_stable_labels() {
    for terminal in DiagnosticMatrixTerminal::all() {
        assert!(!terminal.diagnostic_label().is_empty());
    }
    for unknown in UnknownDiagnosticFacet::all() {
        assert!(!unknown.diagnostic_label().is_empty());
    }
}

#[test]
fn seed_matrix_names_current_regression_cases() {
    let cases = seed_matrix_rows()
        .iter()
        .map(|row| row.case_id)
        .collect::<std::collections::BTreeSet<_>>();
    for expected in [
        "cube_mixed_native_direct_apex",
        "cube_extreme_native_direct_apex",
        "pyramid3_uniform_native_direct_apex",
        "pyramid4_uniform_native_direct_apex",
        "pyramid4_mixed_native_direct_apex",
        "gregory_step_roundtrip_tolerance",
        "cube_extreme_40_legacy_route",
        "cube_extreme_50_legacy_route",
        "cube_extreme_100_legacy_route",
        "cube_tiny_radii_legacy_route",
    ] {
        assert!(
            cases.contains(expected),
            "seed matrix is missing regression case {expected}"
        );
    }
}

#[test]
fn seed_matrix_names_apex_semantic_invariants() {
    let cases = seed_matrix_rows()
        .iter()
        .map(|row| row.case_id)
        .collect::<std::collections::BTreeSet<_>>();
    for expected in [
        "apex_retired_vertex_absence",
        "apex_cap_cylinder_seam_cardinality",
        "apex_support_family_ownership",
        "apex_generated_boundary_replacement",
        "apex_cap_transition_continuity",
    ] {
        assert!(
            cases.contains(expected),
            "seed matrix is missing apex semantic invariant {expected}"
        );
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
struct MatrixSummary {
    rows: usize,
    synthetic_axis: usize,
    executable: usize,
    named_missing_primitive: usize,
    invalid_topology: usize,
    unsupported_primitive: usize,
    repair_required: usize,
    ambiguous_dependency: usize,
    downstream_risk: usize,
    unknown_family_rows: usize,
    unknown_input_fact: usize,
    unknown_fixture_reduction: usize,
    unknown_primitive_family: usize,
    unknown_tolerance_model: usize,
    unpartitioned_unknowns: usize,
}

impl MatrixSummary {
    fn from_rows(rows: &[apex_fillet_diagnostic_matrix_manifest::DiagnosticMatrixRow]) -> Self {
        let mut summary = Self {
            rows: rows.len(),
            ..Self::default()
        };
        for row in rows {
            match row.terminal {
                DiagnosticMatrixTerminal::ExecutableCandidate => summary.executable += 1,
                DiagnosticMatrixTerminal::NamedMissingPrimitive => {
                    summary.named_missing_primitive += 1;
                }
                DiagnosticMatrixTerminal::AmbiguousDependency => {
                    summary.ambiguous_dependency += 1;
                }
                DiagnosticMatrixTerminal::InvalidTopologyReason => summary.invalid_topology += 1,
                DiagnosticMatrixTerminal::ExplicitUnknownFamily => summary.unknown_family_rows += 1,
            }
            if row.fixture_tier == FixtureTier::SyntheticAxis {
                summary.synthetic_axis += 1;
            }
            if row.feasibility.contains("requires-")
                || row
                    .operation_dependency
                    .is_some_and(|dependency| dependency.contains("missing-"))
            {
                summary.unsupported_primitive += 1;
            }
            if row
                .geometry_feasibility
                .is_some_and(|readiness| readiness.contains("repair-required"))
                || row
                    .validity_repair
                    .is_some_and(|readiness| readiness.contains("requires-"))
            {
                summary.repair_required += 1;
            }
            if row
                .operation_dependency
                .is_some_and(|dependency| dependency.contains("ambiguous"))
            {
                summary.ambiguous_dependency += 1;
            }
            if row.downstream_readiness.is_some_and(|readiness| {
                readiness.contains("risk") || readiness.contains("blocked")
            }) {
                summary.downstream_risk += 1;
            }
            match row.unknown_facet {
                None => {}
                Some(UnknownDiagnosticFacet::UnknownInputFact) => summary.unknown_input_fact += 1,
                Some(UnknownDiagnosticFacet::UnknownFixtureReduction) => {
                    summary.unknown_fixture_reduction += 1;
                }
                Some(UnknownDiagnosticFacet::UnknownPrimitiveFamily) => {
                    summary.unknown_primitive_family += 1;
                }
                Some(UnknownDiagnosticFacet::UnknownToleranceModel) => {
                    summary.unknown_tolerance_model += 1;
                }
                Some(UnknownDiagnosticFacet::Unpartitioned) => {
                    summary.unpartitioned_unknowns += 1;
                }
            }
        }
        summary
    }

    fn diagnostic_label(self) -> String {
        format!(
            "apex_fillet_diagnostic_matrix=rows={} synthetic_axis={} executable={} named_missing_primitive={} invalid_topology={} unsupported_primitive={} repair_required={} ambiguous_dependency={} downstream_risk={} unknown_family={} unknowns=input_fact:{} fixture_reduction:{} primitive_family:{} tolerance_model:{} unpartitioned:{}",
            self.rows,
            self.synthetic_axis,
            self.executable,
            self.named_missing_primitive,
            self.invalid_topology,
            self.unsupported_primitive,
            self.repair_required,
            self.ambiguous_dependency,
            self.downstream_risk,
            self.unknown_family_rows,
            self.unknown_input_fact,
            self.unknown_fixture_reduction,
            self.unknown_primitive_family,
            self.unknown_tolerance_model,
            self.unpartitioned_unknowns
        )
    }
}
