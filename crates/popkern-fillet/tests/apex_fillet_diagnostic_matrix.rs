//! Apex/fillet diagnostic classification matrix.
//!
//! This test is intentionally lightweight: it locks the structured manifest
//! and counters that live fixture rows must satisfy as diagnostics mature.

#[path = "support/apex_fillet_diagnostic_matrix_manifest.rs"]
mod apex_fillet_diagnostic_matrix_manifest;

use apex_fillet_diagnostic_matrix_manifest::{
    DiagnosticMatrixTerminal, UnknownDiagnosticFacet, seed_matrix_rows,
};

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
        assert!(
            label.contains("terminal="),
            "{} lacks terminal classification: {label}",
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
        "cube_mixed_native_a2",
        "cube_extreme_native_a2",
        "pyramid3_uniform_native_a2",
        "pyramid4_uniform_native_a2",
        "pyramid4_mixed_native_a2",
        "gregory_step_roundtrip_tolerance",
        "cube_extreme_40_upstream_gap",
        "cube_extreme_50_upstream_gap",
        "cube_extreme_100_upstream_gap",
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
                DiagnosticMatrixTerminal::InvalidTopologyReason => summary.invalid_topology += 1,
                DiagnosticMatrixTerminal::ExplicitUnknownFamily => summary.unknown_family_rows += 1,
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
            "apex_fillet_diagnostic_matrix=rows={} executable={} named_missing_primitive={} invalid_topology={} unsupported_primitive={} repair_required={} ambiguous_dependency={} downstream_risk={} unknown_family={} unknowns=input_fact:{} fixture_reduction:{} primitive_family:{} tolerance_model:{} unpartitioned:{}",
            self.rows,
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
