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

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
struct MatrixSummary {
    rows: usize,
    executable: usize,
    named_missing_primitive: usize,
    invalid_topology: usize,
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
            "apex_fillet_diagnostic_matrix=rows={} executable={} named_missing_primitive={} invalid_topology={} unknown_family={} unknowns=input_fact:{} fixture_reduction:{} primitive_family:{} tolerance_model:{} unpartitioned:{}",
            self.rows,
            self.executable,
            self.named_missing_primitive,
            self.invalid_topology,
            self.unknown_family_rows,
            self.unknown_input_fact,
            self.unknown_fixture_reduction,
            self.unknown_primitive_family,
            self.unknown_tolerance_model,
            self.unpartitioned_unknowns
        )
    }
}
