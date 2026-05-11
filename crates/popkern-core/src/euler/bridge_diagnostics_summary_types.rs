//! Compact bridge diagnostic summary counters.

/// Compact bridge diagnostic summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeDiagnosticSummary {
    /// Number of operations in the graph.
    pub operations: usize,
    /// Number of resources in the graph.
    pub resources: usize,
    /// Number of conflict components in the graph.
    pub conflicts: usize,
    /// Number of readiness entries in the graph.
    pub readiness: usize,
    /// Number of construction candidates in the graph.
    pub construction_candidates: usize,
    /// Number of postcondition predictions in the graph.
    pub postconditions: usize,
    /// Candidates with unknown family.
    pub unknown_candidate_families: usize,
    /// Disjoint cross-loop bridge candidates.
    pub disjoint_cross_loop_bridge_candidates: usize,
    /// Disjoint same-loop bridge candidates.
    pub disjoint_same_loop_bridge_candidates: usize,
    /// Paired same-loop twin exposure candidates.
    pub paired_same_loop_twin_exposure_candidates: usize,
    /// Boundary-side separator bridge candidates.
    pub boundary_side_separator_bridge_candidates: usize,
    /// Candidates with unknown feasibility.
    pub unknown_candidate_feasibility: usize,
    /// Candidates executable by existing operators.
    pub executable_existing_operator_candidates: usize,
    /// Candidates executable after known normalization.
    pub executable_after_normalization_candidates: usize,
    /// Candidates requiring separator subcomplex support.
    pub requires_separator_subcomplex_candidates: usize,
    /// Candidates requiring aggregate authoring support.
    pub requires_custom_aggregate_authoring_candidates: usize,
    /// Candidates requiring cycle regrouping.
    pub requires_cycle_regrouping_candidates: usize,
    /// Candidates rejected as non-manifold.
    pub invalid_non_manifold_candidates: usize,
    /// Unsupported candidate families.
    pub unsupported_unknown_family_candidates: usize,
    /// Unknown construction requirements.
    pub unknown_candidate_requirements: usize,
    /// Live-topology revalidation requirements.
    pub live_topology_revalidation_requirements: usize,
    /// Postcondition verification requirements.
    pub postcondition_verification_requirements: usize,
    /// Candidates carrying a postcondition preview.
    pub candidates_with_postcondition_preview: usize,
    /// Geometry feasibility facts attached to candidates.
    pub geometry_feasibility_facts: usize,
    /// Geometry-ready facts.
    pub geometry_ready_facts: usize,
    /// Unsupported geometry-family facts.
    pub unsupported_geometry_family_facts: usize,
    /// Repair-required geometry facts.
    pub geometry_repair_required_facts: usize,
    /// Not-applicable geometry facts.
    pub geometry_not_applicable_facts: usize,
    /// Unknown geometry facts.
    pub unknown_geometry_feasibility_facts: usize,
    /// Validity/repair facts attached to candidates.
    pub validity_repair_facts: usize,
    /// Candidate inputs valid as-is.
    pub valid_as_input_facts: usize,
    /// Candidate inputs valid after non-mutating normalization.
    pub valid_after_non_mutating_normalization_facts: usize,
    /// Candidate inputs requiring topology repair.
    pub requires_topology_repair_facts: usize,
    /// Candidate inputs requiring geometry repair.
    pub requires_geometry_repair_facts: usize,
    /// Candidate repairs rejected because they would change design intent.
    pub rejected_design_intent_changing_repair_facts: usize,
    /// Not-applicable validity/repair facts.
    pub validity_repair_not_applicable_facts: usize,
    /// Unknown validity/repair facts.
    pub unknown_validity_repair_facts: usize,
    /// Validity facts that need no repair.
    pub no_repair_required_facts: usize,
    /// Validity facts with non-mutating normalization available.
    pub non_mutating_normalization_available_facts: usize,
    /// Validity facts requiring repair.
    pub repair_required_facts: usize,
    /// Validity facts rejected by repair policy.
    pub repair_rejected_facts: usize,
    /// Validity facts with unknown repair status.
    pub unknown_repair_status_facts: usize,
    /// Operation dependency facts attached to candidates.
    pub operation_dependency_facts: usize,
    /// Executable operation dependency facts.
    pub executable_operation_dependency_facts: usize,
    /// Ambiguous operation dependency facts.
    pub ambiguous_operation_dependency_facts: usize,
    /// Invalid operation dependency facts.
    pub invalid_operation_dependency_facts: usize,
    /// Missing-primitive operation dependency facts.
    pub missing_primitive_operation_dependency_facts: usize,
    /// Not-applicable operation dependency facts.
    pub operation_dependency_not_applicable_facts: usize,
    /// Unknown operation dependency facts.
    pub unknown_operation_dependency_facts: usize,
    /// Operation dependency facts satisfied by existing facts/operators.
    pub operation_dependency_satisfied_facts: usize,
    /// Operation dependency facts not required.
    pub operation_dependency_not_required_facts: usize,
    /// Operation dependencies requiring a missing primitive.
    pub operation_dependency_required_missing_primitive_facts: usize,
    /// Operation dependency facts blocked by invalid input.
    pub operation_dependency_blocked_invalid_input_facts: usize,
    /// Operation dependency facts needing ambiguity resolution.
    pub operation_dependency_ambiguous_needs_resolver_facts: usize,
    /// Operation dependency facts with unknown status.
    pub unknown_operation_dependency_status_facts: usize,
    /// Sum of predicted vertex-count deltas.
    pub predicted_delta_vertices: i64,
    /// Sum of predicted edge-count deltas.
    pub predicted_delta_edges: i64,
    /// Sum of predicted face-count deltas.
    pub predicted_delta_faces: i64,
    /// Sum of predicted loop-count deltas.
    pub predicted_delta_loops: i64,
    /// Duplicate half-edge side consumption conflicts.
    pub duplicate_half_edge_consumptions: usize,
    /// Duplicate boundary-edge consumption conflicts.
    pub duplicate_boundary_edge_consumptions: usize,
    /// Shared-start overlap conflicts.
    pub shared_start_overlaps: usize,
    /// Shared-end overlap conflicts.
    pub shared_end_overlaps: usize,
    /// Crossed-endpoint overlap conflicts.
    pub crossed_endpoint_overlaps: usize,
    /// Aligned interval overlap conflicts.
    pub aligned_overlaps: usize,
    /// Nested interval overlap conflicts.
    pub nested_overlaps: usize,
    /// Interior interval overlap conflicts.
    pub interior_overlaps: usize,
    /// Repeated edge occurrence conflicts.
    pub repeated_edge_occurrences: usize,
    /// Already-boundary overlap conflicts.
    pub already_boundary_overlaps: usize,
    /// Paired-side exposure requirements.
    pub paired_side_needs_exposures: usize,
    /// Unsupported-owner conflicts.
    pub unsupported_owner_conflicts: usize,
    /// Stale-handle conflicts.
    pub stale_handle_conflicts: usize,
    /// Unknown conflict kinds.
    pub unknown_conflicts: usize,
    /// Ready operator entries.
    pub ready_operators: usize,
    /// Ready-after-normalization operator entries.
    pub ready_after_normalization_operators: usize,
    /// Missing operator entries.
    pub missing_operators: usize,
    /// Invalid input operator entries.
    pub invalid_input_operators: usize,
    /// Unknown readiness entries.
    pub unknown_readiness: usize,
}
