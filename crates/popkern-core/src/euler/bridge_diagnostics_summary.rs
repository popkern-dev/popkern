use super::{
    ConstructionFamily, ConstructionFeasibility, ConstructionRequirement, GeometryFeasibility,
    OperatorReadinessStatus, ResourceConflictKind, TopologyDiagnosticGraph,
    ValidityRepairReadiness, ValidityRepairStatus,
};

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

impl BridgeDiagnosticSummary {
    pub(super) fn from_graph(graph: &TopologyDiagnosticGraph) -> Self {
        let mut summary = Self::empty();
        summary.operations = graph.operations.len();
        summary.resources = graph.resources.len();
        summary.conflicts = graph.conflicts.len();
        summary.readiness = graph.readiness.len();
        summary.construction_candidates = graph.candidates.len();
        summary.postconditions = graph.postconditions.len();
        for conflict in &graph.conflicts {
            summary.count_conflict(conflict.kind);
        }
        for readiness in &graph.readiness {
            summary.count_readiness(readiness.status);
        }
        for candidate in &graph.candidates {
            summary.count_candidate_family(candidate.family);
            summary.count_candidate_feasibility(candidate.feasibility);
            if candidate.postcondition_preview.is_some() {
                summary.candidates_with_postcondition_preview += 1;
            }
            for requirement in &candidate.requirements {
                summary.count_candidate_requirement(*requirement);
            }
            for fact in &candidate.geometry_facts {
                summary.geometry_feasibility_facts += 1;
                summary.count_geometry_feasibility(fact.feasibility);
            }
            for fact in &candidate.validity_repair_facts {
                summary.validity_repair_facts += 1;
                summary.count_validity_repair_readiness(fact.readiness);
                summary.count_validity_repair_status(fact.status);
            }
        }
        for postcondition in &graph.postconditions {
            summary.predicted_delta_vertices += postcondition.expected_delta.vertices;
            summary.predicted_delta_edges += postcondition.expected_delta.edges;
            summary.predicted_delta_faces += postcondition.expected_delta.faces;
            summary.predicted_delta_loops += postcondition.expected_delta.loops;
        }
        summary
    }

    /// Empty summary for paths that have no diagnostic request yet.
    pub fn empty() -> Self {
        Self {
            operations: 0,
            resources: 0,
            conflicts: 0,
            readiness: 0,
            construction_candidates: 0,
            postconditions: 0,
            unknown_candidate_families: 0,
            disjoint_cross_loop_bridge_candidates: 0,
            disjoint_same_loop_bridge_candidates: 0,
            paired_same_loop_twin_exposure_candidates: 0,
            boundary_side_separator_bridge_candidates: 0,
            unknown_candidate_feasibility: 0,
            executable_existing_operator_candidates: 0,
            executable_after_normalization_candidates: 0,
            requires_separator_subcomplex_candidates: 0,
            requires_custom_aggregate_authoring_candidates: 0,
            requires_cycle_regrouping_candidates: 0,
            invalid_non_manifold_candidates: 0,
            unsupported_unknown_family_candidates: 0,
            unknown_candidate_requirements: 0,
            live_topology_revalidation_requirements: 0,
            postcondition_verification_requirements: 0,
            candidates_with_postcondition_preview: 0,
            geometry_feasibility_facts: 0,
            geometry_ready_facts: 0,
            unsupported_geometry_family_facts: 0,
            geometry_repair_required_facts: 0,
            geometry_not_applicable_facts: 0,
            unknown_geometry_feasibility_facts: 0,
            validity_repair_facts: 0,
            valid_as_input_facts: 0,
            valid_after_non_mutating_normalization_facts: 0,
            requires_topology_repair_facts: 0,
            requires_geometry_repair_facts: 0,
            rejected_design_intent_changing_repair_facts: 0,
            validity_repair_not_applicable_facts: 0,
            unknown_validity_repair_facts: 0,
            no_repair_required_facts: 0,
            non_mutating_normalization_available_facts: 0,
            repair_required_facts: 0,
            repair_rejected_facts: 0,
            unknown_repair_status_facts: 0,
            predicted_delta_vertices: 0,
            predicted_delta_edges: 0,
            predicted_delta_faces: 0,
            predicted_delta_loops: 0,
            duplicate_half_edge_consumptions: 0,
            duplicate_boundary_edge_consumptions: 0,
            shared_start_overlaps: 0,
            shared_end_overlaps: 0,
            crossed_endpoint_overlaps: 0,
            aligned_overlaps: 0,
            nested_overlaps: 0,
            interior_overlaps: 0,
            repeated_edge_occurrences: 0,
            already_boundary_overlaps: 0,
            paired_side_needs_exposures: 0,
            unsupported_owner_conflicts: 0,
            stale_handle_conflicts: 0,
            unknown_conflicts: 0,
            ready_operators: 0,
            ready_after_normalization_operators: 0,
            missing_operators: 0,
            invalid_input_operators: 0,
            unknown_readiness: 0,
        }
    }

    /// Stable compact diagnostic label.
    pub fn diagnostic_label(&self) -> String {
        format!(
            "bridge_diagnostics=operations={} resources={} conflicts={} readiness={} construction_candidates={} postconditions={} conflicts_by_kind=duplicate_half_edge:{} duplicate_boundary_edge:{} shared_start:{} shared_end:{} crossed_endpoint:{} aligned:{} nested:{} interior:{} repeated_edge:{} already_boundary:{} paired_exposure:{} unsupported_owner:{} stale_handle:{} unknown:{} readiness_by_status=ready:{} ready_after_normalization:{} missing_operator:{} invalid_input:{} unknown:{} candidate_families=disjoint_cross_loop:{} disjoint_same_loop:{} paired_same_loop_twin_exposure:{} boundary_side_separator:{} unknown:{} candidate_feasibility=executable_existing:{} executable_after_normalization:{} requires_separator_subcomplex:{} requires_custom_aggregate_authoring:{} requires_cycle_regrouping:{} invalid_non_manifold:{} unsupported_unknown_family:{} unknown:{} candidate_requirements=live_revalidation:{} postcondition_verification:{} unknown:{} postcondition_preview_candidates={} geometry_feasibility=facts:{} ready:{} unsupported_family:{} repair_required:{} not_applicable:{} unknown:{} validity_repair=facts:{} valid_as_input:{} valid_after_non_mutating_normalization:{} requires_topology_repair:{} requires_geometry_repair:{} rejected_design_intent_changing:{} not_applicable:{} unknown:{} status_no_repair_required:{} status_non_mutating_normalization_available:{} status_repair_required:{} status_rejected:{} status_unknown:{} predicted_delta=dV={} dE={} dF={} dL={}",
            self.operations,
            self.resources,
            self.conflicts,
            self.readiness,
            self.construction_candidates,
            self.postconditions,
            self.duplicate_half_edge_consumptions,
            self.duplicate_boundary_edge_consumptions,
            self.shared_start_overlaps,
            self.shared_end_overlaps,
            self.crossed_endpoint_overlaps,
            self.aligned_overlaps,
            self.nested_overlaps,
            self.interior_overlaps,
            self.repeated_edge_occurrences,
            self.already_boundary_overlaps,
            self.paired_side_needs_exposures,
            self.unsupported_owner_conflicts,
            self.stale_handle_conflicts,
            self.unknown_conflicts,
            self.ready_operators,
            self.ready_after_normalization_operators,
            self.missing_operators,
            self.invalid_input_operators,
            self.unknown_readiness,
            self.disjoint_cross_loop_bridge_candidates,
            self.disjoint_same_loop_bridge_candidates,
            self.paired_same_loop_twin_exposure_candidates,
            self.boundary_side_separator_bridge_candidates,
            self.unknown_candidate_families,
            self.executable_existing_operator_candidates,
            self.executable_after_normalization_candidates,
            self.requires_separator_subcomplex_candidates,
            self.requires_custom_aggregate_authoring_candidates,
            self.requires_cycle_regrouping_candidates,
            self.invalid_non_manifold_candidates,
            self.unsupported_unknown_family_candidates,
            self.unknown_candidate_feasibility,
            self.live_topology_revalidation_requirements,
            self.postcondition_verification_requirements,
            self.unknown_candidate_requirements,
            self.candidates_with_postcondition_preview,
            self.geometry_feasibility_facts,
            self.geometry_ready_facts,
            self.unsupported_geometry_family_facts,
            self.geometry_repair_required_facts,
            self.geometry_not_applicable_facts,
            self.unknown_geometry_feasibility_facts,
            self.validity_repair_facts,
            self.valid_as_input_facts,
            self.valid_after_non_mutating_normalization_facts,
            self.requires_topology_repair_facts,
            self.requires_geometry_repair_facts,
            self.rejected_design_intent_changing_repair_facts,
            self.validity_repair_not_applicable_facts,
            self.unknown_validity_repair_facts,
            self.no_repair_required_facts,
            self.non_mutating_normalization_available_facts,
            self.repair_required_facts,
            self.repair_rejected_facts,
            self.unknown_repair_status_facts,
            self.predicted_delta_vertices,
            self.predicted_delta_edges,
            self.predicted_delta_faces,
            self.predicted_delta_loops
        )
    }

    fn count_conflict(&mut self, kind: ResourceConflictKind) {
        match kind {
            ResourceConflictKind::DuplicateHalfEdgeConsumption => {
                self.duplicate_half_edge_consumptions += 1;
            }
            ResourceConflictKind::DuplicateBoundaryEdgeConsumption => {
                self.duplicate_boundary_edge_consumptions += 1;
            }
            ResourceConflictKind::SharedStart => {
                self.shared_start_overlaps += 1;
            }
            ResourceConflictKind::SharedEnd => {
                self.shared_end_overlaps += 1;
            }
            ResourceConflictKind::CrossedEndpoint => {
                self.crossed_endpoint_overlaps += 1;
            }
            ResourceConflictKind::AlignedOverlap => {
                self.aligned_overlaps += 1;
            }
            ResourceConflictKind::NestedOverlap => {
                self.nested_overlaps += 1;
            }
            ResourceConflictKind::InteriorOverlap => {
                self.interior_overlaps += 1;
            }
            ResourceConflictKind::RepeatedEdgeOccurrence => {
                self.repeated_edge_occurrences += 1;
            }
            ResourceConflictKind::AlreadyBoundaryOverlap => {
                self.already_boundary_overlaps += 1;
            }
            ResourceConflictKind::PairedSideNeedsExposure => {
                self.paired_side_needs_exposures += 1;
            }
            ResourceConflictKind::UnsupportedOwner => {
                self.unsupported_owner_conflicts += 1;
            }
            ResourceConflictKind::StaleHandle => {
                self.stale_handle_conflicts += 1;
            }
            ResourceConflictKind::Unknown => {
                self.unknown_conflicts += 1;
            }
        }
    }

    fn count_readiness(&mut self, status: OperatorReadinessStatus) {
        match status {
            OperatorReadinessStatus::Ready => {
                self.ready_operators += 1;
            }
            OperatorReadinessStatus::ReadyAfterNormalization => {
                self.ready_after_normalization_operators += 1;
            }
            OperatorReadinessStatus::MissingOperator => {
                self.missing_operators += 1;
            }
            OperatorReadinessStatus::InvalidInput => {
                self.invalid_input_operators += 1;
            }
            OperatorReadinessStatus::Unknown => {
                self.unknown_readiness += 1;
            }
        }
    }

    fn count_candidate_family(&mut self, family: ConstructionFamily) {
        match family {
            ConstructionFamily::Unknown => self.unknown_candidate_families += 1,
            ConstructionFamily::DisjointCrossLoopBridge => {
                self.disjoint_cross_loop_bridge_candidates += 1;
            }
            ConstructionFamily::DisjointSameLoopBridge => {
                self.disjoint_same_loop_bridge_candidates += 1;
            }
            ConstructionFamily::PairedSameLoopTwinExposure => {
                self.paired_same_loop_twin_exposure_candidates += 1;
            }
            ConstructionFamily::BoundarySideSeparatorBridge => {
                self.boundary_side_separator_bridge_candidates += 1;
            }
        }
    }

    fn count_candidate_feasibility(&mut self, feasibility: ConstructionFeasibility) {
        match feasibility {
            ConstructionFeasibility::Unknown => self.unknown_candidate_feasibility += 1,
            ConstructionFeasibility::ExecutableByExistingOperator => {
                self.executable_existing_operator_candidates += 1;
            }
            ConstructionFeasibility::ExecutableAfterKnownNormalization => {
                self.executable_after_normalization_candidates += 1;
            }
            ConstructionFeasibility::RequiresSeparatorSubcomplex => {
                self.requires_separator_subcomplex_candidates += 1;
            }
            ConstructionFeasibility::RequiresCustomAggregateAuthoring => {
                self.requires_custom_aggregate_authoring_candidates += 1;
            }
            ConstructionFeasibility::RequiresCycleRegrouping => {
                self.requires_cycle_regrouping_candidates += 1;
            }
            ConstructionFeasibility::InvalidNonManifold => {
                self.invalid_non_manifold_candidates += 1;
            }
            ConstructionFeasibility::UnsupportedUnknownFamily => {
                self.unsupported_unknown_family_candidates += 1;
            }
        }
    }

    fn count_candidate_requirement(&mut self, requirement: ConstructionRequirement) {
        match requirement {
            ConstructionRequirement::Unknown => self.unknown_candidate_requirements += 1,
            ConstructionRequirement::LiveTopologyRevalidation => {
                self.live_topology_revalidation_requirements += 1;
            }
            ConstructionRequirement::PostconditionVerification => {
                self.postcondition_verification_requirements += 1;
            }
        }
    }

    fn count_geometry_feasibility(&mut self, feasibility: GeometryFeasibility) {
        match feasibility {
            GeometryFeasibility::Unknown => self.unknown_geometry_feasibility_facts += 1,
            GeometryFeasibility::NotApplicable => self.geometry_not_applicable_facts += 1,
            GeometryFeasibility::GeometryReady => self.geometry_ready_facts += 1,
            GeometryFeasibility::UnsupportedGeometryFamily => {
                self.unsupported_geometry_family_facts += 1;
            }
            GeometryFeasibility::RepairRequiredBeforeConstruction => {
                self.geometry_repair_required_facts += 1;
            }
        }
    }

    fn count_validity_repair_readiness(&mut self, readiness: ValidityRepairReadiness) {
        match readiness {
            ValidityRepairReadiness::Unknown => self.unknown_validity_repair_facts += 1,
            ValidityRepairReadiness::NotApplicable => {
                self.validity_repair_not_applicable_facts += 1
            }
            ValidityRepairReadiness::ValidAsInput => self.valid_as_input_facts += 1,
            ValidityRepairReadiness::ValidAfterNonMutatingNormalization => {
                self.valid_after_non_mutating_normalization_facts += 1;
            }
            ValidityRepairReadiness::RequiresTopologyRepair => {
                self.requires_topology_repair_facts += 1;
            }
            ValidityRepairReadiness::RequiresGeometryRepair => {
                self.requires_geometry_repair_facts += 1;
            }
            ValidityRepairReadiness::RejectedDesignIntentChangingRepair => {
                self.rejected_design_intent_changing_repair_facts += 1;
            }
        }
    }

    fn count_validity_repair_status(&mut self, status: ValidityRepairStatus) {
        match status {
            ValidityRepairStatus::Unknown => self.unknown_repair_status_facts += 1,
            ValidityRepairStatus::NoRepairRequired => self.no_repair_required_facts += 1,
            ValidityRepairStatus::NonMutatingNormalizationAvailable => {
                self.non_mutating_normalization_available_facts += 1;
            }
            ValidityRepairStatus::RepairRequired => self.repair_required_facts += 1,
            ValidityRepairStatus::Rejected => self.repair_rejected_facts += 1,
        }
    }
}
