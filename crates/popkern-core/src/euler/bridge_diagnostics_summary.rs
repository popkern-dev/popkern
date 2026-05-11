use super::{
    BridgeDiagnosticSummary, ConstructionFamily, ConstructionFeasibility, ConstructionRequirement,
    GeometryFeasibility, OperationDependencyReadiness, OperationDependencyStatus,
    OperatorReadinessStatus, ResourceConflictKind, TopologyDiagnosticGraph,
    ValidityRepairReadiness, ValidityRepairStatus,
};

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
            for dependency in &candidate.operation_dependencies {
                summary.operation_dependency_facts += 1;
                summary.count_operation_dependency_readiness(dependency.readiness);
                summary.count_operation_dependency_status(dependency.status);
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
            operation_dependency_facts: 0,
            executable_operation_dependency_facts: 0,
            ambiguous_operation_dependency_facts: 0,
            invalid_operation_dependency_facts: 0,
            missing_primitive_operation_dependency_facts: 0,
            operation_dependency_not_applicable_facts: 0,
            unknown_operation_dependency_facts: 0,
            operation_dependency_satisfied_facts: 0,
            operation_dependency_not_required_facts: 0,
            operation_dependency_required_missing_primitive_facts: 0,
            operation_dependency_blocked_invalid_input_facts: 0,
            operation_dependency_ambiguous_needs_resolver_facts: 0,
            unknown_operation_dependency_status_facts: 0,
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

    fn count_operation_dependency_readiness(&mut self, readiness: OperationDependencyReadiness) {
        match readiness {
            OperationDependencyReadiness::Unknown => self.unknown_operation_dependency_facts += 1,
            OperationDependencyReadiness::NotApplicable => {
                self.operation_dependency_not_applicable_facts += 1;
            }
            OperationDependencyReadiness::Executable => {
                self.executable_operation_dependency_facts += 1;
            }
            OperationDependencyReadiness::Ambiguous => {
                self.ambiguous_operation_dependency_facts += 1;
            }
            OperationDependencyReadiness::Invalid => {
                self.invalid_operation_dependency_facts += 1;
            }
            OperationDependencyReadiness::MissingPrimitive => {
                self.missing_primitive_operation_dependency_facts += 1;
            }
        }
    }

    fn count_operation_dependency_status(&mut self, status: OperationDependencyStatus) {
        match status {
            OperationDependencyStatus::Unknown => {
                self.unknown_operation_dependency_status_facts += 1;
            }
            OperationDependencyStatus::Satisfied => {
                self.operation_dependency_satisfied_facts += 1;
            }
            OperationDependencyStatus::NotRequired => {
                self.operation_dependency_not_required_facts += 1;
            }
            OperationDependencyStatus::RequiredMissingPrimitive => {
                self.operation_dependency_required_missing_primitive_facts += 1;
            }
            OperationDependencyStatus::BlockedInvalidInput => {
                self.operation_dependency_blocked_invalid_input_facts += 1;
            }
            OperationDependencyStatus::AmbiguousNeedsResolver => {
                self.operation_dependency_ambiguous_needs_resolver_facts += 1;
            }
        }
    }
}
