use super::{OperatorReadinessStatus, ResourceConflictKind, TopologyDiagnosticGraph};

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
    /// Number of postcondition predictions in the graph.
    pub postconditions: usize,
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
        summary.postconditions = graph.postconditions.len();
        for conflict in &graph.conflicts {
            summary.count_conflict(conflict.kind);
        }
        for readiness in &graph.readiness {
            summary.count_readiness(readiness.status);
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
            postconditions: 0,
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
            "bridge_diagnostics=operations={} resources={} conflicts={} readiness={} postconditions={} conflicts_by_kind=duplicate_half_edge:{} duplicate_boundary_edge:{} shared_start:{} shared_end:{} crossed_endpoint:{} aligned:{} nested:{} interior:{} repeated_edge:{} already_boundary:{} paired_exposure:{} unsupported_owner:{} stale_handle:{} unknown:{} readiness_by_status=ready:{} ready_after_normalization:{} missing_operator:{} invalid_input:{} unknown:{}",
            self.operations,
            self.resources,
            self.conflicts,
            self.readiness,
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
            self.unknown_readiness
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
}
