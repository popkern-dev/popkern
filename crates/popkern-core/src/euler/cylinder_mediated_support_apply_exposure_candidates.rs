//! Typed construction candidates for cylinder-mediated exposure steps.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct ExposureConstructionCandidate {
    pub(super) family: ExposureConstructionFamily,
    pub(super) feasibility: ExposureConstructionFeasibility,
    pub(super) reason: &'static str,
    pub(super) spans: usize,
    pub(super) half_edges: usize,
    pub(super) expected_edges: i64,
    pub(super) expected_faces: i64,
    pub(super) expected_loops: i64,
    pub(super) semantic_ownership: ExposureSemanticOwnership,
    pub(super) rollback_boundary: ExposureRollbackBoundary,
}

impl ExposureConstructionCandidate {
    pub(super) fn mutation_gate_complete(self) -> bool {
        self.reason != "unknown"
            && (self.feasibility == ExposureConstructionFeasibility::NotApplicable
                || self.rollback_boundary == ExposureRollbackBoundary::CandidateBrepBeforeCommit)
    }

    pub(super) fn diagnostic_label(self) -> String {
        format!(
            "exposure_construction_candidate=family={} feasibility={} reason={} spans={} half_edges={} expected_delta=dV=0 dE={} dF={} dL={} semantic_ownership={} rollback_boundary={}",
            self.family.diagnostic_label(),
            self.feasibility.diagnostic_label(),
            self.reason,
            self.spans,
            self.half_edges,
            self.expected_edges,
            self.expected_faces,
            self.expected_loops,
            self.semantic_ownership.diagnostic_label(),
            self.rollback_boundary.diagnostic_label()
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ExposureConstructionFamily {
    OwnerBoundaryFragmentRetire,
    CylinderResourceBlockerRetire,
}

impl ExposureConstructionFamily {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::OwnerBoundaryFragmentRetire => "owner-boundary-fragment-retire",
            Self::CylinderResourceBlockerRetire => "cylinder-resource-blocker-retire",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ExposureConstructionFeasibility {
    NotApplicable,
    ExecutableByExistingOperator,
}

impl ExposureConstructionFeasibility {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::NotApplicable => "not-applicable",
            Self::ExecutableByExistingOperator => "executable-by-existing-operator",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ExposureSemanticOwnership {
    None,
    OutsideBoundarySideExposed,
    CylinderResourceBoundarySideExposed,
}

impl ExposureSemanticOwnership {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::OutsideBoundarySideExposed => "outside-boundary-side-exposed",
            Self::CylinderResourceBoundarySideExposed => "cylinder-resource-boundary-side-exposed",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ExposureRollbackBoundary {
    NoMutation,
    CandidateBrepBeforeCommit,
}

impl ExposureRollbackBoundary {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::NoMutation => "no-mutation",
            Self::CandidateBrepBeforeCommit => "candidate-brep-before-commit",
        }
    }
}
