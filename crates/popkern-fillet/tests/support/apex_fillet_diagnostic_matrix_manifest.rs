#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct DiagnosticMatrixRow {
    pub(crate) case_id: &'static str,
    pub(crate) fixture_tier: FixtureTier,
    pub(crate) covered_axes: &'static [DiagnosticAxis],
    pub(crate) terminal: DiagnosticMatrixTerminal,
    pub(crate) candidate_family: &'static str,
    pub(crate) feasibility: &'static str,
    pub(crate) geometry_feasibility_required: bool,
    pub(crate) geometry_feasibility: Option<&'static str>,
    pub(crate) validity_repair_required: bool,
    pub(crate) validity_repair: Option<&'static str>,
    pub(crate) operation_dependency_required: bool,
    pub(crate) operation_dependency: Option<&'static str>,
    pub(crate) downstream_readiness_required: bool,
    pub(crate) downstream_readiness: Option<&'static str>,
    pub(crate) expected_reason: &'static str,
    pub(crate) unknown_facet: Option<UnknownDiagnosticFacet>,
    pub(crate) closure_task: Option<&'static str>,
}

impl DiagnosticMatrixRow {
    pub(crate) fn diagnostic_label(self) -> String {
        format!(
            "case={} tier={} axes={} terminal={} candidate_family={} feasibility={} geometry_required={} geometry_feasibility={} validity_repair_required={} validity_repair={} operation_dependency_required={} operation_dependency={} downstream_required={} downstream_readiness={} reason={} unknown_facet={} closure_task={}",
            self.case_id,
            self.fixture_tier.diagnostic_label(),
            axis_list(self.covered_axes),
            self.terminal.diagnostic_label(),
            self.candidate_family,
            self.feasibility,
            self.geometry_feasibility_required,
            self.geometry_feasibility.unwrap_or("none"),
            self.validity_repair_required,
            self.validity_repair.unwrap_or("none"),
            self.operation_dependency_required,
            self.operation_dependency.unwrap_or("none"),
            self.downstream_readiness_required,
            self.downstream_readiness.unwrap_or("none"),
            self.expected_reason,
            self.unknown_facet
                .map(UnknownDiagnosticFacet::diagnostic_label)
                .unwrap_or("none"),
            self.closure_task.unwrap_or("none")
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FixtureTier {
    SeedRegression,
    SyntheticAxis,
    IgnoredZoo,
    CapturedCore,
}

impl FixtureTier {
    pub(crate) fn diagnostic_label(self) -> &'static str {
        match self {
            Self::SeedRegression => "seed-regression",
            Self::SyntheticAxis => "synthetic-axis",
            Self::IgnoredZoo => "ignored-zoo",
            Self::CapturedCore => "captured-core",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum DiagnosticAxis {
    CubeApex,
    PyramidApex,
    GregoryCap,
    GregoryRoundtrip,
    RepeatedApex,
    CylinderMediatedResource,
    BoundarySideSeparator,
    RetiredApexVertex,
    SeamCardinality,
    SupportOwnership,
    GeneratedBoundaryReplacement,
    TransitionContinuity,
    UnitScalePlacement,
    OrientationReversal,
    SplineRepresentation,
    OffsetSurfaceRepresentation,
    SelectionIntent,
    BranchStop,
    VariableLaw,
    OpenShellBoundary,
    HoleAdjacentLoop,
    BooleanSliver,
    ThinWallOffset,
    ExistingBlendAdjacency,
    SeamSingularity,
    SolverConditioning,
    FreeformApproximation,
    RepairPolicyBoundary,
    ZooExport,
    ExtremeRadiusRatio,
}

impl DiagnosticAxis {
    pub(crate) fn all() -> &'static [Self] {
        &[
            Self::CubeApex,
            Self::PyramidApex,
            Self::GregoryCap,
            Self::GregoryRoundtrip,
            Self::RepeatedApex,
            Self::CylinderMediatedResource,
            Self::BoundarySideSeparator,
            Self::RetiredApexVertex,
            Self::SeamCardinality,
            Self::SupportOwnership,
            Self::GeneratedBoundaryReplacement,
            Self::TransitionContinuity,
            Self::UnitScalePlacement,
            Self::OrientationReversal,
            Self::SplineRepresentation,
            Self::OffsetSurfaceRepresentation,
            Self::SelectionIntent,
            Self::BranchStop,
            Self::VariableLaw,
            Self::OpenShellBoundary,
            Self::HoleAdjacentLoop,
            Self::BooleanSliver,
            Self::ThinWallOffset,
            Self::ExistingBlendAdjacency,
            Self::SeamSingularity,
            Self::SolverConditioning,
            Self::FreeformApproximation,
            Self::RepairPolicyBoundary,
            Self::ZooExport,
            Self::ExtremeRadiusRatio,
        ]
    }

    pub(crate) fn diagnostic_label(self) -> &'static str {
        match self {
            Self::CubeApex => "cube-apex",
            Self::PyramidApex => "pyramid-apex",
            Self::GregoryCap => "gregory-cap",
            Self::GregoryRoundtrip => "gregory-roundtrip",
            Self::RepeatedApex => "repeated-apex",
            Self::CylinderMediatedResource => "cylinder-mediated-resource",
            Self::BoundarySideSeparator => "boundary-side-separator",
            Self::RetiredApexVertex => "retired-apex-vertex",
            Self::SeamCardinality => "seam-cardinality",
            Self::SupportOwnership => "support-ownership",
            Self::GeneratedBoundaryReplacement => "generated-boundary-replacement",
            Self::TransitionContinuity => "transition-continuity",
            Self::UnitScalePlacement => "unit-scale-placement",
            Self::OrientationReversal => "orientation-reversal",
            Self::SplineRepresentation => "spline-representation",
            Self::OffsetSurfaceRepresentation => "offset-surface-representation",
            Self::SelectionIntent => "selection-intent",
            Self::BranchStop => "branch-stop",
            Self::VariableLaw => "variable-law",
            Self::OpenShellBoundary => "open-shell-boundary",
            Self::HoleAdjacentLoop => "hole-adjacent-loop",
            Self::BooleanSliver => "boolean-sliver",
            Self::ThinWallOffset => "thin-wall-offset",
            Self::ExistingBlendAdjacency => "existing-blend-adjacency",
            Self::SeamSingularity => "seam-singularity",
            Self::SolverConditioning => "solver-conditioning",
            Self::FreeformApproximation => "freeform-approximation",
            Self::RepairPolicyBoundary => "repair-policy-boundary",
            Self::ZooExport => "zoo-export",
            Self::ExtremeRadiusRatio => "extreme-radius-ratio",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DiagnosticMatrixTerminal {
    ExecutableCandidate,
    NamedMissingPrimitive,
    AmbiguousDependency,
    InvalidTopologyReason,
    ExplicitUnknownFamily,
}

impl DiagnosticMatrixTerminal {
    pub(crate) fn all() -> &'static [Self] {
        &[
            Self::ExecutableCandidate,
            Self::NamedMissingPrimitive,
            Self::AmbiguousDependency,
            Self::InvalidTopologyReason,
            Self::ExplicitUnknownFamily,
        ]
    }

    pub(crate) fn diagnostic_label(self) -> &'static str {
        match self {
            Self::ExecutableCandidate => "executable-candidate",
            Self::NamedMissingPrimitive => "named-missing-primitive",
            Self::AmbiguousDependency => "ambiguous-dependency",
            Self::InvalidTopologyReason => "invalid-topology-reason",
            Self::ExplicitUnknownFamily => "explicit-unknown-family",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum UnknownDiagnosticFacet {
    UnknownInputFact,
    UnknownFixtureReduction,
    UnknownPrimitiveFamily,
    UnknownToleranceModel,
    Unpartitioned,
}

impl UnknownDiagnosticFacet {
    pub(crate) fn all() -> &'static [Self] {
        &[
            Self::UnknownInputFact,
            Self::UnknownFixtureReduction,
            Self::UnknownPrimitiveFamily,
            Self::UnknownToleranceModel,
            Self::Unpartitioned,
        ]
    }

    pub(crate) fn diagnostic_label(self) -> &'static str {
        match self {
            Self::UnknownInputFact => "unknown-input-fact",
            Self::UnknownFixtureReduction => "unknown-fixture-reduction",
            Self::UnknownPrimitiveFamily => "unknown-primitive-family",
            Self::UnknownToleranceModel => "unknown-tolerance-model",
            Self::Unpartitioned => "unpartitioned",
        }
    }
}

pub(crate) fn seed_matrix_rows() -> &'static [DiagnosticMatrixRow] {
    &apex_fillet_diagnostic_matrix_rows::SEED_MATRIX_ROWS
}

pub(crate) fn general_kernel_rows() -> &'static [DiagnosticMatrixRow] {
    &apex_fillet_diagnostic_matrix_general_rows::GENERAL_KERNEL_ROWS
}

fn axis_list(axes: &[DiagnosticAxis]) -> String {
    axes.iter()
        .map(|axis| axis.diagnostic_label())
        .collect::<Vec<_>>()
        .join("+")
}
#[path = "apex_fillet_diagnostic_matrix_general_rows.rs"]
mod apex_fillet_diagnostic_matrix_general_rows;
#[path = "apex_fillet_diagnostic_matrix_rows.rs"]
mod apex_fillet_diagnostic_matrix_rows;
