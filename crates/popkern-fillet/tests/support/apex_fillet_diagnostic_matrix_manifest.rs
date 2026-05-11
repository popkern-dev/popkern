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
    IgnoredZoo,
    CapturedCore,
}

impl FixtureTier {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::SeedRegression => "seed-regression",
            Self::IgnoredZoo => "ignored-zoo",
            Self::CapturedCore => "captured-core",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DiagnosticAxis {
    CubeApex,
    PyramidApex,
    GregoryCap,
    GregoryRoundtrip,
    RepeatedApex,
    CylinderMediatedResource,
    BoundarySideSeparator,
    ZooExport,
    ExtremeRadiusRatio,
}

impl DiagnosticAxis {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::CubeApex => "cube-apex",
            Self::PyramidApex => "pyramid-apex",
            Self::GregoryCap => "gregory-cap",
            Self::GregoryRoundtrip => "gregory-roundtrip",
            Self::RepeatedApex => "repeated-apex",
            Self::CylinderMediatedResource => "cylinder-mediated-resource",
            Self::BoundarySideSeparator => "boundary-side-separator",
            Self::ZooExport => "zoo-export",
            Self::ExtremeRadiusRatio => "extreme-radius-ratio",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DiagnosticMatrixTerminal {
    ExecutableCandidate,
    NamedMissingPrimitive,
    InvalidTopologyReason,
    ExplicitUnknownFamily,
}

impl DiagnosticMatrixTerminal {
    pub(crate) fn all() -> &'static [Self] {
        &[
            Self::ExecutableCandidate,
            Self::NamedMissingPrimitive,
            Self::InvalidTopologyReason,
            Self::ExplicitUnknownFamily,
        ]
    }

    pub(crate) fn diagnostic_label(self) -> &'static str {
        match self {
            Self::ExecutableCandidate => "executable-candidate",
            Self::NamedMissingPrimitive => "named-missing-primitive",
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
    &SEED_MATRIX_ROWS
}

const SEED_MATRIX_ROWS: [DiagnosticMatrixRow; 13] = [
    DiagnosticMatrixRow {
        case_id: "cube_mixed_native_a2",
        fixture_tier: FixtureTier::SeedRegression,
        covered_axes: &[
            DiagnosticAxis::CubeApex,
            DiagnosticAxis::GregoryCap,
            DiagnosticAxis::RepeatedApex,
        ],
        terminal: DiagnosticMatrixTerminal::ExecutableCandidate,
        candidate_family: "native-a2-gregory-repeated-apex",
        feasibility: "executable",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("geometry-ready"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("not-applicable"),
        downstream_readiness_required: true,
        downstream_readiness: Some("export-mesh-ready"),
        expected_reason: "native-a2-completes-without-fallback",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "cube_extreme_native_a2",
        fixture_tier: FixtureTier::SeedRegression,
        covered_axes: &[
            DiagnosticAxis::CubeApex,
            DiagnosticAxis::GregoryCap,
            DiagnosticAxis::RepeatedApex,
            DiagnosticAxis::ExtremeRadiusRatio,
        ],
        terminal: DiagnosticMatrixTerminal::ExecutableCandidate,
        candidate_family: "native-a2-gregory-repeated-apex",
        feasibility: "executable",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("geometry-ready"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("not-applicable"),
        downstream_readiness_required: true,
        downstream_readiness: Some("export-mesh-ready"),
        expected_reason: "native-a2-completes-without-fallback",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "pyramid3_uniform_native_a2",
        fixture_tier: FixtureTier::SeedRegression,
        covered_axes: &[DiagnosticAxis::PyramidApex, DiagnosticAxis::GregoryCap],
        terminal: DiagnosticMatrixTerminal::ExecutableCandidate,
        candidate_family: "native-a2-gregory-anchor-incompatible-pyramid",
        feasibility: "executable",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("geometry-ready"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("patch-filling-executable"),
        downstream_readiness_required: true,
        downstream_readiness: Some("continuity-risk-managed"),
        expected_reason: "native-gregory-used-without-fallback",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "pyramid4_mixed_native_a2",
        fixture_tier: FixtureTier::SeedRegression,
        covered_axes: &[
            DiagnosticAxis::PyramidApex,
            DiagnosticAxis::GregoryCap,
            DiagnosticAxis::RepeatedApex,
        ],
        terminal: DiagnosticMatrixTerminal::ExecutableCandidate,
        candidate_family: "native-a2-partitioned-repeated-apex",
        feasibility: "executable",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("geometry-ready"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("patch-filling-executable"),
        downstream_readiness_required: true,
        downstream_readiness: Some("continuity-risk-managed"),
        expected_reason: "native-a2-completes-without-fallback",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "pyramid4_uniform_native_a2",
        fixture_tier: FixtureTier::SeedRegression,
        covered_axes: &[DiagnosticAxis::PyramidApex, DiagnosticAxis::GregoryCap],
        terminal: DiagnosticMatrixTerminal::ExecutableCandidate,
        candidate_family: "native-a2-partitioned-repeated-apex",
        feasibility: "executable",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("geometry-ready"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("patch-filling-executable"),
        downstream_readiness_required: true,
        downstream_readiness: Some("continuity-risk-managed"),
        expected_reason: "native-a2-completes-without-fallback",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "gregory_step_roundtrip_tolerance",
        fixture_tier: FixtureTier::SeedRegression,
        covered_axes: &[DiagnosticAxis::GregoryCap, DiagnosticAxis::GregoryRoundtrip],
        terminal: DiagnosticMatrixTerminal::ExecutableCandidate,
        candidate_family: "gregory-step-roundtrip",
        feasibility: "executable-with-local-tolerance",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("geometry-ready"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("export-roundtrip-tolerance-check"),
        downstream_readiness_required: true,
        downstream_readiness: Some("export-mesh-ready"),
        expected_reason: "gregory-step-roundtrip-within-local-tolerance",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "cube_extreme_long_interval_resource_refusal",
        fixture_tier: FixtureTier::CapturedCore,
        covered_axes: &[
            DiagnosticAxis::CylinderMediatedResource,
            DiagnosticAxis::ExtremeRadiusRatio,
        ],
        terminal: DiagnosticMatrixTerminal::NamedMissingPrimitive,
        candidate_family: "split-cylinder-resource-merge",
        feasibility: "requires-resource-merge",
        geometry_feasibility_required: false,
        geometry_feasibility: Some("not-applicable"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("missing-resource-merge-primitive"),
        downstream_readiness_required: true,
        downstream_readiness: Some("blocked-missing-primitive"),
        expected_reason: "split-cylinder-resource-merge-required",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "cube_extreme_boundary_side_separator",
        fixture_tier: FixtureTier::CapturedCore,
        covered_axes: &[
            DiagnosticAxis::CylinderMediatedResource,
            DiagnosticAxis::BoundarySideSeparator,
        ],
        terminal: DiagnosticMatrixTerminal::NamedMissingPrimitive,
        candidate_family: "boundary-side-separator-bridge",
        feasibility: "requires-separator-subcomplex",
        geometry_feasibility_required: false,
        geometry_feasibility: Some("not-applicable"),
        validity_repair_required: true,
        validity_repair: Some("valid-as-input"),
        operation_dependency_required: true,
        operation_dependency: Some("missing-boundary-side-separator-subcomplex"),
        downstream_readiness_required: true,
        downstream_readiness: Some("blocked-missing-primitive"),
        expected_reason: "missing-boundary-side-separator-subcomplex",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "pyramid4_zoo_mixed_export_defects",
        fixture_tier: FixtureTier::IgnoredZoo,
        covered_axes: &[DiagnosticAxis::PyramidApex, DiagnosticAxis::ZooExport],
        terminal: DiagnosticMatrixTerminal::InvalidTopologyReason,
        candidate_family: "native-a2-zoo-export-validation",
        feasibility: "invalid-topology",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("repair-required-before-construction"),
        validity_repair_required: true,
        validity_repair: Some("requires-topology-repair"),
        operation_dependency_required: true,
        operation_dependency: Some("blocked-invalid-input"),
        downstream_readiness_required: true,
        downstream_readiness: Some("export-risk-invalid-topology"),
        expected_reason: "retired-apex-support-topology-export-defects",
        unknown_facet: None,
        closure_task: None,
    },
    DiagnosticMatrixRow {
        case_id: "cube_extreme_40_upstream_gap",
        fixture_tier: FixtureTier::IgnoredZoo,
        covered_axes: &[DiagnosticAxis::CubeApex, DiagnosticAxis::ExtremeRadiusRatio],
        terminal: DiagnosticMatrixTerminal::NamedMissingPrimitive,
        candidate_family: "upstream-boundary-injection",
        feasibility: "requires-topology-diagnostic-row",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("repair-required-before-construction"),
        validity_repair_required: true,
        validity_repair: Some("requires-geometry-repair"),
        operation_dependency_required: true,
        operation_dependency: Some("missing-upstream-boundary-injection"),
        downstream_readiness_required: true,
        downstream_readiness: Some("blocked-repair-required"),
        expected_reason: "upstream-gap-before-gregory-ladder",
        unknown_facet: Some(UnknownDiagnosticFacet::UnknownInputFact),
        closure_task: Some("PK-005"),
    },
    DiagnosticMatrixRow {
        case_id: "cube_extreme_50_upstream_gap",
        fixture_tier: FixtureTier::IgnoredZoo,
        covered_axes: &[DiagnosticAxis::CubeApex, DiagnosticAxis::ExtremeRadiusRatio],
        terminal: DiagnosticMatrixTerminal::NamedMissingPrimitive,
        candidate_family: "upstream-boundary-injection",
        feasibility: "requires-topology-diagnostic-row",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("repair-required-before-construction"),
        validity_repair_required: true,
        validity_repair: Some("requires-geometry-repair"),
        operation_dependency_required: true,
        operation_dependency: Some("missing-upstream-boundary-injection"),
        downstream_readiness_required: true,
        downstream_readiness: Some("blocked-repair-required"),
        expected_reason: "upstream-gap-before-gregory-ladder",
        unknown_facet: Some(UnknownDiagnosticFacet::UnknownInputFact),
        closure_task: Some("PK-005"),
    },
    DiagnosticMatrixRow {
        case_id: "cube_extreme_100_upstream_gap",
        fixture_tier: FixtureTier::IgnoredZoo,
        covered_axes: &[DiagnosticAxis::CubeApex, DiagnosticAxis::ExtremeRadiusRatio],
        terminal: DiagnosticMatrixTerminal::NamedMissingPrimitive,
        candidate_family: "upstream-boundary-injection",
        feasibility: "requires-topology-diagnostic-row",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("repair-required-before-construction"),
        validity_repair_required: true,
        validity_repair: Some("requires-geometry-repair"),
        operation_dependency_required: true,
        operation_dependency: Some("missing-upstream-boundary-injection"),
        downstream_readiness_required: true,
        downstream_readiness: Some("blocked-repair-required"),
        expected_reason: "upstream-gap-before-gregory-ladder",
        unknown_facet: Some(UnknownDiagnosticFacet::UnknownInputFact),
        closure_task: Some("PK-005"),
    },
    DiagnosticMatrixRow {
        case_id: "pentagon5_mixed_polar_gap",
        fixture_tier: FixtureTier::IgnoredZoo,
        covered_axes: &[DiagnosticAxis::PyramidApex, DiagnosticAxis::GregoryCap],
        terminal: DiagnosticMatrixTerminal::NamedMissingPrimitive,
        candidate_family: "high-valence-polar-transition",
        feasibility: "requires-geometry-feasibility",
        geometry_feasibility_required: true,
        geometry_feasibility: Some("unsupported-geometry-family"),
        validity_repair_required: true,
        validity_repair: Some("valid-after-non-mutating-normalization"),
        operation_dependency_required: true,
        operation_dependency: Some("missing-high-valence-polar-transition"),
        downstream_readiness_required: true,
        downstream_readiness: Some("blocked-unsupported-geometry"),
        expected_reason: "polar-gap-before-gregory-ladder",
        unknown_facet: Some(UnknownDiagnosticFacet::UnknownInputFact),
        closure_task: Some("PK-005"),
    },
];

fn axis_list(axes: &[DiagnosticAxis]) -> String {
    axes.iter()
        .map(|axis| axis.diagnostic_label())
        .collect::<Vec<_>>()
        .join("+")
}
