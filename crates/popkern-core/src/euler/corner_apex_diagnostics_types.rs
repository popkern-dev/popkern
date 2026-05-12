//! Corner/apex diagnostic data types.
#![allow(missing_docs)]

use super::super::{DiagnosticHandle, DiagnosticTolerance};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CornerApexDiagnosticRequest {
    pub scope: LocalCornerScope,
    pub selected_edges: Vec<LocalBlendEdgeSelection>,
    pub selection_intent: &'static str,
    pub operation_context: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocalCornerScope {
    pub vertex: DiagnosticHandle,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocalBlendEdgeSelection {
    pub edge: DiagnosticHandle,
    pub blend_kind: BlendKind,
    pub radius_nanometers: Option<i64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CornerApexDiagnosticReport {
    pub request: CornerApexDiagnosticRequest,
    pub topology: CornerTopologyRecord,
    pub classification: VertexClassificationFact,
    pub sizing: ApexSizingFact,
    pub trial_interaction: TrialInteractionFact,
    pub support_adequacy: CornerSupportAdequacyFact,
    pub strategy: StrategyReadinessFact,
    pub construction_risk: ConstructionRiskFact,
    pub integration_validation: IntegrationValidationFact,
    pub lifecycle: LifecycleFact,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CornerTopologyRecord {
    pub vertex: DiagnosticHandle,
    pub incident_blend_edges: usize,
    pub incident_non_blend_edges: usize,
    pub surrounding_faces: usize,
    pub face_adjacency_edges: usize,
    pub state: CornerTopologyState,
    pub geometry: CornerGeometryFact,
    pub surface_representation: SurfaceRepresentation,
    pub surface_compatibility: SurfaceCompatibility,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerTopologyState {
    Ready,
    StaleHandle,
    MissingHandle,
    NonManifold,
    Boundary,
    OpenShell,
    WireOrSheet,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CornerGeometryFact {
    pub tangent: CornerGeometryAvailability,
    pub normal: CornerGeometryAvailability,
    pub dihedral: CornerGeometryAvailability,
    pub tolerance: DiagnosticTolerance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerGeometryAvailability {
    Available,
    NotApplicable,
    MissingInputFact,
    BlockedByDependency(GeometryProbeDependency),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GeometryProbeDependency {
    CurveTangent,
    SurfaceNormal,
    DihedralAngle,
    Projection,
    Intersection,
    CenterLocus,
    SupportFootprint,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SurfaceRepresentation {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Spline,
    SweptOrRevolved,
    Offset,
    ExistingBlend,
    UnknownImported,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SurfaceCompatibility {
    Full,
    Conditional(&'static str),
    Incompatible(&'static str),
    Unknown(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BlendKind {
    FilletConstantRadius,
    FilletVariableRadius,
    Chamfer,
    G2Blend,
    Setback,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BlendValence {
    NotBlendVertex,
    FreeEnd,
    ThroughVertex,
    StandardCorner,
    HighValenceCorner,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EdgeConvexity {
    Convex,
    Concave,
    NearTangent,
    VariableNearVertex,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VertexConvexityClass {
    AllConvex,
    AllConcave,
    Mixed,
    NearTangent,
    VariableConvexityStop,
    UnknownInputFact,
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BlendTypeCompatibility {
    Homogeneous,
    ParameterVarying,
    Heterogeneous,
    Incompatible(&'static str),
    UnknownInputFact,
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TopologicalConfiguration {
    SimpleCorner,
    NonManifold,
    BlendOnBlend,
    NestedVertexBlend,
    BoundaryVertex,
    OpenShellOrSheet,
    WireLike,
    UnknownInputFact,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VertexClassificationFact {
    pub valence: BlendValence,
    pub convexity: VertexConvexityClass,
    pub blend_types: BlendTypeCompatibility,
    pub topology: TopologicalConfiguration,
    pub exposure: CornerSupportAdequacyState,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ApexSizingMode {
    ConstantRadius,
    VariableRadiusEndpoint,
    RangeRange,
    RangeAngle,
    AngleRange,
    NotApplicable,
    UnknownInputFact,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ApexSizingState {
    Ready,
    RadiusTooLarge,
    AdjacentApexCollision,
    ZeroOrNegativeRange,
    HighCurvatureRisk,
    MissingSolverFact,
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ApexSizingInput {
    pub mode: ApexSizingMode,
    pub requested_extent_nm: Option<i64>,
    pub available_extent_nm: Option<i64>,
    pub adjacent_clearance_nm: Option<i64>,
    pub curvature_radius_nm: Option<i64>,
    pub tolerance: DiagnosticTolerance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ApexSizingFact {
    pub mode: ApexSizingMode,
    pub state: ApexSizingState,
    pub max_feasible_radius_nm: Option<i64>,
    pub margin_nm: Option<i64>,
    pub tolerance: DiagnosticTolerance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TrialInteractionState {
    WellConditioned,
    NearTangentUpgrade,
    SelfIntersecting,
    Ambiguous,
    MissingIntersectionPrimitive,
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TrialInteractionFact {
    pub state: TrialInteractionState,
    pub dependency: Option<GeometryProbeDependency>,
    pub tolerance: DiagnosticTolerance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerSupportAdequacyState {
    Fits,
    Overflow,
    Cliff,
    HoleOrInnerLoopCollision,
    TrimDomainExit,
    AdjacentBlendCollision,
    InsideTight,
    UnknownInputFact,
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CornerSupportAdequacyFact {
    pub state: CornerSupportAdequacyState,
    pub dependency: Option<GeometryProbeDependency>,
    pub tolerance: DiagnosticTolerance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerApexStrategyFamily {
    CapRunout,
    CliffTermination,
    MutualTrim,
    ThreeFaceBlend,
    VertexPatch,
    HighValencePatch,
    Setback,
    OverflowSmoothing,
    DeferredVertexPatch,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerApexStrategyReadiness {
    Ready,
    ReadyWithManagedRisk,
    RequiresMissingPrimitive(&'static str),
    Invalid(&'static str),
    UnknownInputFact,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StrategyReadinessFact {
    pub family: CornerApexStrategyFamily,
    pub readiness: CornerApexStrategyReadiness,
    pub continuity_target: &'static str,
    pub deferred: bool,
    pub fallback_eligible: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstructionRiskKind {
    IntersectionDegenerate,
    TangencyUnsolved,
    BoundaryCurveMissing,
    SetbackSplitInfeasible,
    BoundaryConditionsIncomplete,
    TwistCompatibility,
    PatchBoundaryIncompatible,
    CurvatureExplosion,
    SelfIntersection,
    VariableRadiusCollapse,
    ChamferRoundMeshing,
    FacetedSmoothMismatch,
    ApproximationResidual,
    NumericalSolverFailure,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstructionRiskState {
    NoRisk,
    ManagedRetry(&'static str),
    RequiresMissingPrimitive(&'static str),
    Invalid(&'static str),
    UnknownInputFact,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ConstructionRiskFact {
    pub kind: ConstructionRiskKind,
    pub state: ConstructionRiskState,
    pub tolerance: DiagnosticTolerance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IntegrationValidationKind {
    Watertightness,
    EdgeOnSurface,
    VertexClosure,
    LocalSelfIntersection,
    AdjacentSurfaceIntersection,
    ToleranceExceeded,
    Continuity,
    SemanticHistory,
    MergeReadiness,
    FaceTracking,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IntegrationValidationState {
    Ready,
    RepairRequired(&'static str),
    DowngradeEligible(&'static str),
    Invalid(&'static str),
    UnknownInputFact,
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct IntegrationValidationFact {
    pub kind: IntegrationValidationKind,
    pub state: IntegrationValidationState,
    pub tolerance: DiagnosticTolerance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LifecycleState {
    ScopedDerivedDiscardable,
    StaleGraph,
    TopologyFingerprintMismatch,
    ReusedAfterMutation,
    ArtifactNeedsCleanup(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticArtifactKind {
    AcceptedGraphFact,
    DiagnosticArtifact,
    DeferredExploration,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LifecycleFact {
    pub state: LifecycleState,
    pub artifact: DiagnosticArtifactKind,
}
