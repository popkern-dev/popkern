//! Topological entities for boundary representation.

/// Edge topology entity.
pub mod edge;
pub mod edge_pcurve_reauthor;
pub mod edge_seam_certificate;
/// Face topology entity.
pub mod face;
/// Half-edge topology entity.
pub mod half_edge;
pub mod handles;
/// Loop topology entity.
pub mod loop_;
/// Semantic ownership role for topology entities.
pub mod role;
/// Shell topology entity.
pub mod shell;
/// Solid topology entity.
pub mod solid;
/// Vertex topology entity.
pub mod vertex;

pub use edge::{CanonicalForm, Edge};
pub use face::Face;
pub use half_edge::HalfEdge;
pub use handles::*;
pub use loop_::Loop;
pub use role::TopologyRole;
pub use shell::Shell;
pub use solid::Solid;
pub use vertex::Vertex;

pub mod brep;
pub use brep::BRep;
pub use brep::{SurfaceBoundaryHintDecline, SurfaceBoundaryHintDeclineReason};
pub mod coherence;
pub use coherence::{
    CoherenceReport, CoherenceValidationScope, L2Disagreement, PendingCoherenceMarker,
};
mod brep_face_use_trim_storage;
mod brep_loop_sampling;
mod validate_coherence;
pub use brep_face_use_trim_storage::FaceUseTrimKey;
mod surface_boundary_audit;
pub use surface_boundary_audit::{
    SurfaceBoundaryAuditEvent, SurfaceBoundaryAuditOutcome, SurfaceBoundarySharedTrimFailureReason,
    SurfaceBoundarySharedTrimValidation,
};

/// Non-manifold radial-cycle helpers for [`BRep`].
mod brep_radial;

/// Non-manifold repair policy and user-facing explanations.
pub mod brep_nonmanifold_policy;
pub use brep_nonmanifold_policy::{
    BrepEdgePairPolicyContext, BrepEdgePairPolicyDecision, BrepEdgePairPolicyDeclineReason,
    BrepEdgePairPolicyOutcome, BrepEdgePolicyState, BrepTopologyRepairPolicy,
    assess_edge_pair_topology_policy,
};

mod edge_tolerance_measure;
mod loop_chart_continuity;
mod pcurve_end_snap;
mod pcurve_image_check;
mod pcurve_loop_winding_repair;
mod pcurve_period_branch;
/// Pcurve recomputation for BRep faces.
mod pcurve_projection_fit;
mod pcurve_projection_walk;
mod pcurve_recompute;
mod pcurve_recompute_analytic;
mod pcurve_recompute_helpers;

/// Identity queries and recovery operations for precision tracking.
mod brep_precision_queries;

mod brep_edge_representatives;
/// Geometric deviation measurement for edges and vertices.
mod brep_measurement;
/// Whole-body orientation reversal (material-side mirror).
mod brep_reverse_orientation;
pub use brep_measurement::{EdgePcurveDeviationReport, HalfEdgePcurveDeviation};
mod brep_face_use_trim_authority;
mod brep_face_use_trim_authority_facts;
mod brep_face_use_trim_authority_repair;
pub use brep_face_use_trim_authority::{
    BrepFaceUseTrimAuthorityFinding, BrepFaceUseTrimAuthorityFindingKind,
    BrepFaceUseTrimAuthorityOptions, BrepFaceUseTrimAuthorityReport,
    analyze_brep_face_use_trim_authority,
};
pub use brep_face_use_trim_authority_repair::{
    BrepFaceUseTrimAuthorityRepairFailure, BrepFaceUseTrimAuthorityRepairReport,
    BrepFaceUseTrimAuthorityRepairResult, repair_brep_face_use_trim_authority_transactionally,
};

/// Transactional construction primitive for shared trim-edge geometry.
mod shared_trim_edge;
pub use shared_trim_edge::{
    SharedTrimEdgeCommit, SharedTrimEdgeFailure, SharedTrimEdgeFailureReason, SharedTrimEdgeInput,
    SharedTrimEdgeValidation, commit_shared_trim_edge_transactionally, validate_shared_trim_edge,
};

mod brep_inner_loop_rehome;
/// Edge and face splitting operations.
/// Sorted-pair serde for handle-keyed maps (JSON-safe, deterministic).
pub(crate) mod serde_handle_map;
pub mod split;
mod split_boundary;
pub(crate) mod split_helpers;
mod split_inner_loop;
pub use split::{SplitEdgeResult, SplitFaceResult};

/// Inner loop (hole) creation for face trimming.
pub mod trim;
pub use trim::AddInnerLoopResult;

/// Vertex-welding primitives for [`BRep`].
mod brep_weld;

pub mod brep_sewing;
/// Transactional BRep sewing for caller-supplied clean edge pairs.
mod brep_sewing_facts;
pub use brep_sewing::{
    BrepSewingDecline, BrepSewingDeclineReason, BrepSewingEdgePair, BrepSewingFailure,
    BrepSewingFailureReason, BrepSewingOptions, BrepSewingProvenance, BrepSewingReport,
    BrepSewingResult, edge_segment, sew_clean_boundary_edge_pairs,
};

/// Read-only near-edge classification for BRep sewing.
pub mod brep_sewing_classification;
pub use brep_sewing_classification::{
    BrepClassifiedEdgePair, BrepEdgePairClassification, BrepEdgePairMeasurement,
    BrepSewingClassificationBucket, BrepSewingClassificationBucketCount,
    BrepSewingClassificationOptions, BrepSewingClassificationReport,
    brep_sewing_classification_bucket_label, classify_boundary_edge_pairs,
};

/// Transactional tolerant snap-and-sew for classified BRep edge pairs.
pub mod brep_sewing_snap;
pub use brep_sewing_snap::{
    BrepSnapSewingDecline, BrepSnapSewingDeclineReason, BrepSnapSewingFailure,
    BrepSnapSewingFailureReason, BrepSnapSewingOptions, BrepSnapSewingReport, BrepSnapSewingResult,
    snap_and_sew_classified_edge_pairs,
};

mod brep_repair_actions;
/// Generic BRep boundary repair pipeline coordinator.
pub mod brep_repair_pipeline;
pub use brep_repair_actions::{
    BrepRepairActionArtifact, BrepRepairActionDecision, BrepRepairActionOutcome,
    BrepRepairActionRecord, BrepRepairActionReport, BrepRepairActionSummary,
    brep_repair_action_outcome_label, brep_repair_action_stage_label,
    build_brep_repair_action_artifact, repair_action_handle_label,
};
pub use brep_repair_pipeline::{
    BrepRepairPipelineFailure, BrepRepairPipelineFailureReason, BrepRepairPipelineOptions,
    BrepRepairPipelineReport, BrepRepairPipelineResult, BrepRepairPipelineStage,
    BrepTopologyOutputPolicy, repair_brep_boundary_edges_transactionally,
};

/// Shared read-only BRep repair knowledge vocabulary.
pub mod brep_repair_knowledge;
mod brep_repair_knowledge_api;
mod brep_repair_knowledge_authority;
mod brep_repair_knowledge_collect;
mod brep_repair_knowledge_diagnostics;
mod brep_repair_knowledge_endpoint;
mod brep_repair_knowledge_export;
mod brep_repair_knowledge_fact_build;
mod brep_repair_knowledge_fact_record;
mod brep_repair_knowledge_graph;
mod brep_repair_knowledge_loop;
mod brep_repair_knowledge_trim;
pub use brep_repair_knowledge::{
    BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION, BrepKnowledgeFingerprint, BrepRepairDiagnostic,
    BrepRepairDiagnosticCode, BrepRepairDiagnosticExport, BrepRepairDiagnosticId,
    BrepRepairDiagnosticSeverity, BrepRepairEvaluationLevel, BrepRepairFact,
    BrepRepairFactDependencies, BrepRepairFactHeader, BrepRepairFactId, BrepRepairFactKind,
    BrepRepairFactState, BrepRepairKnowledgeGraph, BrepRepairKnowledgeOptions,
    BrepRepairKnowledgeSchema, BrepRepairKnowledgeScope, BrepRepairProvenanceLabel,
    BrepRepairScopeExpansion, BrepRepairScopeHandles, BrepRepairScopeSeed,
    BrepRepairSuggestedRepair, BrepRepairTopologyHandle, LocalToleranceEvidence,
    LocalToleranceSource, SemanticFeatureId,
};
pub use brep_repair_knowledge_api::analyze_brep_repair_knowledge;
pub use brep_repair_knowledge_authority::{
    BrepRepairAuthorityFact, BrepRepairAuthorityKind, BrepRepairAuthorityReason,
};
pub use brep_repair_knowledge_collect::{
    collect_brep_adjacent_faces, collect_brep_endpoint_memberships, collect_brep_face_loops,
    collect_brep_loop_adjacencies, collect_brep_loop_facts,
};
pub use brep_repair_knowledge_endpoint::{BrepEndpointLoopMembership, BrepEndpointMembershipFacts};
pub use brep_repair_knowledge_export::{
    BrepRepairDiagnosticRecord, BrepRepairFactRecord, BrepRepairKnowledgeArtifact,
    BrepRepairKnowledgeSchemaRecord, BrepRepairKnowledgeScopeRecord, BrepRepairMeasurementRecord,
    BrepRepairMigrationLabelRecord, BrepRepairToleranceRecord,
    build_brep_repair_knowledge_artifact, build_brep_repair_knowledge_artifact_with_samples,
    render_brep_repair_diagnostic,
};
pub use brep_repair_knowledge_loop::{BrepLoopAdjacencyFacts, BrepLoopFacts, BrepLoopKind};
pub use brep_repair_knowledge_trim::{
    BrepTrimEndpointAgreementFacts, BrepTrimEndpointAgreementStatus,
    collect_brep_trim_endpoint_agreements, collect_brep_trim_endpoint_agreements_with_tolerance,
};

/// Transactional completion of missing BRep edge curves.
pub mod brep_edge_curve_completion;
mod brep_edge_curve_completion_facts;
mod brep_edge_curve_completion_helpers;
pub use brep_edge_curve_completion::{
    BrepEdgeCurveCompletionDecision, BrepEdgeCurveCompletionDeclineReason,
    BrepEdgeCurveCompletionFailure, BrepEdgeCurveCompletionFailureReason,
    BrepEdgeCurveCompletionOptions, BrepEdgeCurveCompletionOutcome, BrepEdgeCurveCompletionReport,
    BrepEdgeCurveCompletionResult, complete_missing_edge_curves_transactionally,
};

/// Transactional pcurve reconstruction from authoritative 3D edge curves.
pub mod brep_curve_trim_reconstruction;
mod brep_curve_trim_reconstruction_domain;
mod brep_curve_trim_reconstruction_facts;
mod brep_curve_trim_reconstruction_helpers;
pub use brep_curve_trim_reconstruction::{
    BrepCurveTrimReconstructionDecision, BrepCurveTrimReconstructionDeclineReason,
    BrepCurveTrimReconstructionFailure, BrepCurveTrimReconstructionFailureReason,
    BrepCurveTrimReconstructionOptions, BrepCurveTrimReconstructionOutcome,
    BrepCurveTrimReconstructionReport, BrepCurveTrimReconstructionResult,
    reconstruct_curve_backed_trims_transactionally,
};

/// Direct edge-curve to adjacent face-surface diagnostics and repair.
pub mod brep_edge_surface_compatibility;
mod brep_edge_surface_compatibility_boundary;
mod brep_edge_surface_compatibility_candidate_facts;
mod brep_edge_surface_compatibility_commit;
mod brep_edge_surface_compatibility_decline;
mod brep_edge_surface_compatibility_diagnostics;
mod brep_edge_surface_compatibility_helpers;
mod brep_edge_surface_compatibility_pcurve_trace;
mod brep_edge_surface_compatibility_repair;
mod brep_edge_surface_compatibility_scope;
mod brep_edge_surface_compatibility_shared_facts;
mod brep_edge_surface_compatibility_single_surface;
mod brep_edge_surface_compatibility_trim_candidate;
mod brep_edge_surface_compatibility_vertex_candidate;
mod brep_edge_surface_compatibility_vertex_candidate_geometry;
mod brep_edge_surface_compatibility_vertex_star;
mod brep_edge_surface_compatibility_vertex_star_repair;
pub use brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityDeclineReason, BrepEdgeSurfaceCompatibilityFinding,
    BrepEdgeSurfaceCompatibilityOptions, BrepEdgeSurfaceCompatibilityRepairDecision,
    BrepEdgeSurfaceCompatibilityRepairFailure, BrepEdgeSurfaceCompatibilityRepairFailureReason,
    BrepEdgeSurfaceCompatibilityRepairOutcome, BrepEdgeSurfaceCompatibilityRepairReport,
    BrepEdgeSurfaceCompatibilityRepairResult, BrepEdgeSurfaceCompatibilityReport,
    analyze_brep_edge_surface_compatibility,
    repair_brep_edge_surface_compatibility_transactionally,
};
pub use brep_edge_surface_compatibility_diagnostics::{
    BrepEdgeSurfaceBoundaryDeclineEvidence, BrepEdgeSurfaceBoundaryHintEvidence,
    BrepEdgeSurfaceKind, BrepEdgeSurfaceRepairEvidence,
};
pub use brep_edge_surface_compatibility_scope::analyze_scoped_brep_edge_surface_compatibility;
pub use brep_edge_surface_compatibility_vertex_candidate::{
    BrepEdgeSurfaceVertexCandidateDiagnostic, BrepEdgeSurfaceVertexCandidateOutcome,
    BrepEdgeSurfaceVertexCandidateSource,
};
pub use brep_edge_surface_compatibility_vertex_star::{
    BrepEdgeSurfaceVertexStarDecision, BrepEdgeSurfaceVertexStarDeclineReason,
    BrepEdgeSurfaceVertexStarOutcome, BrepEdgeSurfaceVertexStarReport,
};

/// Read-only degenerate edge diagnostics.
pub mod brep_degenerate_edges;
pub use brep_degenerate_edges::{
    BrepDegenerateEdgeFinding, BrepDegenerateEdgeKind, BrepDegenerateEdgeOptions,
    BrepDegenerateEdgeReport, analyze_brep_degenerate_edges,
};

/// Read-only BRep export-quality diagnostics.
pub mod brep_export_quality;
mod brep_export_quality_scope;
mod brep_fidelity_blocker_ownership;
pub use brep_export_quality::{
    BrepBoundaryEdgeArtifact, BrepExportBlockerSummary, BrepExportQualityArtifact,
    BrepExportQualityOptions, BrepExportQualityReport, BrepExportQualitySummary,
    BrepExportReadiness, BrepSurfacelessFaceRecord, analyze_brep_export_quality,
    analyze_brep_export_quality_with_local_diagnostics,
    analyze_brep_export_quality_with_pcurve_fidelity, build_brep_export_quality_artifact,
};
pub use brep_export_quality_scope::{
    BrepExportQualityScope, BrepScopedExportQualityReport, analyze_scoped_brep_export_quality,
};
pub use brep_fidelity_blocker_ownership::{
    BrepFidelityBlockerOwnership, BrepFidelityBlockerOwnershipReport,
    BrepFidelityBlockerSideOwnership, classify_brep_fidelity_blocker_ownership,
};

/// Read-only loop and shell diagnostics for repair decisions.
pub mod brep_shell_diagnostics;
mod brep_shell_diagnostics_facts;
mod brep_shell_diagnostics_helpers;
pub use brep_shell_diagnostics::{
    BrepLoopShellDiagnosticReport, BrepTopologyDiagnosticFinding, BrepTopologyDiagnosticKind,
    BrepTopologyDiagnosticSeverity, analyze_brep_loop_shell_validity,
};

/// Trim and pcurve diagnostics plus transactional repair.
pub mod brep_trim_repair;
pub use brep_trim_repair::{
    BrepTrimDiagnosticReport, BrepTrimFinding, BrepTrimFindingKind, BrepTrimRepairFailure,
    BrepTrimRepairFailureReason, BrepTrimRepairOptions, BrepTrimRepairReport, BrepTrimRepairResult,
    analyze_brep_trims, analyze_scoped_brep_trims, repair_brep_trims_transactionally,
};
