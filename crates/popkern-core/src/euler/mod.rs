//! Euler operators — the only way to modify BRep topology.
//!
//! Each operator maintains the Euler-Poincaré invariant by construction.

mod boundary_side_separator_bridge;
mod bridge_diagnostics;
mod construction_residual_retire;
mod contract;
mod cylinder_mediated_support;
mod cylinder_mediated_support_apply;
mod cylinder_mediated_support_partition;
mod face_boundary_fragment;
mod face_boundary_retire;
mod face_family_bridge;
mod face_merge;
mod face_merge_split;
mod kemr;
mod kfmrh;
mod loop_bridge;
mod loop_bridge_apply;
mod mef;
mod mev;
mod mvfs;
mod same_loop_boundary_exposure;
mod same_loop_bridge_dissolve;
mod same_loop_face_bridge;
mod sweep;
pub use boundary_side_separator_bridge::{
    BoundarySideSeparatorBridgePlan, BoundarySideSeparatorBridgeRequest,
    BoundarySideSeparatorBridgeTopologyDelta,
};
pub use bridge_diagnostics::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticReport, BridgeDiagnosticRequest,
    BridgeDiagnosticSummary, ConstructionCandidate, ConstructionCandidateId, ConstructionFamily,
    ConstructionFeasibility, ConstructionRequirement, DiagnosticHandle, DiagnosticOperation,
    DiagnosticOperationId, DiagnosticOperationKind, DiagnosticOperationPurpose,
    DiagnosticOperatorFamily, DiagnosticOwner, DiagnosticResource, DiagnosticResourceId,
    DiagnosticResourceKind, DiagnosticResourceStatus, OperatorReadiness, OperatorReadinessStatus,
    PostconditionPrediction, PostconditionPredictionId, ResourceConflict, ResourceConflictId,
    ResourceConflictKind, ResourceConflictStatus, TopologyDiagnosticGraph,
    TopologyGraphFingerprint,
};
pub use construction_residual_retire::ConstructionResidualRetireResult;
pub use cylinder_mediated_support::{
    CylinderMediatedOutsideIntervalInteraction, CylinderMediatedSupportContinuationAggregateDelta,
    CylinderMediatedSupportContinuationBatchPlan, CylinderMediatedSupportContinuationDelta,
    CylinderMediatedSupportContinuationPlan, CylinderMediatedSupportContinuationRequest,
};
pub use cylinder_mediated_support_apply::{
    CylinderMediatedContextOwnerFace, CylinderMediatedOwnerRepairSpan,
    CylinderMediatedSupportContinuationBatchResult,
    CylinderMediatedSupportContinuationPostconditions,
    CylinderMediatedSupportContinuationSegmentAuthoringPlan,
    CylinderMediatedSupportContinuationSegmentRequirement, CylinderMediatedUnsupportedOwnerFace,
};
pub use cylinder_mediated_support_partition::{
    CylinderMediatedOutsideIntervalPartition, CylinderMediatedOutsideIntervalPartitionSegment,
    CylinderMediatedOutsideLoopPartition,
};
pub use face_boundary_fragment::BoundaryFragmentRetireResult;
pub use face_boundary_retire::BoundaryFaceRetireResult;
pub use face_family_bridge::{
    FaceFamilyBridgeCycleRequest, FaceFamilyBridgeCycleResult, FaceFamilyBridgePlan,
    FaceFamilyBridgeResult,
};
pub use face_merge::FaceMergeResult;
pub use face_merge_split::{FaceMergeSplitRequest, FaceMergeSplitResult};
pub use kemr::{KemrResult, MekrResult};
pub use kfmrh::MfkrhResult;
pub use loop_bridge::{TwoLoopBridgePlan, TwoLoopBridgeRequest};
pub use loop_bridge_apply::TwoLoopBridgeResult;
pub use mef::MefResult;
pub use mev::MevResult;
pub use mvfs::MvfsResult;
pub use same_loop_boundary_exposure::SameLoopBoundaryExposureResult;
pub use same_loop_face_bridge::{
    SameLoopFaceBridgePlan, SameLoopFaceBridgeRequest, SameLoopFaceBridgeResult,
};
pub use sweep::ExtrudeResult;
pub use sweep::RevolveResult;
