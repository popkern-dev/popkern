//! Durable surface-boundary audit records.

use serde::{Deserialize, Serialize};

use crate::geometry::edge_triple::BoundarySide;
use crate::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle, SurfaceBoundaryHintDeclineReason};

/// Durable outcome recorded by construction or repair when auditing
/// canonical surface-boundary evidence for a half-edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SurfaceBoundaryAuditOutcome {
    /// Existing boundary hint was revalidated against current trim evidence.
    FreshHint,
    /// Existing boundary hint no longer matched current trim/topology evidence.
    StaleHintCleared,
    /// A previous decline was revisited by a later audit pass.
    ReauditDecline,
    /// Exact shared trim edge was validated and committed.
    SharedTrimAccepted,
    /// Exact shared trim edge was declined by validation.
    SharedTrimDeclined,
    /// Pcurve classification accepted a canonical boundary.
    ClassifiedHint,
    /// Pcurve classification declined canonical boundary ownership.
    ClassifiedDecline,
}

/// Stable copy of shared-trim validation failure reasons for BRep audit
/// records. This avoids making audit consumers depend on the construction
/// primitive's internal error type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SurfaceBoundarySharedTrimFailureReason {
    /// Edge, half-edge, loop, face, vertex, or surface topology was missing.
    MissingTopology,
    /// The supplied half-edges did not both belong to the target edge.
    HalfEdgeEdgeMismatch,
    /// The supplied half-edges were not twins.
    HalfEdgesNotTwins,
    /// Curve endpoints did not agree with edge topology.
    EndpointMismatch,
    /// One or both pcurves did not agree with the 3-D curve on their surfaces.
    CurveSurfaceMismatch,
    /// The two pcurve/surface images did not describe the same 3-D locus.
    SideAgreementMismatch,
    /// Repaired clone failed BRep validation.
    ValidationFailed,
}

/// Numeric validation evidence from a shared trim-edge audit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SurfaceBoundarySharedTrimValidation {
    /// Maximum sampled curve-to-pcurve/surface deviation.
    pub max_curve_surface_deviation: f64,
    /// Maximum sampled side-A to side-B pcurve/surface deviation.
    pub max_side_agreement_deviation: f64,
    /// Maximum endpoint gap against topology vertices.
    pub max_endpoint_gap: f64,
}

/// One durable surface-boundary audit event for diagnostics and repair
/// knowledge promotion.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SurfaceBoundaryAuditEvent {
    /// Half-edge being audited.
    pub half_edge: HalfEdgeHandle,
    /// Owning edge when available.
    pub edge: Option<EdgeHandle>,
    /// Owning face of `half_edge` when available.
    pub face: Option<FaceHandle>,
    /// Adjacent twin face when available.
    pub twin_face: Option<FaceHandle>,
    /// Audit outcome.
    pub outcome: SurfaceBoundaryAuditOutcome,
    /// Boundary side accepted or evaluated, when known.
    pub boundary: Option<BoundarySide>,
    /// Decline reason from canonical-boundary classification, when present.
    pub decline_reason: Option<SurfaceBoundaryHintDeclineReason>,
    /// Shared trim-edge failure reason, when shared trim validation declined.
    pub shared_trim_failure: Option<SurfaceBoundarySharedTrimFailureReason>,
    /// Shared trim-edge numeric validation evidence, when measured.
    pub shared_trim_validation: Option<SurfaceBoundarySharedTrimValidation>,
    /// Local tolerance used by the audit, when applicable.
    pub tolerance: Option<f64>,
    /// Construction or repair stage that recorded the audit.
    pub provenance: String,
}
