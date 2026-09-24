//! Blend-decline projection: the per-variant and per-edge records of
//! [`FeatureFailure::BlendDeclined`](crate::failure::FeatureFailure).

use popkern_blend::error::{BlendError, BlendFailure};
use popkern_blend::remediation::Remediation;
use popkern_core::naming::EntityName;
use popkern_core::topology::{BRep, EdgeHandle};
use serde::{Deserialize, Serialize};

/// Payload of [`FeatureFailure::BlendDeclined`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlendDeclined {
    /// The `BlendError` display text.
    pub message: String,
    /// Engine remediation for the whole decline.
    pub remediation: RemediationRecord,
    /// Variant-specific data, tagged by `blend_kind`.
    #[serde(flatten)]
    pub decline: BlendDeclineKind,
}

/// Variant-specific payload of a blend decline.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "blend_kind")]
pub enum BlendDeclineKind {
    /// Two fillet strips overlap mid-face.
    FilletInterference {
        /// Spec index of the first fillet.
        fillet_a: usize,
        /// Spec index of the second fillet.
        fillet_b: usize,
    },
    /// A face-blend attach declined.
    AttachDeclined {
        /// Surgery stage that declined.
        mechanism: String,
        /// Measured evidence.
        detail: String,
    },
    /// A parameter was invalid.
    InvalidParameter {
        /// Kernel spec path.
        parameter: String,
        /// Rendered value.
        value: String,
        /// Reason.
        reason: String,
    },
    /// Every edge failed pre-analysis.
    NoValidEdges {
        /// Detail.
        detail: String,
    },
    /// A core topology error.
    Topology {
        /// Detail.
        detail: String,
    },
    /// Some edges failed; the operation declined transactionally.
    PartialFailure {
        /// Distinct failed edges.
        failed: usize,
        /// Edges attempted.
        total: usize,
        /// One record per failed edge.
        failures: Vec<BlendFailureRecord>,
    },
    /// Cancelled by the caller.
    Cancelled,
}

/// One failed edge of a partial blend failure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlendFailureRecord {
    /// Debug form of the kernel edge handle.
    pub edge: String,
    /// The edge's stable name, when the body carries one.
    pub edge_name: Option<EntityName>,
    /// Why this edge failed, tagged by `reason_kind`.
    #[serde(flatten)]
    pub reason: BlendFailureReason,
    /// Engine remediation for this edge.
    pub remediation: RemediationRecord,
}

/// Per-edge failure reason (projection of `BlendFailure`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "reason_kind")]
pub enum BlendFailureReason {
    /// The radius does not fit.
    RadiusTooLarge {
        /// Requested radius.
        requested: f64,
        /// Largest radius the engine estimates would fit.
        max_radius: Option<f64>,
        /// Where it was detected.
        location: Option<[f64; 3]>,
    },
    /// Surface offset failed.
    OffsetFailed {
        /// Detail.
        reason: String,
        /// Where it was detected.
        location: Option<[f64; 3]>,
    },
    /// A vertex blend failed.
    VertexBlendFailed {
        /// Debug form of the vertex handle.
        vertex: String,
        /// Detail.
        reason: String,
        /// Where it was detected.
        location: Option<[f64; 3]>,
    },
    /// Trim intersection failed.
    IntersectionFailed {
        /// Detail.
        reason: String,
        /// Where it was detected.
        location: Option<[f64; 3]>,
    },
    /// The edge's topology is invalid for blending.
    InvalidTopology {
        /// Detail.
        reason: String,
        /// Where it was detected.
        location: Option<[f64; 3]>,
    },
    /// The rolling-ball contact solver did not converge.
    ContactSolverFailed {
        /// Detail.
        reason: String,
        /// Spine parameter at the breakdown.
        spine_parameter: Option<f64>,
        /// Where it was detected.
        location: Option<[f64; 3]>,
    },
    /// The edge is already blended.
    AlreadyBlended {
        /// Where.
        location: Option<[f64; 3]>,
    },
    /// The edge is a seam.
    SeamEdge {
        /// Where.
        location: Option<[f64; 3]>,
    },
    /// The constructed blend surface self-intersects.
    SelfIntersectingBlend {
        /// Detail.
        reason: String,
        /// Where it was detected.
        location: Option<[f64; 3]>,
    },
    /// The served blend territory has an edge off its faces or an
    /// over-ceiling band.
    InconsistentGeometry {
        /// Detail.
        reason: String,
        /// Where the worst edge is.
        location: Option<[f64; 3]>,
    },
}

/// Projection of `popkern_blend::remediation::Remediation`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RemediationRecord {
    /// `RemediationAction` variant name.
    pub action: String,
    /// Tolerance parameter to adjust, if any.
    pub parameter: Option<String>,
    /// Suggested value, if any.
    pub value: Option<f64>,
    /// Why it should help.
    pub rationale: String,
}

impl From<&Remediation> for RemediationRecord {
    fn from(r: &Remediation) -> Self {
        RemediationRecord {
            action: format!("{:?}", r.action),
            parameter: r.parameter.clone(),
            value: r.value,
            rationale: r.rationale.clone(),
        }
    }
}

fn point(p: &Option<popkern_core::geometry::Point3>) -> Option<[f64; 3]> {
    p.map(|p| [p.x, p.y, p.z])
}

impl From<&BlendFailure> for BlendFailureReason {
    fn from(f: &BlendFailure) -> Self {
        match f {
            BlendFailure::RadiusTooLarge {
                requested,
                max_radius,
                location,
            } => BlendFailureReason::RadiusTooLarge {
                requested: *requested,
                max_radius: *max_radius,
                location: point(location),
            },
            BlendFailure::OffsetFailed { reason, location } => BlendFailureReason::OffsetFailed {
                reason: reason.clone(),
                location: point(location),
            },
            BlendFailure::VertexBlendFailed {
                vertex,
                reason,
                location,
            } => BlendFailureReason::VertexBlendFailed {
                vertex: format!("{vertex:?}"),
                reason: reason.clone(),
                location: point(location),
            },
            BlendFailure::IntersectionFailed { reason, location } => {
                BlendFailureReason::IntersectionFailed {
                    reason: reason.clone(),
                    location: point(location),
                }
            }
            BlendFailure::InvalidTopology { reason, location } => {
                BlendFailureReason::InvalidTopology {
                    reason: reason.clone(),
                    location: point(location),
                }
            }
            BlendFailure::ContactSolverFailed {
                reason,
                spine_parameter,
                location,
            } => BlendFailureReason::ContactSolverFailed {
                reason: reason.clone(),
                spine_parameter: *spine_parameter,
                location: point(location),
            },
            BlendFailure::AlreadyBlended { location } => BlendFailureReason::AlreadyBlended {
                location: point(location),
            },
            BlendFailure::SeamEdge { location } => BlendFailureReason::SeamEdge {
                location: point(location),
            },
            BlendFailure::SelfIntersectingBlend { reason, location } => {
                BlendFailureReason::SelfIntersectingBlend {
                    reason: reason.clone(),
                    location: point(location),
                }
            }
            BlendFailure::InconsistentGeometry { reason, location } => {
                BlendFailureReason::InconsistentGeometry {
                    reason: reason.clone(),
                    location: point(location),
                }
            }
        }
    }
}

impl BlendFailureRecord {
    /// Project one per-edge failure; `brep` supplies the stable name
    /// when the edge carries one.
    pub fn new(edge: EdgeHandle, failure: &BlendFailure, brep: Option<&BRep>) -> Self {
        let edge_name = brep.and_then(|b| {
            b.edges
                .get(edge)
                .and_then(|e| e.name)
                .and_then(|nid| b.name_pool.get(nid).cloned())
        });
        BlendFailureRecord {
            edge: format!("{edge:?}"),
            edge_name,
            reason: BlendFailureReason::from(failure),
            remediation: RemediationRecord::from(&failure.remediation_hint()),
        }
    }
}

impl BlendDeclined {
    /// Project a blend decline; `brep` is the body the blend ran on
    /// (for stable edge names in `failures[]`).
    pub fn project(error: &BlendError, brep: Option<&BRep>) -> Self {
        let decline = match error {
            BlendError::FilletInterference { fillet_a, fillet_b } => {
                BlendDeclineKind::FilletInterference {
                    fillet_a: *fillet_a,
                    fillet_b: *fillet_b,
                }
            }
            BlendError::AttachDeclined { mechanism, detail } => BlendDeclineKind::AttachDeclined {
                mechanism: mechanism.clone(),
                detail: detail.clone(),
            },
            BlendError::InvalidParameter {
                parameter,
                value,
                reason,
            } => BlendDeclineKind::InvalidParameter {
                parameter: parameter.clone(),
                value: value.clone(),
                reason: reason.clone(),
            },
            BlendError::NoValidEdges(s) => BlendDeclineKind::NoValidEdges { detail: s.clone() },
            BlendError::Topology(c) => BlendDeclineKind::Topology {
                detail: c.to_string(),
            },
            BlendError::PartialFailure {
                failed,
                total,
                failures,
                ..
            } => BlendDeclineKind::PartialFailure {
                failed: *failed,
                total: *total,
                failures: failures
                    .iter()
                    .map(|(h, f)| BlendFailureRecord::new(*h, f, brep))
                    .collect(),
            },
            BlendError::Cancelled => BlendDeclineKind::Cancelled,
        };
        BlendDeclined {
            message: error.to_string(),
            remediation: RemediationRecord::from(&error.remediation_hint()),
            decline,
        }
    }
}
