//! Downstream readiness facts for bridge diagnostic candidates.

use crate::euler::{
    ConstructionFamily, ConstructionFeasibility, DownstreamReadiness, DownstreamReadinessFact,
    DownstreamReadinessStatus,
};

pub(super) fn downstream_facts_for_candidate(
    family: ConstructionFamily,
    feasibility: ConstructionFeasibility,
) -> Vec<DownstreamReadinessFact> {
    match feasibility {
        ConstructionFeasibility::ExecutableByExistingOperator
        | ConstructionFeasibility::ExecutableAfterKnownNormalization => {
            vec![ready_for_family(family)]
        }
        ConstructionFeasibility::RequiresSeparatorSubcomplex
        | ConstructionFeasibility::RequiresCustomAggregateAuthoring
        | ConstructionFeasibility::RequiresCycleRegrouping => {
            vec![blocked_by_missing_primitive(family)]
        }
        ConstructionFeasibility::InvalidNonManifold => vec![invalid_for_export_or_mesh(family)],
        ConstructionFeasibility::Unknown | ConstructionFeasibility::UnsupportedUnknownFamily => {
            vec![unknown()]
        }
    }
}

fn ready_for_family(family: ConstructionFamily) -> DownstreamReadinessFact {
    DownstreamReadinessFact {
        readiness: DownstreamReadiness::Ready,
        status: DownstreamReadinessStatus::Ready,
        reason: match family {
            ConstructionFamily::DisjointCrossLoopBridge => {
                "cross-loop-bridge-downstream-consumers-ready"
            }
            ConstructionFamily::DisjointSameLoopBridge => {
                "same-loop-bridge-downstream-consumers-ready"
            }
            ConstructionFamily::PairedSameLoopTwinExposure => {
                "same-loop-exposure-downstream-consumers-ready"
            }
            ConstructionFamily::BoundarySideSeparatorBridge => {
                "boundary-side-separator-downstream-consumers-ready"
            }
            ConstructionFamily::Unknown => "unknown-candidate-downstream-not-classified",
        },
    }
}

fn blocked_by_missing_primitive(family: ConstructionFamily) -> DownstreamReadinessFact {
    DownstreamReadinessFact {
        readiness: DownstreamReadiness::BlockedByMissingPrimitive,
        status: DownstreamReadinessStatus::Blocked,
        reason: match family {
            ConstructionFamily::BoundarySideSeparatorBridge => {
                "boundary-side-separator-downstream-blocked-by-missing-primitive"
            }
            ConstructionFamily::DisjointCrossLoopBridge => {
                "cross-loop-bridge-downstream-blocked-by-missing-primitive"
            }
            ConstructionFamily::DisjointSameLoopBridge => {
                "same-loop-bridge-downstream-blocked-by-missing-primitive"
            }
            ConstructionFamily::PairedSameLoopTwinExposure => {
                "same-loop-exposure-downstream-blocked-by-missing-primitive"
            }
            ConstructionFamily::Unknown => "unknown-candidate-downstream-not-classified",
        },
    }
}

fn invalid_for_export_or_mesh(family: ConstructionFamily) -> DownstreamReadinessFact {
    DownstreamReadinessFact {
        readiness: DownstreamReadiness::InvalidForExportOrMesh,
        status: DownstreamReadinessStatus::Blocked,
        reason: match family {
            ConstructionFamily::DisjointCrossLoopBridge => {
                "cross-loop-bridge-downstream-invalid-for-export-or-mesh"
            }
            ConstructionFamily::DisjointSameLoopBridge => {
                "same-loop-bridge-downstream-invalid-for-export-or-mesh"
            }
            ConstructionFamily::BoundarySideSeparatorBridge => {
                "boundary-side-separator-downstream-invalid-for-export-or-mesh"
            }
            ConstructionFamily::PairedSameLoopTwinExposure => {
                "same-loop-exposure-downstream-invalid-for-export-or-mesh"
            }
            ConstructionFamily::Unknown => "unknown-candidate-downstream-not-classified",
        },
    }
}

fn unknown() -> DownstreamReadinessFact {
    DownstreamReadinessFact {
        readiness: DownstreamReadiness::Unknown,
        status: DownstreamReadinessStatus::Unknown,
        reason: "candidate-downstream-readiness-not-classified",
    }
}
