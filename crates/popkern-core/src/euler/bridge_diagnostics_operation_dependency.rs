//! Neighboring operation dependency facts for bridge diagnostic candidates.

use crate::euler::{
    ConstructionFamily, OperationDependency, OperationDependencyReadiness,
    OperationDependencyStatus,
};

pub(super) fn dependencies_for_candidate(family: ConstructionFamily) -> Vec<OperationDependency> {
    match family {
        ConstructionFamily::DisjointCrossLoopBridge => vec![
            executable(
                "bridge-surface-projection-dependency-ready",
                "surface-projection",
            ),
            executable("bridge-section-trim-dependency-ready", "section-trim"),
            executable("bridge-filling-dependency-ready", "face-filling"),
            executable(
                "bridge-history-generated-entities-expected",
                "generated-entity-history",
            ),
            executable(
                "bridge-history-modified-supports-expected",
                "modified-entity-history",
            ),
        ],
        ConstructionFamily::DisjointSameLoopBridge => vec![
            executable(
                "same-loop-bridge-surface-projection-dependency-ready",
                "surface-projection",
            ),
            executable(
                "same-loop-bridge-section-trim-dependency-ready",
                "section-trim",
            ),
            executable(
                "same-loop-bridge-history-generated-entities-expected",
                "generated-entity-history",
            ),
            executable(
                "same-loop-bridge-history-modified-supports-expected",
                "modified-entity-history",
            ),
        ],
        ConstructionFamily::PairedSameLoopTwinExposure => vec![
            executable(
                "same-loop-boundary-exposure-dependency-satisfied",
                "same-loop-boundary-exposure",
            ),
            executable(
                "same-loop-boundary-exposure-history-reclassification-expected",
                "reclassified-boundary-history",
            ),
        ],
        ConstructionFamily::BoundarySideSeparatorBridge => vec![
            missing(
                "boundary-side-separator-subcomplex-primitive-missing",
                "boundary-side-separator-subcomplex",
            ),
            executable(
                "boundary-side-separator-projection-dependency-ready",
                "surface-projection",
            ),
            executable(
                "boundary-side-separator-history-reclassification-expected",
                "reclassified-boundary-history",
            ),
            executable(
                "boundary-side-separator-history-generated-entities-expected",
                "generated-entity-history",
            ),
        ],
        ConstructionFamily::Unknown => vec![unknown()],
    }
}

fn executable(reason: &'static str, primitive: &'static str) -> OperationDependency {
    OperationDependency {
        readiness: OperationDependencyReadiness::Executable,
        status: OperationDependencyStatus::Satisfied,
        reason,
        primitive,
    }
}

fn missing(reason: &'static str, primitive: &'static str) -> OperationDependency {
    OperationDependency {
        readiness: OperationDependencyReadiness::MissingPrimitive,
        status: OperationDependencyStatus::RequiredMissingPrimitive,
        reason,
        primitive,
    }
}

fn unknown() -> OperationDependency {
    OperationDependency {
        readiness: OperationDependencyReadiness::Unknown,
        status: OperationDependencyStatus::Unknown,
        reason: "operation-dependency-not-classified",
        primitive: "unknown",
    }
}
