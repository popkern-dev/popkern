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
        ],
        ConstructionFamily::PairedSameLoopTwinExposure => vec![executable(
            "same-loop-boundary-exposure-dependency-satisfied",
            "same-loop-boundary-exposure",
        )],
        ConstructionFamily::BoundarySideSeparatorBridge => vec![
            missing(
                "boundary-side-separator-subcomplex-primitive-missing",
                "boundary-side-separator-subcomplex",
            ),
            executable(
                "boundary-side-separator-projection-dependency-ready",
                "surface-projection",
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
