//! Neighboring operation dependency facts for bridge diagnostic candidates.

use crate::euler::{
    ConstructionFamily, OperationDependency, OperationDependencyReadiness,
    OperationDependencyStatus,
};

pub(super) fn dependencies_for_candidate(family: ConstructionFamily) -> Vec<OperationDependency> {
    let (readiness, status, reason, primitive) = match family {
        ConstructionFamily::DisjointCrossLoopBridge
        | ConstructionFamily::DisjointSameLoopBridge => (
            OperationDependencyReadiness::NotApplicable,
            OperationDependencyStatus::NotRequired,
            "bridge-candidate-has-no-neighboring-operation-dependency",
            "none",
        ),
        ConstructionFamily::PairedSameLoopTwinExposure => (
            OperationDependencyReadiness::Executable,
            OperationDependencyStatus::Satisfied,
            "same-loop-boundary-exposure-dependency-satisfied",
            "same-loop-boundary-exposure",
        ),
        ConstructionFamily::BoundarySideSeparatorBridge => (
            OperationDependencyReadiness::MissingPrimitive,
            OperationDependencyStatus::RequiredMissingPrimitive,
            "boundary-side-separator-subcomplex-primitive-missing",
            "boundary-side-separator-subcomplex",
        ),
        ConstructionFamily::Unknown => (
            OperationDependencyReadiness::Unknown,
            OperationDependencyStatus::Unknown,
            "operation-dependency-not-classified",
            "unknown",
        ),
    };
    vec![OperationDependency {
        readiness,
        status,
        reason,
        primitive,
    }]
}
