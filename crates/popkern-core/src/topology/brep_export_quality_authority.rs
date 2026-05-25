//! Export-quality authority summary helpers.

use super::BrepExportAuthoritySummary;
use crate::topology::{
    BRep, BrepRepairAuthorityKind, BrepRepairAuthorityReason, BrepRepairFactKind,
    BrepRepairKnowledgeGraph, BrepRepairKnowledgeScope, BrepRepairScopeExpansion,
    BrepRepairScopeSeed, BrepRepairTopologyHandle, BrepTrimRepairOptions, EdgeHandle,
    LocalToleranceSource,
};
use std::collections::BTreeSet;

pub(super) fn export_authority_summary(
    brep: &BRep,
    trim: BrepTrimRepairOptions,
) -> BrepExportAuthoritySummary {
    let graph = export_authority_graph(brep, trim);
    let mut summary = BrepExportAuthoritySummary {
        facts: 0,
        stale_pcurves: 0,
        conflicting_sources: 0,
        missing_sources: 0,
        degenerate_topology: 0,
    };
    for fact in graph
        .facts
        .iter()
        .filter(|fact| fact.header.kind == BrepRepairFactKind::Authority)
    {
        let Some(authority) = fact.authority.as_ref() else {
            continue;
        };
        summary.facts += 1;
        if authority.kind == BrepRepairAuthorityKind::StalePcurve {
            summary.stale_pcurves += 1;
        }
        if authority.kind == BrepRepairAuthorityKind::SideLocalSurfaceMismatch {
            summary.conflicting_sources += 1;
        }
        if authority.kind == BrepRepairAuthorityKind::NoAuthoritativeSource
            && authority.reason != BrepRepairAuthorityReason::MissingSurfaceBoundaryEvidence
        {
            summary.missing_sources += 1;
        }
        if authority.kind == BrepRepairAuthorityKind::AuthoritativeDegenerateTopology {
            summary.degenerate_topology += 1;
        }
    }
    summary
}

fn export_authority_graph(brep: &BRep, trim: BrepTrimRepairOptions) -> BrepRepairKnowledgeGraph {
    let mut scope = BrepRepairKnowledgeScope::seeds_only(
        brep.edges
            .keys()
            .map(BrepRepairScopeSeed::Edge)
            .collect::<Vec<_>>(),
    );
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    BrepRepairKnowledgeGraph::build_with_local_tolerance(
        brep,
        scope,
        trim.endpoint_tolerance,
        LocalToleranceSource::RepairPolicy,
        trim.samples,
    )
}

pub(super) fn authoritative_degenerate_edges(
    brep: &BRep,
    trim: BrepTrimRepairOptions,
) -> BTreeSet<EdgeHandle> {
    export_authority_graph(brep, trim)
        .facts
        .iter()
        .filter(|fact| {
            fact.header.kind == BrepRepairFactKind::Authority
                && fact.authority.as_ref().is_some_and(|authority| {
                    authority.kind == BrepRepairAuthorityKind::AuthoritativeDegenerateTopology
                })
        })
        .flat_map(|fact| fact.handles.iter())
        .filter_map(|handle| match *handle {
            BrepRepairTopologyHandle::Edge(edge) => Some(edge),
            _ => None,
        })
        .collect()
}
