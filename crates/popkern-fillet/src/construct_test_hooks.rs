//! Test-only entry points into the construct pipeline (PK-097):
//! PnS setback candidates and the setback-topology commit switches.

/// PK-097: the full PnS setback candidate (patches + surfaces).
pub fn pns_setback_candidate_full(
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    vertex_position: popkern_core::geometry::Point3,
) -> Result<crate::construct::pns_setback_candidate::PnsSetbackCandidate, String> {
    crate::construct::pns_setback_candidate::try_pns_setback_candidate(
        fillet_surfaces,
        vertex_position,
    )
}

/// PK-097 Phase A: run the PnS setback candidate (G1 cap on the
/// setback border) and return its quality report.
pub fn pns_setback_candidate(
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    vertex_position: popkern_core::geometry::Point3,
) -> Result<crate::construct::pns_setback_candidate::PnsSetbackReport, String> {
    crate::construct::pns_setback_candidate::try_pns_setback_candidate(
        fillet_surfaces,
        vertex_position,
    )
    .map(|c| c.report)
}

/// Force COMMIT mode on (tests; production defaults to on).
pub fn enable_setback_topology_commit() {
    crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT
        .store(true, std::sync::atomic::Ordering::Relaxed);
    crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT_SET
        .store(true, std::sync::atomic::Ordering::Relaxed);
}

/// Force COMMIT mode off (tests needing the legacy corner path).
pub fn disable_setback_topology_commit() {
    crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT
        .store(false, std::sync::atomic::Ordering::Relaxed);
    crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT_SET
        .store(true, std::sync::atomic::Ordering::Relaxed);
}

/// Read the last WIP checkpoint report (integration tests).
pub fn take_setback_topology_checkpoint() -> Option<String> {
    crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_LAST
        .lock()
        .unwrap_or_else(|p| p.into_inner())
        .take()
}
