use popkern_core::topology::classify_brep_fidelity_blocker_ownership;

#[test]
fn gregory_surface_boundary_knowledge_classifies_remaining_blockers() {
    let construction = super::mixed_radius_gregory_cube();
    let ownership = classify_brep_fidelity_blocker_ownership(&construction, 1e-6);

    // Census provenance: authored at 21 blockers; c3eb03ea (PK-095
    // tessellation hardening) repaired the six setback-cluster edges
    // (40v1–42v1, 48v1–50v1) — the same subset-verified 21→15 drop
    // re-pinned in cap_family_step_acceptance (PK-116). This binary
    // was outside the sweep window until --no-fail-fast (PK-118).
    // Dropped 15 -> 14 with PK-124's exact-boundary seam fidelity
    // package (canonical shared-trim curve + reprojected support
    // pcurves), matching the cap_family_step_acceptance census.
    assert_eq!(
        ownership.blockers, 14,
        "exact-boundary Gregory construction should expose the current structured fidelity blocker class: {ownership:#?}"
    );
    assert_eq!(
        ownership.production_owned, ownership.blockers,
        "remaining blockers should be production-owned, not hidden on construction residual topology: {ownership:#?}"
    );
    assert_eq!(ownership.construction_residual_owned, 0, "{ownership:#?}");
}
