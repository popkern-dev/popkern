use popkern_core::topology::classify_brep_fidelity_blocker_ownership;

#[test]
fn gregory_surface_boundary_knowledge_classifies_remaining_blockers() {
    let construction = super::mixed_radius_gregory_cube();
    let ownership = classify_brep_fidelity_blocker_ownership(&construction, 1e-6);

    assert_eq!(
        ownership.blockers, 18,
        "exact-boundary Gregory construction should expose the current structured fidelity blocker class: {ownership:#?}"
    );
    assert_eq!(
        ownership.production_owned, ownership.blockers,
        "remaining blockers should be production-owned, not hidden on construction residual topology: {ownership:#?}"
    );
    assert_eq!(ownership.construction_residual_owned, 0, "{ownership:#?}");
}
