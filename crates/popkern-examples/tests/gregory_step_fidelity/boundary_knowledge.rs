use popkern_core::topology::classify_brep_fidelity_blocker_ownership;

#[test]
fn gregory_surface_boundary_knowledge_has_no_remaining_blockers() {
    let construction = super::mixed_radius_gregory_cube();
    let ownership = classify_brep_fidelity_blocker_ownership(&construction, 1e-6);

    assert_eq!(
        ownership.blockers, 0,
        "exact-boundary Gregory construction should not leave remaining fidelity blockers requiring boundary-knowledge classification: {ownership:#?}"
    );
    assert!(ownership.edges.is_empty(), "{ownership:#?}");
}
