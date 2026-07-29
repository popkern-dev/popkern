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
    // Dropped 14 -> 9 with PK-125 (ribbon orientation continuity +
    // corner-compatible spokes): five pcurve-residual blockers on the
    // cap seams resolved outright.
    //
    // The remaining 9 ROOT-CAUSED (PK-132, 2026-07-18): each is a
    // CUBE-EDGE STUB between two different-radius fillet setbacks
    // whose loop wires the stub to the fillet CYLINDER face — and the
    // cube edge lies exactly (√2−1)·r off that cylinder (measured
    // 0.0331/0.0414/0.0497 = (√2−1)·{0.08,0.10,0.12}, the axis sits
    // r√2 from the edge). No pcurve can trace a curve that is not on
    // the surface, so no authoring or recompute pass can close these
    // (measured: image-aware recompute clears and rebuilds the
    // pcurves, count stays 9). The fix is CONSTRUCTION: mixed-radius
    // corner transition geometry (the stub must border a transition
    // patch or a properly trimmed neighbor, never the cylinder).
    //
    // Dropped 9 -> 0 with PK-153: the V-free DefaultOn n>=3 serving
    // IS that construction fix — the consolidated shell (winding
    // reconciliation + seam welds + rebuilt support rings) retires
    // the cube-edge stubs entirely, so no edge borders a surface it
    // does not lie on. The class is empty by construction.
    // PK-205: the class rose 0 -> 3 in the dca764a4 window (the
    // in-place V-free DefaultOn shipped glued surfaceless residual
    // membranes) and returned to 0 with the PK-205/PK-206 repair —
    // the serving-shape gate keeps membrane-free bodies, so no
    // residual-owned blockers remain either.
    assert_eq!(
        ownership.blockers, 0,
        "PK-205 census: fidelity blocker class must stay empty: {ownership:#?}"
    );
    assert_eq!(ownership.production_owned, 0, "{ownership:#?}");
    assert_eq!(ownership.construction_residual_owned, 0, "{ownership:#?}");
}
