//! Rim transverse-G2 evidence gate (PK-519): a `g2` fillet request on
//! the pentagon flagship (base 1.0, height 1.2, r 0.16 — the widest
//! serving radius class) must MEASURE G2 through the dressup evidence
//! instrument: every gated rim seam under the 5e-2 curvature bar and
//! `achieved == G2`, no seam left unmeasured.
//!
//! This is the promoted PK-519 diagnostic probe. Baseline before the
//! dual-jet corner substrate: the ten corner-adjacent rim quarters
//! measured 0.9599 (the corner piece's shared cells could not carry
//! two strips' cross fields); landed: 4.88e-2 worst. A regression in
//! the substrate — promotion, claim routing, or the reconciles —
//! reads here as a named seam over the bar, not as a gallery
//! surprise.

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn pentagon_flagship_g2_request_measures_g2_on_every_rim_quarter() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| {
            a.1.position
                .z
                .partial_cmp(&b.1.position.z)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let pre_faces: std::collections::BTreeSet<_> = pyramid.faces.keys().collect();
    let result = fillet_edges(&pyramid, &specs, &options, tol, None).expect("pipeline");
    let brep = &result.brep;
    brep.validate().expect("valid");
    let blend_faces: Vec<_> = brep
        .faces
        .keys()
        .filter(|f| !pre_faces.contains(f))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        brep,
        &blend_faces,
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        100,
    );
    eprintln!(
        "achieved={:?} worst_g0={:.3e} worst_g1={:.3e} worst_g2={:.3e} seams={} sharp={} unmeasured={}",
        evidence.achieved,
        evidence.worst_g0_gap,
        evidence.worst_g1_normal,
        evidence.worst_g2_residual,
        evidence.seams.len(),
        evidence.sharp_seams,
        evidence.unmeasured_seams
    );
    // Per-seam listing so a regression names its seam.
    let gates = ContinuityGates::commercial();
    for s in evidence.seams.iter().filter(|s| s.gated) {
        assert!(
            s.worst_g2_residual <= gates.g2_curvature,
            "rim/fan seam {:?} over the G2 bar: {:.3e} (g0 {:.3e}, g1 {:.3e})",
            s.edge,
            s.worst_g2_residual,
            s.worst_g0_gap,
            s.worst_g1_normal
        );
    }
    assert_eq!(evidence.unmeasured_seams, 0, "every seam must measure");
    assert_eq!(
        evidence.achieved,
        ContinuityLevel::G2,
        "a g2 request on the flagship must MEASURE G2 (worst g2 {:.3e})",
        evidence.worst_g2_residual
    );
    assert!(!evidence.downgraded(), "no silent downgrade");
}
