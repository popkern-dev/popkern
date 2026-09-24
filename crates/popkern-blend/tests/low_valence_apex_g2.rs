//! Low-valence pyramid apexes serve the setback G2 cap (PK-518).
//!
//! Before PK-518 a G2 request on a triangular or square pyramid apex
//! (h = 1.2, r = 0.16 — constructor and boolean-carved alike; the
//! two measure identically) never reached the setback cap: the N = 4
//! KP20 cap refused the regular valence outright, and the N = 3
//! strip-conflict proxy declined every G2 strip for being slightly
//! flatter than circular. Both fell to an unpromoted N-quad cap
//! whose centroid pole sat level with its rim (strip↔cap seams at
//! 86–88°, filed as "sharp" by the evidence classifier) and the base
//! runoff edges sagged 1.0e-4 / 8.0e-4 off their strips, so every such
//! fillet read achieved G0.
//!
//! Each fixture must now serve the 4n-piece setback cap, measure
//! achieved G2 on its live seams with no tolerance-widened seam, and
//! validate fully (self-intersection included).

use popkern_blend::{
    ContinuityGates, ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
    measure_blend_continuity_evidence,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn apex_g2(n: usize) {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(n, 1.0, 1.2, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    let specs: Vec<FilletSpec> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, edge)| {
            let he = &pyramid.half_edges[edge.half_edge];
            (he.origin == apex || pyramid.half_edges[he.twin].origin == apex).then_some(eh)
        })
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    assert_eq!(specs.len(), n);
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = fillet_edges(&pyramid, &specs, &options, tol, None).expect("apex fillet serves");
    let setback: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|cp| cp.setback_cap)
        .collect();
    assert_eq!(
        setback.len(),
        4 * n,
        "n={n}: the 4n-piece setback cap serves"
    );
    assert!(
        setback
            .iter()
            .all(|cp| cp.achieved_continuity == ContinuityLevel::G2),
        "n={n}: every cap piece records achieved G2"
    );
    let blend_faces: Vec<_> = result
        .brep
        .faces
        .keys()
        .filter(|face| !pyramid.faces.contains_key(*face))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        &result.brep,
        &blend_faces,
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        33,
    );
    assert_eq!(
        evidence.achieved,
        ContinuityLevel::G2,
        "n={n}: live seams measure G2 (worst g0 {:.3e}, g1 {:.3e}, g2 {:.3e})",
        evidence.worst_g0_gap,
        evidence.worst_g1_normal,
        evidence.worst_g2_residual
    );
    assert_eq!(
        evidence.tolerance_widened_seams, 0,
        "n={n}: no seam passes only through a widened recorded tolerance"
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "n={n}: body must fully validate: {errors:?}"
    );
}

#[test]
fn triangular_pyramid_apex_serves_the_setback_g2_cap() {
    apex_g2(3);
}

#[test]
fn square_pyramid_apex_serves_the_setback_g2_cap() {
    apex_g2(4);
}

#[test]
fn pentagonal_pyramid_apex_still_serves_the_setback_g2_cap() {
    apex_g2(5);
}
