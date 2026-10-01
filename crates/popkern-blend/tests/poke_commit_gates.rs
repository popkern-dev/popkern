//! PK-228/PK-229 — reachability and serving pins for the standard
//! N-quad Gregory poke commit (`apply_three_quad_decomposition_to_blend`).
//!
//! A poke-reachability probe sweep over the full workspace (2026-08-02, PK-228)
//! measured exactly which production classes still commit the
//! standard poke after PK-222 (exact-boundary), PK-225/PK-226
//! (setback fan), and the sphere/PnS rungs claimed the historical
//! cube-corner fixtures. PK-229 then closed the last default-options
//! serving class:
//!
//! * **Uniform-radius cylindrical N = 4 apexes** (pyramid apex
//!   slant-edge blends) were the ONLY class whose poke commit reached
//!   serving output under default options, at sub-commercial quality
//!   (cap_outer g0 = 2.18e-1 — the apex is sphere-owned by the
//!   syntactic pre-screen, so the setback rung is SKIPPED at N = 4
//!   (min_n = 5), yet `try_spherical_blend` never finds the exact
//!   rolling-ball sphere that exists at regular pyramid apexes: its
//!   circumcenter solve goes singular on co-circular samples of one
//!   boundary arc, the corner falls to Gregory, and the n < 5 poke
//!   acceptance ships it unmeasured). PK-229 promotes this class to
//!   the direct-apex route by default (`direct_apex_promotion.rs`
//!   uniform arm); the exact-boundary Coons partition serves at
//!   cap_outer g0 = 3.7e-8 / g1 = 4.3e-8 — pinned below at the
//!   commercial bars.
//! * **Explicit sub-patch family flags** (biquintic-full,
//!   biseptic, two-stripe, K-P binonic, K-P polar) carve the uniform
//!   arm out so the flag-gated promotion ladder keeps its reaching
//!   fixture — the poke still commits and serves there, pinned below
//!   at its measured class (PK-230 tracks migrating the flag-family
//!   dispatch pins onto this fixture).
//! * **Apex-bridge route** (direct-apex `CapKind::Gregory` plan →
//!   `assemble_gregory_cap_surface`): commits during construction on
//!   promoted apexes, but every measured serving output is claimed by
//!   the exact-boundary Coons partition (pyramid pins below) or
//!   discarded with a declining transaction (unfillet restore
//!   trials). Construction-transient, never measured serving.
//! * **Non-default flag stacks**: cube_extreme [0.02, 0.10, 0.20]
//!   commits via the vertex-blend route only when `use_kp_binonic_gregory`
//!   opts high-asymmetry trihedra out of direct-apex promotion —
//!   pinned at its measured class by kp_cap_outer_boundary_gates.
//!
//! These pins keep the serving truth honest: route provenance,
//! committed topology, and the measured cap-outer quality class. If
//! a future rung claims (or loses) either class, the asserts fire
//! and the class must be re-measured — not silently absorbed the way
//! the r8_c7 three-quad pin was (vestigial from PK-222 until PK-226
//! caught it).

use popkern_blend::types::{ApexConstructionRoute, CornerMethod};
use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, measure_cap_outer_boundary, pick_vertex_with_valence};

fn uniform_pyramid_apex_blend(
    options: &FilletOptions,
) -> (
    popkern_blend::types::BlendResult,
    popkern_core::topology::VertexHandle,
) {
    let pyramid = BRep::square_pyramid(2.0, 1.5, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyramid, 4).expect("pyramid apex");
    let edges = incident_edges(&pyramid, apex);
    assert_eq!(edges.len(), 4);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(&pyramid, &specs, options, Tolerance::default(), None)
        .expect("uniform pyramid apex blend");
    (result, apex)
}

/// PK-615 serve pin: the uniform-radius square-pyramid apex blend (4
/// slant edges only) under PURE DEFAULT options serves the setback
/// rung. PK-229 promoted the class to the direct-apex route, whose
/// exact-boundary Coons cap pinches where its four patches meet
/// (curvature growing faster than 1/distance toward the centre); the
/// setback rung serves the class certified G1 inside the fairness
/// bound.
#[test]
fn uniform_pyramid_apex_serves_the_setback_rung_by_default() {
    let (result, apex) = uniform_pyramid_apex_blend(&FilletOptions::default());
    let route = result
        .repair_report
        .apex_routes
        .iter()
        .find(|r| r.vertex == apex)
        .expect("apex route recorded");
    assert_eq!(route.route, ApexConstructionRoute::VertexBlend);
    let caps: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == apex)
        .collect();
    assert!(
        !caps.is_empty() && caps.iter().all(|p| p.setback_cap),
        "the setback rung serves: {caps:?}"
    );
    result.brep.validate().expect("valid BRep");
}

/// PK-229 serve pin, on request: the exact-boundary direct-apex route
/// (`use_direct_apex_construction`) serves 4 exact-boundary Coons
/// sector caps at commercial cap-outer quality (measured 2026-08-02:
/// cap_outer g0 = 3.671e-8, g1 = 4.270e-8).
///
/// Before PK-229 this class was the ONLY default-options serving
/// consumer of the N-quad poke, at sub-commercial quality
/// (g0 = 2.18e-1) — see the module doc for the root-cause chain.
#[test]
fn uniform_pyramid_apex_serves_direct_apex_exact_boundary_on_request() {
    let (result, apex) = uniform_pyramid_apex_blend(&FilletOptions {
        use_direct_apex_construction: true,
        ..FilletOptions::default()
    });

    let route = result
        .repair_report
        .apex_routes
        .iter()
        .find(|r| r.vertex == apex)
        .expect("apex route recorded");
    assert_eq!(
        route.route,
        ApexConstructionRoute::DirectApex,
        "the requested direct-apex route serves the uniform N=4 apex"
    );

    let caps: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == apex)
        .collect();
    assert_eq!(caps.len(), 4, "one sector cap per slant fillet");
    for patch in &caps {
        assert_eq!(
            patch.method,
            CornerMethod::Coons,
            "serving cap must be the exact-boundary Coons partition"
        );
        assert_eq!(patch.n_sides, 4);
        assert!(patch.exact_boundary_cap && !patch.setback_cap);
        let face = result.brep.faces.get(patch.face).expect("cap face");
        assert!(
            matches!(
                face.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Nurbs(_))
            ),
            "exact-boundary sector cap must carry a NURBS surface"
        );
    }
    assert_eq!(route.exact_boundary_report.committed_faces.len(), 4);
    result.brep.validate().expect("valid BRep");

    let faces: Vec<_> = caps.iter().map(|p| p.face).collect();
    let metrics = measure_cap_outer_boundary(&result.brep, &faces, 100);
    eprintln!(
        "[poke-commit] uniform pyramid apex (requested direct-apex serve): \
         cap_outer g0={:.3e} g1={:.3e} edges={}",
        metrics.worst_g0_gap, metrics.worst_g1_normal, metrics.edges_measured
    );
    assert!(metrics.edges_measured >= 4);
    assert!(
        metrics.worst_g0_gap < 5.0e-3,
        "uniform pyramid serving cap left the commercial G0 bar \
         (measured 3.671e-8): {:.3e}",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_normal < 5.0e-2,
        "uniform pyramid serving cap left the commercial G1 bar \
         (measured 4.270e-8): {:.3e}",
        metrics.worst_g1_normal
    );
}

/// PK-229 carve-out pin: an explicit sub-patch family flag
/// (two-stripe here) keeps the uniform pyramid apex on the flag-gated
/// route, where the 4-quad Gregory poke still commits and serves —
/// this is the flag-family dispatch ladder's reaching fixture
/// (PK-230). Pinned at its measured class (2026-08-02: cap_outer
/// g0 = 2.182e-1, g1 = 9.820e-1; bounds sit ~1.4x above so a silent
/// further collapse fails while the recorded class keeps serving).
/// The cap-outer displacement is surface-family-independent — it is
/// authored by the vertex-blend boundary anchors, not the patch family —
/// so the two-stripe rung's spoke-tangency win (min-cos 0.996 vs
/// 0.142) does not move this metric.
#[test]
fn uniform_pyramid_apex_subpatch_family_flag_keeps_poke_commit_reachable() {
    let options = FilletOptions {
        use_two_stripe_gregory: true,
        ..FilletOptions::default()
    };
    let (result, apex) = uniform_pyramid_apex_blend(&options);

    let route = result
        .repair_report
        .apex_routes
        .iter()
        .find(|r| r.vertex == apex)
        .expect("apex route recorded");
    assert_eq!(
        route.route,
        ApexConstructionRoute::VertexBlend,
        "explicit sub-patch family flags must keep the vertex-blend \
         dispatch (and its promotion ladder) reachable"
    );

    let caps: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == apex)
        .collect();
    assert_eq!(
        caps.len(),
        4,
        "vertex-blend dispatch must commit the 4-quad poke (one patch per sector)"
    );
    for patch in &caps {
        assert_eq!(patch.method, CornerMethod::Gregory);
        assert_eq!(patch.n_sides, 4);
        assert!(!patch.setback_cap && !patch.exact_boundary_cap);
    }
    result.brep.validate().expect("valid BRep");

    let faces: Vec<_> = caps.iter().map(|p| p.face).collect();
    let metrics = measure_cap_outer_boundary(&result.brep, &faces, 100);
    eprintln!(
        "[poke-commit] uniform pyramid apex (two-stripe carve-out, vertex-blend 4-quad): \
         cap_outer g0={:.3e} g1={:.3e} edges={}",
        metrics.worst_g0_gap, metrics.worst_g1_normal, metrics.edges_measured
    );
    assert!(metrics.edges_measured >= 4);
    assert!(
        metrics.worst_g0_gap < 3.0e-1,
        "carve-out poke commit collapsed past its measured G0 class \
         (measured 2.182e-1): {:.3e}",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_normal < 1.4,
        "carve-out poke commit collapsed past its measured G1 class \
         (measured 9.820e-1): {:.3e}",
        metrics.worst_g1_normal
    );
}

/// Mixed-radius square-pyramid apex [0.07, 0.09, 0.11, 0.13]: mixed
/// radii auto-promote to direct apex. The Gregory-cap bridge DOES
/// commit a 4-quad poke during construction (PK-228 probe measured
/// the apex-bridge accept), but the exact-boundary candidate claims
/// the SERVING output: 4 Coons sector caps, `exact_boundary_cap`,
/// committed at commercial quality (measured g0 = 4.4e-8,
/// g1 = 5.7e-6). This pin records that the apex-bridge poke commit
/// is construction-transient on served output — if it ever starts
/// serving (exact-boundary regression), the method/flag asserts
/// fire and the class must be re-measured.
///
/// The N = 3 apex-bridge commit is likewise transient-only
/// (unfillet_corner_group's declining recordless restore re-fillets
/// a restored cube and discards the result) — see the PK-228 ledger
/// entry.
#[test]
fn mixed_pyramid_apex_serves_exact_boundary_over_apex_bridge_poke() {
    let pyramid = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyramid, 4).expect("pyramid apex");
    let edges = incident_edges(&pyramid, apex);
    assert_eq!(edges.len(), 4);
    let radii = [0.07_f64, 0.09, 0.11, 0.13];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    // Dispatch-preference instrument: exact-boundary must beat the
    // apex-bridge poke at the ROUTING layer. The serving retry is
    // pinned off — with true defaults this fixture's candidate flags
    // a strip×cap self-intersection and re-routes vertex-blend (PK-245
    // default-serve contract, pinned in
    // exact_boundary_corner_patch_ownership).
    let options = FilletOptions {
        serving_self_intersection_retry: false,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&pyramid, &specs, &options, Tolerance::default(), None)
        .expect("mixed pyramid apex blend");

    let route = result
        .repair_report
        .apex_routes
        .iter()
        .find(|r| r.vertex == apex)
        .expect("apex route recorded");
    assert_eq!(
        route.route,
        ApexConstructionRoute::DirectApex,
        "mixed-radius pyramid apex must auto-promote to direct apex"
    );

    let caps: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == apex)
        .collect();
    assert_eq!(
        caps.len(),
        4,
        "apex bridge must commit the 4-quad poke (one patch per sector)"
    );
    for patch in &caps {
        assert_eq!(
            patch.method,
            CornerMethod::Coons,
            "serving cap must be the exact-boundary Coons partition, \
             not the transient apex-bridge Gregory poke"
        );
        assert_eq!(patch.n_sides, 4);
        assert!(patch.exact_boundary_cap && !patch.setback_cap);
        let face = result.brep.faces.get(patch.face).expect("cap face");
        assert!(
            matches!(
                face.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Nurbs(_))
            ),
            "exact-boundary sector cap must carry a NURBS surface"
        );
    }
    assert_eq!(
        route.exact_boundary_report.committed_faces.len(),
        4,
        "exact-boundary transaction must commit all 4 sector faces"
    );
    result.brep.validate().expect("valid BRep");

    let faces: Vec<_> = caps.iter().map(|p| p.face).collect();
    let metrics = measure_cap_outer_boundary(&result.brep, &faces, 100);
    eprintln!(
        "[poke-commit] mixed pyramid apex (serving caps): \
         cap_outer g0={:.3e} g1={:.3e} edges={}",
        metrics.worst_g0_gap, metrics.worst_g1_normal, metrics.edges_measured
    );
    assert!(metrics.edges_measured >= 4);
    // Measured 2026-08-02: g0 = 4.387e-8, g1 = 5.677e-6 — pinned at
    // the user-visible commercial bars.
    assert!(
        metrics.worst_g0_gap < 5.0e-3,
        "exact-boundary serving cap left the commercial G0 bar: {:.3e}",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_normal < 5.0e-2,
        "exact-boundary serving cap left the commercial G1 bar: {:.3e}",
        metrics.worst_g1_normal
    );
}
