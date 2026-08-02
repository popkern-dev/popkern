//! PK-228 — reachability pins for the standard N-quad Gregory poke
//! commit (`apply_three_quad_decomposition_to_blend`).
//!
//! A PK_POKE_PROBE sweep over the full workspace (2026-08-02)
//! measured exactly which production classes still commit the
//! standard poke after PK-222 (exact-boundary), PK-225/PK-226
//! (setback fan), and the sphere/PnS rungs claimed the historical
//! cube-corner fixtures:
//!
//! * **Legacy route** (`blend_insert` rational dispatch, which
//!   accepts n < 5 commits unconditionally): uniform-radius pyramid
//!   apex slant-edge blends — the ONLY class whose poke commit
//!   reaches SERVING output under default options. Uniform radii
//!   never promote to direct apex, the pyramid apex is not
//!   sphere-owned (strip axes converge above the apex), and the
//!   setback/PnS rungs decline — the poke serves, sub-commercially
//!   (PK-229).
//! * **Apex-bridge route** (direct-apex `CapKind::Gregory` plan →
//!   `assemble_gregory_cap_surface`): commits during construction on
//!   auto-promoted mixed-radius apexes, but every measured serving
//!   output is claimed by the exact-boundary Coons partition (mixed
//!   pyramid) or discarded with a declining transaction (unfillet
//!   restore trials). Construction-transient, never measured serving.
//! * **Non-default flag stacks**: cube_extreme [0.02, 0.10, 0.20]
//!   commits via the legacy route only when `use_kp_binonic_gregory`
//!   opts high-asymmetry trihedra out of direct-apex promotion —
//!   pinned at its measured class by kp_cap_outer_boundary_gates.
//!
//! These pins keep the serving truth honest: route provenance,
//! committed topology, and the measured cap-outer quality class. If
//! a future rung claims (or loses) either class, the asserts fire
//! and the class must be re-measured — not silently absorbed the way
//! the r8_c7 three-quad pin was (vestigial from PK-222 until PK-226
//! caught it).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::types::{ApexConstructionRoute, CornerMethod};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, measure_cap_outer_boundary, pick_vertex_with_valence};

/// Uniform-radius square-pyramid apex blend (4 slant edges only):
/// the legacy dispatch commits a 4-quad Gregory poke under default
/// options — the ONLY serving poke-commit class measured by the
/// PK-228 sweep. Uniform radii never auto-promote to direct apex,
/// the pyramid apex is not sphere-owned, and setback/PnS decline.
///
/// The unconditional n < 5 acceptance ships this cap UNMEASURED at
/// dispatch time; the measured class is SUB-COMMERCIAL (cap_outer
/// g0 = 2.18e-1, g1 coplanarity = 9.82e-1; validation_sweep spoke
/// min-cos 0.142 at the default rung vs 0.996 under the two-stripe
/// flag on the same fixture). The bounds below pin the measured
/// class against silent further collapse — the commercial-grade
/// serve is PK-229's scope.
#[test]
fn uniform_pyramid_apex_commits_legacy_n_quad_gregory_poke() {
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
    let result = fillet_edges(
        &pyramid,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("uniform pyramid apex blend");

    // Route provenance: uniform radii never auto-promote, so the
    // apex must be served by the legacy dispatch.
    let route = result
        .repair_report
        .apex_routes
        .iter()
        .find(|r| r.vertex == apex)
        .expect("apex route recorded");
    assert_eq!(
        route.route,
        ApexConstructionRoute::Legacy,
        "uniform pyramid apex must be served by the legacy dispatch"
    );

    let caps: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == apex)
        .collect();
    assert_eq!(
        caps.len(),
        4,
        "legacy dispatch must commit the 4-quad poke (one patch per sector)"
    );
    for patch in &caps {
        assert_eq!(patch.method, CornerMethod::Gregory);
        assert_eq!(patch.n_sides, 4);
        assert!(!patch.setback_cap && !patch.exact_boundary_cap);
        let face = result.brep.faces.get(patch.face).expect("cap face");
        assert!(
            matches!(
                face.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Gregory(_))
            ),
            "poke sub-face must carry a Gregory surface"
        );
    }
    result.brep.validate().expect("valid BRep");

    let faces: Vec<_> = caps.iter().map(|p| p.face).collect();
    let metrics = measure_cap_outer_boundary(&result.brep, &faces, 100);
    eprintln!(
        "[poke-commit] uniform pyramid apex (legacy 4-quad): cap_outer g0={:.3e} g1={:.3e} edges={}",
        metrics.worst_g0_gap, metrics.worst_g1_coplanarity, metrics.edges_measured
    );
    assert!(metrics.edges_measured >= 4);
    // Measured class 2026-08-02: g0 = 2.182e-1, g1 = 9.820e-1
    // (sub-commercial; see module doc and PK-229). Guard bounds sit
    // ~1.4x above the measured worst so a silent further collapse
    // fails while the recorded class keeps serving.
    assert!(
        metrics.worst_g0_gap < 3.0e-1,
        "uniform pyramid poke commit collapsed past its measured G0 class \
         (measured 2.182e-1): {:.3e}",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_coplanarity < 1.4,
        "uniform pyramid poke commit collapsed past its measured G1 class \
         (measured 9.820e-1): {:.3e}",
        metrics.worst_g1_coplanarity
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
    let result = fillet_edges(
        &pyramid,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
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
        metrics.worst_g0_gap, metrics.worst_g1_coplanarity, metrics.edges_measured
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
        metrics.worst_g1_coplanarity < 5.0e-2,
        "exact-boundary serving cap left the commercial G1 bar: {:.3e}",
        metrics.worst_g1_coplanarity
    );
}
