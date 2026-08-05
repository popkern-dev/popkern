//! PK-232: uniform oblique trihedral apexes serve a smooth cap
//! under default options.
//!
//! Two defects stacked to ship NO cap at all on
//! `n_gonal_pyramid(3, 1.0, 0.8)` with its slant edges filleted at
//! uniform r = 0.08:
//!
//! * `kef` removed half-edge pairs without re-pointing the endpoint
//!   vertices' `outgoing` handles. Tier-2 anchor collection starts
//!   its vertex-star walk at `outgoing` and returned false before
//!   inspecting a single face, so `split_corner_triangles` saw one
//!   anchor where two exist and the whole smooth blend fell to
//!   Mitred (`VertexModeFallback`, zero corner patches).
//! * With anchors repaired, the legacy N = 3 poke serves only
//!   g0 = 7.2e-2 (recorded declines): the corner-split chords are
//!   straight lines on the planar supports, so no boundary sphere
//!   exists for `try_spherical_blend` (fit declines honestly at
//!   max_dev 4.9e-2) and the Gregory poke is the ceiling. The
//!   direct-apex welded route serves the same apex at g0 = 4.9e-8
//!   under pure defaults, so the uniform promotion arm extends to
//!   N = 3 gated on the axis-obliqueness discriminator
//!   (mutually-orthogonal cube-class corners stay legacy for their
//!   exact rolling-ball sphere serve — pinned in
//!   `setback_n3_corner_serve.rs` / `conic_cube_corner_dispatch.rs`).

#[path = "kp_gate_helpers.rs"]
mod helpers;

use helpers::{incident_edges, measure_cap_outer_boundary};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};
use popkern_fillet::BlendWarning;
use popkern_fillet::fillet_edges;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};

fn oblique_apex_body() -> (BRep, VertexHandle, Vec<EdgeHandle>) {
    let tol = Tolerance::default();
    let brep = BRep::n_gonal_pyramid(3, 1.0, 0.8, tol).unwrap();
    let apex = brep
        .vertices
        .iter()
        .max_by(|(_, a), (_, b)| a.position.z.total_cmp(&b.position.z))
        .map(|(h, _)| h)
        .unwrap();
    let edges = incident_edges(&brep, apex);
    assert_eq!(edges.len(), 3, "trihedral apex must have 3 edges");
    (brep, apex, edges)
}

fn oblique_apex_result(options: &FilletOptions) -> popkern_fillet::BlendResult {
    let (brep, _apex, edges) = oblique_apex_body();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.08),
        })
        .collect();
    fillet_edges(&brep, &specs, options, Tolerance::default(), None)
        .expect("oblique apex fillet succeeds")
}

fn direct_apex_committed_faces(
    result: &popkern_fillet::BlendResult,
) -> Vec<popkern_core::topology::FaceHandle> {
    result
        .repair_report
        .apex_routes
        .iter()
        .flat_map(|r| r.exact_boundary_report.committed_faces.iter().copied())
        .collect()
}

fn assert_no_mitred_fallback(result: &popkern_fillet::BlendResult) {
    assert!(
        !result
            .edge_warnings
            .iter()
            .any(|(_, w)| matches!(w, BlendWarning::VertexModeFallback { .. })),
        "smooth blend must not fall back to Mitred, warnings: {:?}",
        result.edge_warnings
    );
}

/// The class serve: welded direct-apex cap at the commercial bars
/// under pure `FilletOptions::default()` (measured g0 4.9e-8 /
/// g1 4.5e-7 at promotion time).
#[test]
fn oblique_uniform_trihedral_apex_serves_direct_apex_by_default() {
    let result = oblique_apex_result(&FilletOptions::default());
    assert_no_mitred_fallback(&result);

    let committed = direct_apex_committed_faces(&result);
    assert!(
        committed.len() >= 3,
        "the promoted direct-apex route must commit the sector complex, got {committed:?}"
    );

    let metrics = measure_cap_outer_boundary(&result.brep, &committed, 100);
    assert!(
        metrics.edges_measured >= 3,
        "all cap-to-cylinder seams must be measured, got {}",
        metrics.edges_measured
    );
    assert!(
        metrics.worst_g0_gap <= 5.0e-3,
        "cap_outer_g0 {:.3e} exceeds the commercial gate 5e-3",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_coplanarity <= 5.0e-2,
        "cap_outer_g1 {:.3e} exceeds the commercial gate 5e-2",
        metrics.worst_g1_coplanarity
    );

    // PK-220: the cap's INTERNAL seams meet the same 5e-2 G1 bar on
    // the served body. The commit-time diagnostic recorded 5.885e-1
    // here — a premature-measurement artifact (mid-pipeline pcurves
    // mis-resolve the pole half of each spoke; identical control
    // points measure 3.4e-3 served). The serving-boundary refresh
    // stamps `served_mid_seam_normal`; a missing stamp is a failure.
    let mut served_pinned = false;
    for route in &result.repair_report.apex_routes {
        for diagnostic in &route.exact_boundary_report.diagnostics {
            if diagnostic.code != "exact_boundary_cap.mid_seam_normal" {
                continue;
            }
            let detail = diagnostic.detail_code.as_deref().unwrap_or_default();
            let served: f64 = detail
                .split(';')
                .find_map(|part| {
                    part.strip_prefix("served_mid_seam_normal")?
                        .strip_prefix('=')?
                        .parse()
                        .ok()
                })
                .unwrap_or(f64::NAN);
            assert!(
                served.is_finite() && served <= 5.0e-2,
                "served_mid_seam_normal {served:.3e} exceeds the crate G1 bar 5e-2 ({detail})"
            );
            served_pinned = true;
        }
    }
    assert!(
        served_pinned,
        "the served oblique cap must record a served mid-seam measurement"
    );

    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "default serve must validate clean, got {errors:?}"
    );
}

/// Explicit legacy surface-family flags stand the N = 3 uniform
/// promotion down (the legacy ladder stays measurable on its
/// reaching fixture) and the corner serves smooth Gregory patches.
/// PK-239: the sphere-ownership prediction is measurably FALSE here
/// (uniform cylinder strips on mutually oblique axes — no
/// rolling-ball sphere exists, PK-232's recorded 4.9e-2 fit
/// decline), so the corner rides the boundary-ring carve inside the
/// legacy route: the mitred chord class this pin used to ceiling at
/// 7.75e-2 (collect-all check 9, PK-235/PK-236) is replaced by an
/// exact carved ring and the serve is validate_all CLEAN.
#[test]
fn legacy_family_flag_stands_promotion_down_but_keeps_smooth_cap() {
    let options = FilletOptions {
        use_two_stripe_gregory: true,
        ..FilletOptions::default()
    };
    let result = oblique_apex_result(&options);
    assert_no_mitred_fallback(&result);
    assert!(
        direct_apex_committed_faces(&result).is_empty(),
        "explicit family flags must keep the apex on the legacy route"
    );
    assert!(
        !result.corner_patches.is_empty(),
        "the legacy route must still serve corner patches (kef anchor repair)"
    );

    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "legacy-flag ring serve must validate clean, got {errors:?}"
    );
    assert!(
        result
            .repair_report
            .apex_diagnostics
            .iter()
            .any(|d| matches!(
                d,
                popkern_fillet::ApexConstructionDiagnostic::ExactBoundaryCap(fact)
                    if fact.code == "boundary_ring.fact"
                        && fact.stage == "boundary_ring.applied"
            )),
        "the oblique legacy corner must serve through the boundary ring carve (PK-239)"
    );
}

/// The discriminator boundary: the orthogonal (cube-class) N = 3
/// uniform corner must NOT promote — its legacy carve serves the
/// exact rolling-ball sphere, which no fitted cap can beat. The
/// sphere-quality truth itself is pinned in
/// `setback_n3_corner_serve.rs` and `conic_cube_corner_dispatch.rs`;
/// this pin asserts the routing stays legacy.
#[test]
fn orthogonal_cube_corner_does_not_promote_to_direct_apex() {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let vertex = brep.vertices.keys().next().unwrap();
    let edges = incident_edges(&brep, vertex);
    assert_eq!(edges.len(), 3);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.15),
        })
        .collect();
    let result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
        .expect("cube corner fillet succeeds");
    assert!(
        direct_apex_committed_faces(&result).is_empty(),
        "orthogonal cube corner must stay on the legacy sphere-serving route"
    );
    assert!(
        result
            .corner_patches
            .iter()
            .any(|cp| cp.method == popkern_fillet::types::CornerMethod::RollingBallSphere),
        "cube corner must keep the rolling-ball sphere method"
    );
}
