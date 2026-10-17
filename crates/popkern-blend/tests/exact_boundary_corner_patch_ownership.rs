//! PK-119 semantic-ownership gate: accepted exact-boundary cap
//! commits must surface as `CornerPatchInfo` records.
//!
//! The direct-apex partition construction (production route for
//! mixed-radius N∈{3,4} apexes since the 2026-07-11 promotion)
//! commits per-sector NURBS faces. Before PK-119 those faces were
//! recorded only on `repair_report.apex_routes[..]
//! .exact_boundary_report.committed_faces`, so corner-patch
//! consumers (continuity bookkeeping, zebra-style fixtures, GUI
//! selectors) needed route-report archaeology to find the cap.
//! This gate locks the ownership contract:
//!
//! * every committed exact-boundary face has a `CornerPatchInfo`
//!   with `exact_boundary_cap = true` and `method = Coons`;
//! * the records carry the acceptance-gated continuity claim
//!   (at least G1 — the candidate gates cap/support normal
//!   coplanarity ≤ 5e-2 along every accepted seam);
//! * no stray `exact_boundary_cap` records exist for faces the
//!   route did not commit.

use popkern_blend::{
    ContinuityLevel, CornerMethod, ExactBoundaryCandidateOutcome, FilletOptions, FilletRadius,
    FilletSpec, fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// pyramid4_mixed exact-boundary construction (the PK-119
/// acceptance fixture): mixed radii promote the apex to the
/// direct-apex partition route. The serving retry is pinned OFF
/// here because this gate measures the CONSTRUCTION's ownership
/// contract; the default serve for this fixture is direct and clean
/// (see `pyramid4_mixed_default_serve_is_direct_and_clean`).
#[test]
fn pyramid4_mixed_default_options_reports_exact_boundary_corner_patches() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(4, 1.0, 0.8, tol).expect("pyramid fixture");
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.z > 0.1)
        .map(|(h, _)| h)
        .expect("pyramid apex");
    let slant_edges: Vec<_> = {
        let mut seen = std::collections::BTreeSet::new();
        pyramid
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .filter(|e| seen.insert(*e))
            .collect()
    };
    let radii = [0.04_f64, 0.06, 0.05, 0.07];
    let specs: Vec<_> = slant_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        serving_self_intersection_retry: false,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&pyramid, &specs, &options, tol, None)
        .expect("pyramid4_mixed fillet (construction instrument)");

    let committed: Vec<_> = result
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| {
            route.exact_boundary_report.candidate_status == ExactBoundaryCandidateOutcome::Accepted
        })
        .flat_map(|route| route.exact_boundary_report.committed_faces.iter().copied())
        .collect();
    assert!(
        !committed.is_empty(),
        "pyramid4_mixed under default options must commit exact-boundary sector faces \
         (mixed-radius N=4 promotion regressed?): routes={:?} warnings={:?}",
        result
            .repair_report
            .apex_routes
            .iter()
            .map(|route| (
                route.route,
                route.exact_boundary_report.candidate_status,
                route.fallback_reason.clone(),
            ))
            .collect::<Vec<_>>(),
        result
            .edge_warnings
            .iter()
            .map(|(_, w)| w)
            .collect::<Vec<_>>(),
    );

    for face in &committed {
        let patch = result
            .corner_patches
            .iter()
            .find(|p| p.face == *face)
            .unwrap_or_else(|| {
                panic!(
                    "committed exact-boundary face {face:?} has no CornerPatchInfo — \
                     consumers are back to route-report archaeology (PK-119)"
                )
            });
        assert!(
            patch.exact_boundary_cap,
            "committed face {face:?} record must carry exact-boundary ownership"
        );
        assert_eq!(
            patch.method,
            CornerMethod::Coons,
            "exact-boundary cap method must reflect the Coons-wedge construction"
        );
        assert_eq!(
            patch.vertex, apex,
            "cap record must reference the apex vertex"
        );
        assert_eq!(patch.n_sides, 4, "N=4 apex record");
        assert!(
            patch.achieved_continuity >= ContinuityLevel::G1,
            "acceptance gates G1 (normal deficit <= 5e-2); record claims {:?}",
            patch.achieved_continuity
        );
    }

    for patch in &result.corner_patches {
        if patch.exact_boundary_cap {
            assert!(
                committed.contains(&patch.face),
                "exact_boundary_cap record {:?} does not match any committed face",
                patch.face
            );
        }
    }
}

/// PK-253 default serving contract for pyramid4_mixed: the historical
/// "boundary participation" crossing was one-sided tolerance HOVER of
/// the cap's bridge boundary over the neighbour facet plane near the
/// shared raw-edge corner (window minimum at machine zero — no sign
/// change anywhere). Under the sign-change-honest fast path the corner
/// serves DIRECT with no repair machinery engaged.
#[test]
fn pyramid4_mixed_default_serve_is_direct_and_clean() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(4, 1.0, 0.8, tol).expect("pyramid fixture");
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.z > 0.1)
        .map(|(h, _)| h)
        .expect("pyramid apex");
    let slant_edges: Vec<_> = {
        let mut seen = std::collections::BTreeSet::new();
        pyramid
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .filter(|e| seen.insert(*e))
            .collect()
    };
    let radii = [0.04_f64, 0.06, 0.05, 0.07];
    let specs: Vec<_> = slant_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let result = fillet_edges(&pyramid, &specs, &FilletOptions::default(), tol, None)
        .expect("pyramid4_mixed fillet under true default options");
    // PK-597: the direct-apex cap creases across its own seams (certified
    // cap-cap normal sine 4.8e-1, body G0), so the corner takes the
    // serving exclusion re-run and the setback rung serves it (certified
    // G1, g1 2.8e-9). The only retry cause is the crease.
    let retries = &result.repair_report.serving_self_intersection_retries;
    assert!(
        !retries.is_empty()
            && retries
                .iter()
                .all(|retry| retry.detail.contains("served cap-cap crease")),
        "the crease is the only retry cause: {retries:?}"
    );
    assert!(
        result
            .corner_patches
            .iter()
            .filter(|p| p.vertex == apex)
            .all(|p| p.setback_cap),
        "the setback rung serves the corner: {:?}",
        result.corner_patches
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "default serve must validate clean: {errors:?}"
    );
}
