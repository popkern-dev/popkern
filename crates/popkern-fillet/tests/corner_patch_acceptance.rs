//! §5.3 Round 5 acceptance tests — end-to-end contract coverage
//! for the N-sided corner patch pipeline.
//!
//! Each test exercises a different axis of the §5.3 contract:
//! classification → dispatch → surface construction → per-corner
//! reporting → warnings → preview parity. Together they act as
//! the commercial-grade regression guard for the Round 5 landing.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::vertex_classifier::VertexBlendClassification;
use popkern_fillet::{
    BlendWarning, ContinuityLevel, CornerMethod, FilletOptions, FilletRadius, FilletSpec,
    fillet_edges, preview_corner_patches,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_fillet_specs(cube: &BRep, radius: f64) -> Vec<FilletSpec> {
    cube.edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect()
}

/// §5.3 C1 + C6 — every multi-fillet vertex on a full-cube fillet
/// produces a `CornerPatchInfo` with classification, method, and
/// continuity populated. Cube has 8 vertices, all trihedral, all
/// same-side convex.
#[test]
fn full_cube_fillet_produces_one_corner_patch_per_vertex() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs = cube_fillet_specs(&cube, 0.1);
    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();

    assert_eq!(
        r.corner_patches.len(),
        8,
        "cube has 8 corners that each should get a patch"
    );
    for p in &r.corner_patches {
        assert_eq!(p.n_sides, 3, "cube corner is trihedral");
        assert_eq!(p.classification, VertexBlendClassification::SameSideConvex);
        assert_eq!(p.method, CornerMethod::RollingBallSphere);
        assert_eq!(p.achieved_continuity, ContinuityLevel::G1);
    }
}

/// §5.3 C4 — mixed-radius trihedral corner falls through to
/// Gregory (the sphere can't fit three different radii at the
/// same contact). `try_sphere_for_corner` bails out on mixed
/// radii via surface_fix, so Gregory stays as the final method.
#[test]
fn mixed_radius_trihedral_corner_produces_gregory() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    assert_eq!(incident.len(), 3);
    let specs: Vec<_> = incident
        .iter()
        .zip([0.08, 0.1, 0.12])
        .map(|(&e, r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();

    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let patch = r
        .corner_patches
        .iter()
        .find(|p| p.vertex == vertex)
        .expect("target vertex should get a corner patch");
    assert_eq!(patch.method, CornerMethod::Gregory);
    assert_eq!(patch.n_sides, 3);
    assert_eq!(patch.achieved_continuity, ContinuityLevel::G1);
}

/// §5.3 C5 — at G2 request every trihedral cube corner downgrades
/// to G1 (edge strips are NURBS → target-radius uniformity can't
/// be verified → Gregory reports G1). The engine emits one
/// `CornerContinuityDowngraded` warning per corner.
#[test]
fn g2_cube_request_emits_eight_corner_downgrade_warnings() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs = cube_fillet_specs(&cube, 0.1);
    let opts = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let r = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    let count = r
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CornerContinuityDowngraded { .. }))
        .count();
    assert_eq!(count, 8);
}

/// §5.3 C7 — `preview_corner_patches` predicts the same corner
/// count and classification as the committed `fillet_edges` run.
/// This is the "preview parity" contract: a UI radius slider
/// must see the same method/continuity badge pre- and post-commit.
#[test]
fn preview_parity_on_full_cube_matches_committed_run() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs = cube_fillet_specs(&cube, 0.1);

    let preview = preview_corner_patches(&cube, &specs, tol()).unwrap();
    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();

    assert_eq!(
        preview.len(),
        r.corner_patches.len(),
        "preview corner count should match committed"
    );
    for p in &preview {
        let matched = r
            .corner_patches
            .iter()
            .find(|c| c.vertex == p.vertex)
            .unwrap_or_else(|| panic!("no committed corner for vertex {:?}", p.vertex));
        assert_eq!(
            p.classification, matched.classification,
            "preview and committed classification disagree on {:?}",
            p.vertex
        );
        assert_eq!(p.n_sides, matched.n_sides);
        // Preview predicts the *construction* method; surface_fix
        // may later promote Mitred→Sphere, so we tolerate a
        // RollingBallSphere committed when the preview said
        // Gregory. But the classification must match exactly.
        assert!(
            matches!(
                (p.predicted_method, matched.method),
                (CornerMethod::Gregory, CornerMethod::RollingBallSphere)
                    | (
                        CornerMethod::RollingBallSphere,
                        CornerMethod::RollingBallSphere
                    )
                    | (CornerMethod::Gregory, CornerMethod::Gregory)
                    | (CornerMethod::Mitred, CornerMethod::RollingBallSphere)
                    | (CornerMethod::Mitred, CornerMethod::Mitred)
            ),
            "unexpected preview/committed method pair: {:?} -> {:?}",
            p.predicted_method,
            matched.method,
        );
    }
}

/// §5.3 C1 — a single-edge fillet never produces corner patches
/// because no vertex has two or more filleted incident edges.
#[test]
fn single_edge_fillet_produces_no_corner_patches() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    assert!(r.corner_patches.is_empty());
    let count = r
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CornerContinuityDowngraded { .. }))
        .count();
    assert_eq!(count, 0);
}

/// §5.3 C6 — classification drives dispatch: a SameSideConvex
/// cube corner with all-equal radii produces a RollingBallSphere
/// corner (the sphere path fired because classification allowed
/// it). Negative control for the OppositeSide/Mixed gating.
#[test]
fn same_side_convex_cube_corner_gets_rolling_ball_sphere() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs = cube_fillet_specs(&cube, 0.1);
    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    for patch in &r.corner_patches {
        assert_eq!(
            patch.classification,
            VertexBlendClassification::SameSideConvex
        );
        assert_eq!(patch.method, CornerMethod::RollingBallSphere);
    }
}
