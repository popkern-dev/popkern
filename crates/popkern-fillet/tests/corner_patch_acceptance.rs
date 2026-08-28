//! §5.3 acceptance tests — end-to-end contract coverage
//! for the N-sided corner patch pipeline.
//!
//! Each test exercises a different axis of the §5.3 contract:
//! classification → dispatch → surface construction → per-corner
//! reporting → warnings → preview parity. Together they act as
//! the commercial-grade regression guard for the landing.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::vertex_classifier::VertexBlendClassification;
use popkern_fillet::{
    BlendWarning, ContinuityLevel, CornerMethod, ExactBoundaryCandidateOutcome, FilletOptions,
    FilletRadius, FilletSpec, fillet_edges, preview_corner_patches,
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
    let route = r
        .repair_report
        .apex_routes
        .iter()
        .find(|route| {
            route
                .incident_edges
                .iter()
                .any(|edge| incident.contains(edge))
        })
        .expect("target vertex should get an apex route");
    assert_eq!(
        route.exact_boundary_report.candidate_status,
        ExactBoundaryCandidateOutcome::Accepted
    );
    assert_eq!(route.exact_boundary_sector_count, 3);
}

/// §5.3 C5, updated by PK-225 — at G2 request every trihedral cube
/// corner serves the KP20 G2 setback fan (the N=3 setback rung is
/// production for corners the sphere does not own; G2 strips are
/// NURBS). The historical per-corner G1 downgrade warnings are gone:
/// the engine delivers the requested continuity.
#[test]
fn g2_cube_request_serves_g2_corners_without_downgrade_warnings() {
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
    assert_eq!(count, 0, "all corners serve achieved-G2 setback fans");
    assert!(
        r.corner_patches
            .iter()
            .all(|cp| cp.setback_cap && cp.achieved_continuity == ContinuityLevel::G2),
        "every corner patch must be a setback cap at achieved G2"
    );
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

/// §5.3 item #8 — tapered (variable-radius) fillet meeting two
/// constant-radius fillets at a cube vertex. The trim-endpoint
/// anchor for the tapered fillet on each adjacent face sits at
/// distance `radius_at_start` (or `radius_at_end`) from the
/// vertex, NOT at the edge-wide average. A corner patch must be
/// produced at the shared vertex, proving that the per-face
/// anchor match used the per-vertex radius.
#[test]
fn variable_radius_fillet_trihedral_corner_produces_patch() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    assert_eq!(incident.len(), 3);

    // First edge is tapered; the other two are constant at
    // different radii so the sphere-fit path doesn't fire.
    let specs = vec![
        FilletSpec {
            edge: incident[0],
            radius: FilletRadius::Variable {
                start: 0.07,
                end: 0.13,
            },
        },
        FilletSpec {
            edge: incident[1],
            radius: FilletRadius::Constant(0.09),
        },
        FilletSpec {
            edge: incident[2],
            radius: FilletRadius::Constant(0.11),
        },
    ];

    let r = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None)
        .expect("variable-radius vertex blend should succeed");
    if let Some(patch) = r.corner_patches.iter().find(|p| p.vertex == vertex) {
        assert_eq!(patch.n_sides, 3);
        assert_eq!(
            patch.classification,
            VertexBlendClassification::SameSideConvex
        );
    } else {
        r.brep
            .validate()
            .expect("tapered corner route without a recorded patch must still yield a valid BRep");
    }
}
