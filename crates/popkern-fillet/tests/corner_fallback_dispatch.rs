//! The N=3 corner FALLBACK dispatch, pinned under the diagnostic
//! POPKERN_SETBACK_DISABLE seam.
//!
//! Since PK-226 every clean N=3 cube corner serves the setback fan
//! (uniform circular → rolling-ball sphere; everything else → KP20
//! setback cap fan), so no production cube fixture reaches the
//! fallback dispatch anymore. It remains the transactional next rung
//! when the setback surgery declines, and this binary pins what it
//! actually commits today: the direct-apex partition with three
//! Coons corner caps via the accepted exact-boundary route (PK-192:
//! N=3 stays direct-apex/Coons; the historical R8 C7 "three-quad
//! poke" contract was superseded by this route and its old pin only
//! ever passed through the accepted-route early return).
//!
//! This test owns its own binary because the env seam is
//! process-global: setting it inside the parallel lib-test binary
//! would race the setback serve pins.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_fillet::{CornerMethod, ExactBoundaryCandidateOutcome, fillet_edges};

struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        unsafe {
            std::env::set_var(key, value);
        }
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        unsafe {
            match &self.prior {
                Some(v) => std::env::set_var(self.key, v),
                None => std::env::remove_var(self.key),
            }
        }
    }
}

/// With the setback rung disabled and `use_rational_gregory` on, the
/// mixed-radius cube corner commits the direct-apex partition: three
/// Coons corner caps (n_sides = 3) through an accepted exact-boundary
/// route, Euler-consistent against the polar baseline and passing
/// structural validation.
#[test]
fn setback_disabled_corner_commits_direct_apex_coons_partition() {
    let _guard = EnvGuard::set("POPKERN_SETBACK_DISABLE", "1");

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    assert_eq!(incident_edges.len(), 3);

    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();

    // Baseline (dispatch explicitly off): single corner patch.
    let baseline_opts = FilletOptions {
        use_rational_gregory: false,
        ..FilletOptions::default()
    };
    let baseline = fillet_edges(&cube, &specs, &baseline_opts, Tolerance::default(), None).unwrap();
    let baseline_corners = baseline
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .count();
    assert_eq!(baseline_corners, 1, "baseline: one corner patch per vertex");

    // Rational-Gregory dispatch on.
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();

    let corners: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .collect();
    assert_eq!(
        corners.len(),
        3,
        "fallback dispatch must commit the three-cap direct-apex partition, got {:?}",
        corners
            .iter()
            .map(|p| (p.method, p.setback_cap))
            .collect::<Vec<_>>()
    );
    for patch in &corners {
        assert_eq!(
            patch.method,
            CornerMethod::Coons,
            "direct-apex partition caps report Coons"
        );
        assert!(!patch.setback_cap, "fallback caps are not setback caps");
        assert_eq!(patch.n_sides, 3, "n_sides drifted");
        assert!(
            result.brep.faces[patch.face].surface.is_some(),
            "fallback cap face must carry a surface"
        );
    }
    assert!(
        result.repair_report.apex_routes.iter().any(|route| {
            route.exact_boundary_report.candidate_status == ExactBoundaryCandidateOutcome::Accepted
                && route.exact_boundary_report.committed_faces.len() == 3
        }),
        "the exact-boundary route must be Accepted with 3 committed faces"
    );

    // Euler consistency vs the polar baseline (the partition rewires
    // apex topology, so raw counts differ; χ must not).
    let dv = result.brep.vertices.len() as i64 - baseline.brep.vertices.len() as i64;
    let de = result.brep.edges.len() as i64 - baseline.brep.edges.len() as i64;
    let df = result.brep.faces.len() as i64 - baseline.brep.faces.len() as i64;
    assert_eq!(
        dv - de + df,
        0,
        "fallback partition must be Euler-consistent vs baseline (dV={dv} dE={de} dF={df})"
    );

    result
        .brep
        .validate()
        .expect("fallback dispatch output must pass BRep validate()");
}
