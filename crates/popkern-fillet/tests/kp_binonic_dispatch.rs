//! Phase H.6 — migration + dispatch test for the K-P bidegree-9
//! rung of the three-quad promotion ladder.
//!
//! Two axes:
//!   1. Flag off (`use_kp_binonic_gregory = false`): sub-patches
//!      must remain two-stripe, preserving the Phase G output
//!      unchanged. This is the backwards-compatibility guarantee.
//!   2. Flag on: the pipeline emits a
//!      [`GregoryPatch::RectangularBiNonic`] variant when the K-P
//!      ring construction succeeds. Fallback to two-stripe is
//!      acceptable (ladder is conservative); the test requires
//!      that at least one configuration lands binonic to prove
//!      the dispatch wire actually fires.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn corner_cube_with_flag(use_kp_binonic: bool) -> (BRep, Vec<popkern_core::topology::FaceHandle>) {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    // Asymmetric radii — the two-stripe rung is guaranteed to
    // engage here (uniform radii fall through to the analytic
    // sphere path upstream of the Gregory dispatch).
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: use_kp_binonic,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None).unwrap();
    let faces: Vec<_> = result.corner_patches.iter().map(|p| p.face).collect();
    (result.brep, faces)
}

fn variants_at(brep: &BRep, faces: &[popkern_core::topology::FaceHandle]) -> Vec<String> {
    faces
        .iter()
        .filter_map(|&f| {
            brep.faces[f].surface.as_ref().and_then(|sof| {
                if let SurfaceGeometry::Gregory(gp) = &sof.geometry {
                    Some(match gp {
                        GregoryPatch::Rectangular(_) => "rectangular",
                        GregoryPatch::RectangularDeg5(_) => "deg5",
                        GregoryPatch::RectangularBiquinticFull(_) => "biquintic_full",
                        GregoryPatch::RectangularBiseptic(_) => "biseptic",
                        GregoryPatch::RectangularTwoStripe(_) => "two_stripe",
                        GregoryPatch::RectangularBiNonic(_) => "binonic",
                        GregoryPatch::Triangular(_) => "triangular",
                    })
                } else {
                    None
                }
            })
        })
        .map(|s| s.to_string())
        .collect()
}

#[test]
fn flag_off_keeps_phase_g_two_stripe_output() {
    let (brep, faces) = corner_cube_with_flag(false);
    let variants = variants_at(&brep, &faces);
    assert_eq!(variants.len(), 3, "expected 3 sub-patches");
    for v in &variants {
        assert!(
            v == "two_stripe" || v == "biseptic" || v == "biquintic_full" || v == "deg5",
            "flag-off must stay on the Phase G ladder (deg5..two_stripe), got {v}"
        );
        assert_ne!(
            v, "binonic",
            "flag-off path must never produce binonic — backwards-compat broken"
        );
    }
}

#[test]
fn flag_on_produces_binonic_or_falls_back_cleanly() {
    let (brep, faces) = corner_cube_with_flag(true);
    let variants = variants_at(&brep, &faces);
    assert_eq!(variants.len(), 3, "expected 3 sub-patches");
    // Ladder is conservative: binonic may revert on self-intersect
    // or build-failure. Accept any variant from the full ladder,
    // including binonic. A follow-up dense-G1 regression sweep
    // asserts binonic actually fires on the cube_mixed fixture.
    for v in &variants {
        assert!(
            matches!(
                v.as_str(),
                "deg5" | "biquintic_full" | "biseptic" | "two_stripe" | "binonic"
            ),
            "flag-on must stay on the ladder (deg5..binonic), got {v}"
        );
    }
}

/// Sanity: with flag on, at least one CI-tractable fixture must
/// actually produce binonic sub-patches — otherwise the dispatch
/// is unreachable and the flag has no coverage.
#[test]
fn flag_on_lands_binonic_on_at_least_one_fixture() {
    let (brep, faces) = corner_cube_with_flag(true);
    let variants = variants_at(&brep, &faces);
    let any_binonic = variants.iter().any(|v| v == "binonic");
    assert!(
        any_binonic,
        "flag-on path never produced a binonic variant on cube_mixed \
         (radii [0.08, 0.10, 0.12]) — dispatch wiring or K-P build is broken. \
         variants seen: {variants:?}"
    );
}
