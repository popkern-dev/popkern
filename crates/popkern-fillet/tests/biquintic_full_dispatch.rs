//! §5.3 item #3 extended — Karciauskas-Peters biquintic-full
//! Gregory dispatch pin.
//!
//! PK-230: migrated from the mixed-radius cube corner, where the flag
//! has had zero effect since PK-226 made mixed cube corners
//! setback-served (the old pin passed through an early-return escape
//! hatch — zero flag coverage). The measured reaching fixture is the
//! uniform square-pyramid apex (r = 0.05):
//!
//! * Flag OFF is the pure-default configuration
//!   (`use_rational_gregory`/`use_deg5_gregory` are default-on), and
//!   defaults PROMOTE uniform cylindrical N=4 apexes to the
//!   direct-apex route (PK-229) — the serve carries no rectangular
//!   Gregory sub-patches.
//! * Flag ON is an explicit legacy-family request: the uniform
//!   promotion arm stands down (`explicit_subpatch_family_requested`
//!   carve-out) and the apex serves the legacy ladder as 4
//!   `RectangularBiquinticFull` sub-patches (measured min/mean-cos
//!   0.942/0.975, validation sweep 2026-08-11).
//!
//! Mixed-radius corners cannot host this pin: the deg5-only baseline
//! IS the default flag set, and defaults promote mixed N∈{3,4}
//! apexes (PK-099) — cube_extreme reached the ladder only while the
//! pre-PK-230 grazing false positive made the serving retry re-route
//! it to legacy, which the honest validator retired.
//!
//! Acceptance rule (PK-230): this pin must FAIL, not skip, if the
//! flag stops reaching its construction on this fixture.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn incident_edges(
    brep: &BRep,
    vertex: popkern_core::topology::VertexHandle,
) -> Vec<popkern_core::topology::EdgeHandle> {
    brep.half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn fillet_pyramid4_uniform(use_biquintic_full: bool) -> BlendResult {
    let brep = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = brep
        .vertices
        .keys()
        .find(|&v| incident_edges(&brep, v).len() == 4)
        .expect("square_pyramid valence-4 apex");
    let edges = incident_edges(&brep, vertex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: use_biquintic_full,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &opts, Tolerance::default(), None)
        .expect("fillet_edges should succeed on the uniform pyramid4 apex")
}

fn count_variant(result: &BlendResult, want: fn(&GregoryPatch) -> bool) -> usize {
    result
        .corner_patches
        .iter()
        .filter_map(|p| result.brep.faces.get(p.face))
        .filter_map(|f| f.surface.as_ref())
        .filter(|sof| matches!(&sof.geometry, SurfaceGeometry::Gregory(gp) if want(gp)))
        .count()
}

fn is_biquintic_full(gp: &GregoryPatch) -> bool {
    matches!(gp, GregoryPatch::RectangularBiquinticFull(_))
}

/// The flag flips the uniform pyramid4 apex between the promoted
/// direct-apex serve (no rectangular Gregory sub-patches) and the
/// legacy 4-patch biquintic-full ladder. Locks both directions of the
/// dispatch.
#[test]
fn biquintic_full_flag_engages_pyramid4_uniform_apex() {
    let result_default = fillet_pyramid4_uniform(false);
    let result_full = fillet_pyramid4_uniform(true);

    let gregory_baseline = count_variant(&result_default, |_| true);
    assert_eq!(
        gregory_baseline, 0,
        "the default flag set serves the promoted direct-apex route on \
         pyramid4_uniform (no rectangular Gregory sub-patches) — if this changed, \
         re-measure the reaching matrix"
    );
    let biquintic_full_count = count_variant(&result_full, is_biquintic_full);
    assert_eq!(
        biquintic_full_count, 4,
        "the biquintic flag must serve the pyramid4_uniform apex as 4 \
         RectangularBiquinticFull sub-patches"
    );
}

/// The biquintic-full opt-in must not break the BRep's validation
/// invariants on the reaching fixture. Uses the cheap `validate()`
/// path (structural + geometric, no self-intersection).
#[test]
fn biquintic_full_dispatch_produces_valid_brep() {
    fillet_pyramid4_uniform(true)
        .brep
        .validate()
        .expect("biquintic-full pyramid4_uniform blend must yield a valid BRep");
}
