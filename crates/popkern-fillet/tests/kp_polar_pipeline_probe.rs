//! Phase I "first lap" pipeline probe — diagnostic for the polar
//! dispatch reach.
//!
//! Runs pentagon5 fixtures through `fillet_edges` with the full rung
//! stack enabled (deg5 → binonic → polar) and reports whether the
//! three-quad dispatch reaches the polar engine. The answer as of
//! 2026-04-21 is **no** on every pyramid geometry tested: the
//! upstream vertex-blend pipeline bails to a single polar-NURBS
//! corner patch for N ≥ 5 before the three-quad path engages.
//!
//! # Finding (2026-04-21)
//!
//! For every tested geometry (pentagonal pyramid, heights 0.8–4.0 ×
//! base, radii 0.02–0.07 uniform + mixed):
//! - 1 corner patch (not 5), with `face_surface = Cylinder` — the
//!   kef-merged blend face handle resolves to a fillet-strip
//!   cylinder, not a Gregory sub-face, confirming the three-quad
//!   poke never ran.
//! - 5 `FilletOverlap` warnings (adjacent sub-fillet pairs at the
//!   apex — expected for any N ≥ 3 vertex blend, not a failure).
//! - `POPKERN_KP_POLAR_TRACE=1` prints nothing, confirming the
//!   polar promote function was never called.
//! - `CornerContinuityDowngraded` warning: "N = 5 (>= 5): polar-
//!   NURBS corner patch interior quality may be degraded; true
//!   rational Gregory for N >= 5 is pending."
//!
//! # Implication
//!
//! Phase I's bi-6 polar engine is proven on synthetic
//! `PolarBoundaryInput` (see `kp_polar_bi6_apply_tests` in
//! `popkern-core`), but in the fillet pipeline it's behind an
//! upstream N ≥ 5 gap: the blend-boundary extraction / three-quad
//! poke / self-intersection path bails out before invoking it.
//! Phase II must either (a) fix the upstream pipeline to produce
//! a valid N-gon blend face at N ≥ 5 apex vertices, or (b) wire
//! a direct polar-only bypass that constructs the blend boundary
//! from the fillet ring without going through the three-quad
//! scaffold. Same gap afflicts the binonic path on N ≥ 5 — this
//! is a pre-Phase-I upstream limitation, surfaced but not caused
//! by Phase I.
//!
//! All tests `#[ignore]` — run manually with
//! `cargo test --test kp_polar_pipeline_probe -- --ignored --nocapture`.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn opts_full_stack() -> FilletOptions {
    FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    }
}

fn probe_pyramid(n: usize, base_r: f64, height: f64, radii: &[f64], tag: &str) -> (usize, usize) {
    let pyr = BRep::n_gonal_pyramid(n, base_r, height, Tolerance::default()).unwrap();
    let Some(vertex) = pick_vertex_with_valence(&pyr, n) else {
        eprintln!("[{tag}] no N={n} apex");
        return (0, 0);
    };
    let edges = incident_edges(&pyr, vertex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    match fillet_edges(&pyr, &specs, &opts_full_stack(), Tolerance::default(), None) {
        Ok(result) => {
            let gregory_patches = result
                .corner_patches
                .iter()
                .filter(|p| {
                    result
                        .brep
                        .faces
                        .get(p.face)
                        .and_then(|f| f.surface.as_ref())
                        .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Gregory(_)))
                })
                .count();
            eprintln!(
                "[{tag}] corner_patches={} gregory_surface={} warnings={}",
                result.corner_patches.len(),
                gregory_patches,
                result.edge_warnings.len(),
            );
            (result.corner_patches.len(), gregory_patches)
        }
        Err(e) => {
            eprintln!("[{tag}] ERR {e:?}");
            (0, 0)
        }
    }
}

/// Sweep pyramid proportions + radii combinations. If any row
/// reports `gregory_surface=5`, the three-quad dispatch reached
/// the polar engine on that fixture — commercial gate now
/// exercise-able. As of 2026-04-21 no row passes.
#[ignore]
#[test]
fn probe_pentagon5_geometry_sweep() {
    let radii_uniform_005 = [0.05_f64; 5];
    let radii_mixed = [0.04_f64, 0.06, 0.05, 0.07, 0.05];
    let radii_tiny = [0.02_f64; 5];

    let mut any_success = false;
    for (base_r, height, tag) in [
        (1.0, 0.8, "tall=0.8×base — pinched"),
        (1.0, 1.5, "tall=1.5×base"),
        (1.0, 2.5, "tall=2.5×base"),
        (1.0, 4.0, "tall=4.0×base — blunt"),
        (2.0, 2.0, "big base=2, h=2"),
    ] {
        eprintln!("\n--- geom: base_r={base_r} height={height} [{tag}] ---");
        for (radii, rtag) in [
            (&radii_tiny[..], "r=0.02 uniform"),
            (&radii_uniform_005[..], "r=0.05 uniform"),
            (&radii_mixed[..], "r=mixed"),
        ] {
            let (_, gregory) = probe_pyramid(5, base_r, height, radii, rtag);
            if gregory == 5 {
                any_success = true;
            }
        }
    }
    if any_success {
        eprintln!("\n[RESULT] polar dispatch reachable on at least one fixture");
    } else {
        eprintln!(
            "\n[RESULT] polar dispatch unreachable on every fixture — \
             upstream N ≥ 5 gap (see module header)"
        );
    }
}
