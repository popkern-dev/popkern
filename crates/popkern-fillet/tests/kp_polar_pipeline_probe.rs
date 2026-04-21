//! Phase I "first lap" pipeline probe — diagnostic for the polar
//! dispatch reach.
//!
//! Runs pentagon5 fixtures through `fillet_edges` with the full
//! rung stack enabled (deg5 → binonic → polar) and reports
//! whether the three-quad dispatch reaches the polar engine.
//!
//! # Status (as of the N-HE poke fix)
//!
//! The upstream topology gap identified at 2026-04-21 — the
//! blend face having N half-edges instead of 2N at sharp N ≥ 5
//! apex geometries — is closed. `poke_n_gon_face_into_n_quads`
//! in `crates/popkern-fillet/src/construct/poke_hex.rs` now
//! accepts either shape for N ≥ 5: on N-HE input it splits each
//! arc at `t=0.5` via `BRep::split_edge` to synthesize the
//! missing bridge-endpoint vertices, then runs the classical
//! `mev + (N-1) × mef` algorithm unchanged.
//!
//! On some pyramid fixtures the three-quad dispatch now reaches
//! the polar promote function (confirmed by
//! `POPKERN_KP_POLAR_TRACE=1` traces firing — e.g.,
//! `[kp_polar_bi6 N=5] fallback @ non-two-stripe input`). The
//! *next* gap is downstream: the promotion ladder
//! (biquintic_full → biseptic → two_stripe) doesn't fully prep
//! N=5 sub-patches into the two-stripe form that polar consumes,
//! so polar still rejects at the "non-two-stripe input" check.
//! That's a separate phase of work.
//!
//! # Historical root cause (documented for future reference)
//!
//! Temporary tracing at the `blend_insert::insert_vertex_blend`
//! three-quad gate showed every other prerequisite passing, with
//! one hard failure:
//!
//! ```text
//! [insert_vertex_blend] n_sides=5 bound_some=true blend_face_he_count=5 (need 2n=10)
//! ```
//!
//! At N=3 (cube corner) adjacent-fillet bridges on the cube
//! faces are separate BRep edges that survive the kef / kev Euler
//! ops — the merged blend face has 2N = 6 half-edges (3 arcs +
//! 3 bridges). At N=5 pyramidal apex the slant-face geometry is
//! sharper; the fillet trim-curves converge tightly enough that
//! the "bridge" region is absorbed into the fillet-arc edges
//! during kef merging — the merged blend face has only N = 5
//! half-edges (all fillet arcs, no topologically-distinct
//! bridges). The abstract `BlendBoundary` still carries 10
//! segments (5 arcs + 5 bridge curves on the supporting slant
//! faces) because those are geometric curves, not topological
//! edges — but the BRep face loop has lost the distinction.
//!
//! Same gap afflicted the binonic path on N ≥ 5 — one poke
//! upgrade fixes both rungs.
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
    let half_angle_deg = (base_r / height).atan().to_degrees();
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
                "[{tag} apex={half_angle_deg:.0}°] corner_patches={} gregory_surface={} warnings={}",
                result.corner_patches.len(),
                gregory_patches,
                result.edge_warnings.len(),
            );
            (result.corner_patches.len(), gregory_patches)
        }
        Err(e) => {
            eprintln!("[{tag} apex={half_angle_deg:.0}°] ERR {e:?}");
            (0, 0)
        }
    }
}

/// Sweep pyramid proportions + radii combinations. Reports the
/// state of the dispatch chain at each row:
/// - `corner_patches=5 gregory_surface=5`: three-quad dispatch
///   reached polar AND all sub-faces committed as Gregory.
///   Commercial gate can be flipped on for that fixture.
/// - `corner_patches=5 gregory_surface=N<5`: poke unlocked the
///   three-quad topology but downstream promotion / surface
///   assignment left some sub-faces non-Gregory. Polar may have
///   fired or not; check stderr for `[kp_polar_bi6 N=5]` traces.
/// - `corner_patches=1`: pipeline fell through to single-face
///   polar-NURBS. Upstream gap still blocks the dispatch for
///   that fixture.
#[ignore]
#[test]
fn probe_pentagon5_geometry_sweep() {
    let radii_uniform_005 = [0.05_f64; 5];
    let radii_mixed = [0.04_f64, 0.06, 0.05, 0.07, 0.05];
    let radii_tiny = [0.02_f64; 5];

    let mut fully_passing = false;
    let mut any_five_subfaces = false;
    // Full sweep: from pathologically-sharp (h=0.5, apex half-angle
    // 63°) through commercial-grade CAD (h=3.0, 18°) to near-flat
    // (h=10.0, 6°). Apex half-angle is what actually stresses the
    // fillet pipeline — sharper = more fillet overlap near C,
    // blunter = gentler cross-derivative magnitudes.
    for (base_r, height, tag) in [
        (1.0, 0.5, "pinched"),
        (1.0, 0.8, "tall=0.8×base"),
        (1.0, 1.5, "moderate"),
        (1.0, 2.0, "moderate-blunt"),
        (1.0, 3.0, "commercial-gentle"),
        (1.0, 5.0, "blunt"),
        (1.0, 10.0, "very blunt"),
    ] {
        eprintln!("\n--- geom: base_r={base_r} height={height} [{tag}] ---");
        for (radii, rtag) in [
            (&radii_tiny[..], "r=0.02 uniform"),
            (&radii_uniform_005[..], "r=0.05 uniform"),
            (&radii_mixed[..], "r=mixed"),
        ] {
            let (corner_patches, gregory) = probe_pyramid(5, base_r, height, radii, rtag);
            if corner_patches == 5 {
                any_five_subfaces = true;
            }
            if gregory == 5 {
                fully_passing = true;
            }
        }
    }
    if fully_passing {
        eprintln!(
            "\n[RESULT] polar dispatch fully reaches bi-9 Gregory output on at least one fixture"
        );
    } else if any_five_subfaces {
        eprintln!(
            "\n[RESULT] poke unblocked three-quad topology (5 sub-faces produced on some fixtures); \
             downstream promotion ladder and surface-assignment still gap at N=5"
        );
    } else {
        eprintln!("\n[RESULT] polar dispatch unreachable on every fixture");
    }
}
