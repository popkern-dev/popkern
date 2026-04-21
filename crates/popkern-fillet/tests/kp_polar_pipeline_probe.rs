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
//! # Root cause (localised 2026-04-21)
//!
//! Temporary tracing at the `blend_insert::insert_vertex_blend`
//! three-quad gate showed every other prerequisite passing, with
//! one hard failure:
//!
//! ```text
//! [insert_vertex_blend] n_sides=5 bound_some=true blend_face_he_count=5 (need 2n=10)
//! ```
//!
//! `compute_blend_boundary` succeeds and returns the expected
//! `segments.len() == 2 * n_sides` (10 alternating
//! `FilletEnd`/`FaceEdge`), but the BRep **blend face's outer
//! loop only has N half-edges**, not 2N. The
//! `face_has_2n_outer_half_edges` precondition on
//! `poke_n_gon_face_into_n_quads` rejects, so the three-quad
//! dispatch (and transitively the polar dispatch) never fires.
//!
//! Why: at N=3 (cube corner) adjacent-fillet bridges on the cube
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
//! # Phase II fix paths
//!
//! a. **Topology surgery** — after kef / kev, detect N-half-edge
//!    blend faces (N ≥ 5) and insert N synthetic midpoint
//!    vertices on each fillet arc via a mid-edge vertex Euler
//!    op. Gives 2N half-edges, lets `poke_n_gon_face_into_n_quads`
//!    run as-is. Midpoint placement should be the parameter
//!    corresponding to the arc's midpoint in the
//!    `BlendBoundary`'s `FilletEnd` curve.
//!
//! b. **Native N-gon poke** — generalise
//!    `poke_n_gon_face_into_n_quads` to handle N-half-edge outer
//!    loops. Each half-edge is a fillet-arc to be split at its
//!    midpoint parameter; the poke inserts N midpoint vertices
//!    plus N spokes to the central point, producing N quads.
//!    More invasive but avoids the synthetic-vertex coupling to
//!    specific Euler ops.
//!
//! c. **Polar-only bypass** — sidestep the three-quad scaffold
//!    entirely for N ≥ 5. Build `PolarBoundaryInput` directly
//!    from the fillet-strip surfaces + the
//!    `compute_blend_boundary` output (which is correctly 2N);
//!    construct the sub-faces via a dedicated N-fan Euler op
//!    sequence; assign the N bi-9 grids. Bypasses the 2N-HE
//!    precondition entirely.
//!
//! Same gap afflicts the binonic path on N ≥ 5 — this is a
//! pre-Phase-I upstream limitation, surfaced but not caused by
//! Phase I.
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
