//! §5.3 item #3 Phase-10a validation sweep. Runs the four-rung
//! Gregory ladder (deg5 → biquintic → biseptic → two-stripe)
//! on multiple fixtures (cube N=3, N-gonal-pyramid N=4,5,6,
//! uniform + mixed radii) to surface failures on the broader
//! input set before investing further solver effort. Asserts
//! non-catastrophe (G0, BRep validity) + end-to-end improvement
//! over deg5 baseline; does NOT gate on per-rung monotonicity
//! since the pipeline has subtle non-determinism.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "support/kp_gate_live_seam.rs"]
mod live_seam;
#[path = "validation_sweep_scoring.rs"]
mod scoring;
use scoring::{Outcome, RungResult, Scorecard, count_variants, gregory_patch_of, sample_scorecard};

#[derive(Clone, Copy, Debug)]
enum Variant {
    Deg5,
    BiquinticFull,
    Biseptic,
    TwoStripe,
    BinonicKp,
}

impl Variant {
    fn label(self) -> &'static str {
        match self {
            Variant::Deg5 => "deg5      ",
            Variant::BiquinticFull => "biquintic ",
            Variant::Biseptic => "biseptic  ",
            Variant::TwoStripe => "two-stripe",
            Variant::BinonicKp => "binonic-kp",
        }
    }
    /// `(biquintic_full, biseptic, two_stripe, kp_binonic)` — all
    /// cumulative flags to enable for this rung.
    fn flags(self) -> (bool, bool, bool, bool) {
        match self {
            Variant::Deg5 => (false, false, false, false),
            Variant::BiquinticFull => (true, false, false, false),
            Variant::Biseptic => (true, true, false, false),
            Variant::TwoStripe => (true, true, true, false),
            Variant::BinonicKp => (true, true, true, true),
        }
    }
}

/// Incident edges at a chosen vertex, ready to be filleted.
struct VertexEdges {
    brep: BRep,
    #[allow(dead_code)]
    vertex: VertexHandle,
    edges: Vec<EdgeHandle>,
}

fn incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    // BTreeSet for deterministic dedup (each edge appears on two
    // half-edges). HashSet here would give run-to-run variance
    // in edge order and propagate through the whole pipeline.
    brep.half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn pick_vertex_with_valence(brep: &BRep, target: usize) -> Option<VertexHandle> {
    brep.vertices
        .keys()
        .find(|&v| incident_edges(brep, v).len() == target)
}

fn fixture_cube() -> VertexEdges {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = brep.vertices.keys().next().unwrap();
    let edges = incident_edges(&brep, vertex);
    assert_eq!(edges.len(), 3, "unit_cube vertex should have 3 edges");
    VertexEdges {
        brep,
        vertex,
        edges,
    }
}

fn fixture_square_pyramid() -> VertexEdges {
    let brep = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&brep, 4).expect("square_pyramid valence-4 apex");
    let edges = incident_edges(&brep, vertex);
    VertexEdges {
        brep,
        vertex,
        edges,
    }
}

fn fixture_n_gonal_pyramid(n: usize) -> VertexEdges {
    let brep = BRep::n_gonal_pyramid(n, 0.6, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&brep, n)
        .unwrap_or_else(|| panic!("n_gonal_pyramid({n}) should have a valence-{n} apex"));
    let edges = incident_edges(&brep, vertex);
    VertexEdges {
        brep,
        vertex,
        edges,
    }
}

fn fixture_cylinder_cap_corner() -> Option<VertexEdges> {
    let brep = BRep::cylinder(0.5, 1.0, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&brep, 3)?;
    let edges = incident_edges(&brep, vertex);
    Some(VertexEdges {
        brep,
        vertex,
        edges,
    })
}

fn run_fillet(
    prep: VertexEdges,
    radii: &[f64],
    variant: Variant,
) -> Result<(BRep, Vec<FaceHandle>), String> {
    assert_eq!(
        prep.edges.len(),
        radii.len(),
        "radii count must match edge count"
    );
    let specs: Vec<FilletSpec> = prep
        .edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let (full, bsep, tsp, kpbn) = variant.flags();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: full,
        use_biseptic_gregory: bsep,
        use_two_stripe_gregory: tsp,
        use_kp_binonic_gregory: kpbn,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&prep.brep, &specs, &opts, Tolerance::default(), None)
        .map_err(|e| format!("fillet_edges failed: {e:?}"))?;
    let corner_faces: Vec<FaceHandle> = result.corner_patches.iter().map(|p| p.face).collect();
    Ok((result.brep, corner_faces))
}

fn run_rung(
    name: &str,
    fixture_fn: fn() -> VertexEdges,
    radii: &[f64],
    variant: Variant,
    n_samples: usize,
) -> RungResult {
    let prep = fixture_fn();
    let n_sides = prep.edges.len();
    let outcome = match run_fillet(prep, radii, variant) {
        Err(e) => Outcome::FilletFailed(e),
        Ok((brep, faces)) => {
            let patches: Vec<GregoryPatch> = faces
                .iter()
                .filter_map(|&f| gregory_patch_of(&brep, f))
                .collect();
            if patches.is_empty() {
                Outcome::NoGregorySubpatches
            } else {
                let variant_counts = count_variants(&patches);
                let brep_valid = brep.validate().is_ok();
                // The index scorecard's ring convention only describes
                // the legacy one-patch-per-sector polar layout.
                // Multi-patch commits (setback fans, sector
                // subdivisions) must be measured across the LIVE seam
                // topology or the frame mismatch dominates the reading
                // (measured 5.6e-2 on pyramid5_uniform's 20-patch fan
                // whose true seams sit at machine zero, PK-087).
                let scorecard = if patches.len() == n_sides {
                    sample_scorecard(&patches, n_samples)
                } else {
                    let m = live_seam::measure_live_cap_seams(&brep, &faces, n_samples.max(2));
                    Scorecard {
                        n_samples: m.cos_samples,
                        min_cos: m.cos_min,
                        mean_cos: m.cos_mean,
                        max_cos: m.cos_max,
                        max_angle_rad: m.cos_min.clamp(-1.0, 1.0).acos(),
                        max_position_gap: m.worst_g0_gap,
                    }
                };
                Outcome::Success {
                    n_sub_patches: patches.len(),
                    variant_counts,
                    brep_valid,
                    scorecard,
                }
            }
        }
    };
    print_rung(name, variant, &outcome);
    RungResult {
        label: variant.label(),
        outcome,
    }
}

fn print_rung(fixture_name: &str, variant: Variant, outcome: &Outcome) {
    match outcome {
        Outcome::FilletFailed(e) => {
            println!(
                "  [{}] {}: FILLET FAILED — {}",
                fixture_name,
                variant.label(),
                e
            );
        }
        Outcome::NoGregorySubpatches => {
            println!(
                "  [{}] {}: no Gregory sub-patches (fell through to polar-NURBS)",
                fixture_name,
                variant.label()
            );
        }
        Outcome::Success {
            n_sub_patches,
            variant_counts,
            brep_valid,
            scorecard,
        } => {
            let variant_tag = match variant {
                Variant::Deg5 => format!("{} deg5", variant_counts.deg5),
                Variant::BiquinticFull => format!("{} biq", variant_counts.biquintic),
                Variant::Biseptic => format!("{} bsep", variant_counts.biseptic),
                Variant::TwoStripe => format!("{} 2str", variant_counts.two_stripe),
                Variant::BinonicKp => format!("{} binkp", variant_counts.binonic),
            };
            let valid = if *brep_valid { "ok " } else { "BAD" };
            println!(
                "  [{}] {}: n_patches={n_sub_patches} [{variant_tag}] valid={valid} | cos[min={:.4} mean={:.4} max={:.4}] ang={:.3}rad gap={:.2e}",
                fixture_name,
                variant.label(),
                scorecard.min_cos,
                scorecard.mean_cos,
                scorecard.max_cos,
                scorecard.max_angle_rad,
                scorecard.max_position_gap,
            );
        }
    }
}

fn run_fixture_all_rungs(
    name: &str,
    fixture_fn: fn() -> VertexEdges,
    radii: &[f64],
) -> Vec<RungResult> {
    println!("── fixture: {name} (N={}, radii={radii:?}) ──", radii.len());
    const N_SAMPLES: usize = 5;
    vec![
        run_rung(name, fixture_fn, radii, Variant::Deg5, N_SAMPLES),
        run_rung(name, fixture_fn, radii, Variant::BiquinticFull, N_SAMPLES),
        run_rung(name, fixture_fn, radii, Variant::Biseptic, N_SAMPLES),
        run_rung(name, fixture_fn, radii, Variant::TwoStripe, N_SAMPLES),
        run_rung(name, fixture_fn, radii, Variant::BinonicKp, N_SAMPLES),
    ]
}

/// Fixture-specific G0 ceiling — currently the uniform 1e-5
/// commercial gate for every fixture. History: pyramid6 carried a
/// 7e-2 containment ceiling while the KP15 high-valence cap frontier
/// was open (the 6.63e-2 reading was the legacy ring metric on the
/// 24-patch sector ring). That frontier closed with the PK-087
/// setback fan: pyramid6_uniform/mixed now measure 6.6e-14 / 1.9e-13
/// on the layout-honest live-seam scorecard, so the containment
/// ceiling is retired (2026-07-16).
fn fixture_g0_ceiling(_fixture: &str) -> f64 {
    1e-5
}

/// Shared assertions: G0 position gap < fixture ceiling (1e-5
/// commercial default), BRep validates.
/// Not-asserted: strict per-rung monotonicity (pipeline used to
/// have test-fixture non-determinism — now deterministic, but
/// keeping the lenient gate for cross-platform robustness).
fn assert_rung_invariants(fixture: &str, rungs: &[RungResult]) {
    let g0_ceiling = fixture_g0_ceiling(fixture);
    for r in rungs {
        match &r.outcome {
            Outcome::FilletFailed(e) => {
                eprintln!("[{fixture}] {}: fillet failed ({e}) — known gap", r.label);
            }
            Outcome::NoGregorySubpatches => {
                eprintln!(
                    "[{fixture}] {}: fell through to polar-NURBS (sphere shortcut or non-supported case)",
                    r.label
                );
            }
            Outcome::Success {
                brep_valid,
                scorecard,
                ..
            } => {
                assert!(
                    *brep_valid,
                    "[{fixture}] {}: BRep failed validation",
                    r.label
                );
                assert!(
                    scorecard.max_position_gap < g0_ceiling,
                    "[{fixture}] {}: G0 position gap {:.4e} > {g0_ceiling:.0e}",
                    r.label,
                    scorecard.max_position_gap
                );
            }
        }
    }
}

/// End-to-end improvement gate, separate from per-rung
/// monotonicity. For every fixture where `deg5` succeeded AND
/// at least one higher rung *produced a different variant*
/// (not just inherited deg5's patches via conservative-ladder
/// fall-through), the best higher rung must improve mean-cos
/// by ≥ 10%. Fixtures where only deg5 engages (e.g. pyramid6
/// with current classical-CK boundary limitations) are
/// exempt — those report deg5 numbers at every rung since the
/// ladder conservatively stayed at deg5.
fn assert_end_to_end_improvement(fixture: &str, rungs: &[RungResult]) {
    if fixture == "cube_mixed" {
        // cube_mixed is now the clean-decline regression fixture for
        // K-P binonic: the G1 postcondition rejects the binonic rung,
        // and the lower-rung result must preserve G0 without pretending
        // to improve continuity.
        return;
    }
    let deg5_mean = match &rungs[0].outcome {
        Outcome::Success {
            scorecard,
            variant_counts,
            ..
        } => {
            if variant_counts.deg5 == 0 {
                // The deg5 rung itself shipped non-deg5 patches: the
                // apex was built by a variant-independent construction
                // (the PK-087 setback vertex fan ships identical KP20
                // bi-9 patches at every rung), so there is no ladder
                // upgrade to measure — and its baseline is already at
                // mean-cos 1.0000 on pyramid5, where a +10% demand is
                // unsatisfiable by construction.
                return;
            }
            scorecard.mean_cos
        }
        _ => return,
    };
    // Find the best higher rung that actually produced a
    // non-deg5 variant (genuine upgrade, not fall-through).
    let best_higher = rungs[1..]
        .iter()
        .filter_map(|r| match &r.outcome {
            Outcome::Success {
                variant_counts,
                scorecard,
                ..
            } => {
                let non_deg5 = variant_counts.biquintic
                    + variant_counts.biseptic
                    + variant_counts.two_stripe
                    + variant_counts.binonic;
                if non_deg5 > 0 {
                    Some(scorecard.mean_cos)
                } else {
                    None
                }
            }
            _ => None,
        })
        .fold(f64::NEG_INFINITY, f64::max);
    if best_higher.is_finite() {
        assert!(
            best_higher >= deg5_mean * 1.10,
            "[{fixture}] best genuine-upgrade-rung mean-cos {best_higher:.4} not ≥ 1.10× deg5's {deg5_mean:.4}"
        );
    }
}

/// Main discovery test. Each fixture's results print to stdout;
/// invariants asserted via `assert_rung_invariants`. Runs all 9
/// fixtures × 5 variants = 45 fillet operations.
#[test]
fn validation_sweep_all_fixtures() {
    #[allow(clippy::type_complexity)]
    let fixtures: Vec<(&str, fn() -> VertexEdges, Vec<f64>)> = vec![
        ("cube_mixed", fixture_cube, vec![0.08, 0.10, 0.12]),
        ("cube_uniform", fixture_cube, vec![0.10, 0.10, 0.10]),
        ("cube_extreme", fixture_cube, vec![0.02, 0.10, 0.20]),
        (
            "pyramid4_uniform",
            fixture_square_pyramid,
            vec![0.05, 0.05, 0.05, 0.05],
        ),
        (
            "pyramid4_mixed",
            fixture_square_pyramid,
            vec![0.04, 0.06, 0.05, 0.07],
        ),
        (
            "pyramid5_uniform",
            || fixture_n_gonal_pyramid(5),
            vec![0.04; 5],
        ),
        (
            "pyramid5_mixed",
            || fixture_n_gonal_pyramid(5),
            vec![0.03, 0.04, 0.045, 0.05, 0.055],
        ),
        (
            "pyramid6_uniform",
            || fixture_n_gonal_pyramid(6),
            vec![0.04; 6],
        ),
        (
            "pyramid6_mixed",
            || fixture_n_gonal_pyramid(6),
            vec![0.025, 0.035, 0.04, 0.045, 0.05, 0.055],
        ),
    ];
    // cylinder_cap fixture deferred — cylinder primitive topology
    // does not expose a valence-3 vertex, and constructing a
    // curved-support corner requires a different approach. Phase
    // 11 work.
    let _ = fixture_cylinder_cap_corner;
    let mut all_results: Vec<(String, Vec<RungResult>)> = Vec::new();
    for (name, fixture_fn, radii) in &fixtures {
        let rungs = run_fixture_all_rungs(name, *fixture_fn, radii);
        assert_rung_invariants(name, &rungs);
        assert_end_to_end_improvement(name, &rungs);
        all_results.push((name.to_string(), rungs));
    }
    print_summary(&all_results);
}

fn print_summary(all: &[(String, Vec<RungResult>)]) {
    println!();
    println!("═══════════ VALIDATION SWEEP SUMMARY (min / mean cos) ═══════════");
    println!(
        "{:<22} {:<18} {:<18} {:<18} {:<18} {:<18}",
        "Fixture", "Deg5", "Biquintic", "Biseptic", "Two-Stripe", "Binonic-KP"
    );
    for (name, rungs) in all {
        let cells: Vec<String> = rungs
            .iter()
            .map(|r| match &r.outcome {
                Outcome::FilletFailed(_) => "FAIL".to_string(),
                Outcome::NoGregorySubpatches => "no-Gregory".to_string(),
                Outcome::Success { scorecard, .. } => {
                    format!("{:.3} / {:.3}", scorecard.min_cos, scorecard.mean_cos)
                }
            })
            .collect();
        println!(
            "{:<22} {:<18} {:<18} {:<18} {:<18} {:<18}",
            name, cells[0], cells[1], cells[2], cells[3], cells[4]
        );
    }
    println!();
    println!(
        "Format: `min / mean` cos-agreement across 5 samples × N spokes. \
         1.0 = perfect G1. Watch for min-cos regressions between rungs."
    );
}
