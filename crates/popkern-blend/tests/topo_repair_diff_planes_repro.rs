//! Task #136 regression gate for the different-planes chord-pair
//! repair.
//!
//! Pyramid apex-only fillets leave cylinder faces with chord pairs
//! through an off-cylinder vertex (the cap centre or unfilleted
//! base corner). The different-planes repair branch
//! (`dispatch_overrides::TOPO_REPAIR_DIFF`) once regressed N≥8 pyramids (max
//! tess-nm divergence 0.12 → 1.06); the anchor on-cylinder
//! precondition in `try_repair_chord_pair_different_planes` fixed
//! the class, and `t_136_d_no_regression_from_diff_planes` pins it.
//! The mechanism-characterisation probes that diagnosed it were
//! retired with PK-181/PK-183 (findings recorded in the fix-side
//! doc comments and the gate below).

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn fillet_apex_only(n: usize, height: f64, radius: f64) -> BRep {
    let pyr = BRep::n_gonal_pyramid(n, 1.0, height, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyr, n).expect("apex");
    let edges = incident_edges(&pyr, apex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let opts = FilletOptions::default();
    fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None)
        .unwrap()
        .brep
}

fn worst_tess_nm_divergence(brep: &BRep) -> f64 {
    brep.find_tessellation_non_manifold_edges(32, 1e-9)
        .iter()
        .map(|r| r.max_divergence)
        .fold(0.0, f64::max)
}

fn run_with_flags(n: usize, height: f64, radius: f64, repair: bool, diff: bool) -> BRep {
    let _g_repair = repair.then(|| popkern_blend::dispatch_overrides::TOPO_REPAIR.set_guard(true));
    let _g_diff = diff.then(|| popkern_blend::dispatch_overrides::TOPO_REPAIR_DIFF.set_guard(true));
    fillet_apex_only(n, height, radius)
}

/// Regression gate for the #136 chord-chain anchor-pollution
/// fix.  Asserts that enabling the diff-planes branch never
/// makes tess-nm divergence worse than the same-plane-only
/// baseline on N=4..8 pyramid apex-only fillets.
///
/// Pre-fix this would catch N=7 jumping 0.12 → 0.93 and N=8
/// 0.032 → 0.93.  Post-fix (anchor on-cyl precondition in
/// `try_repair_chord_pair_different_planes`) all values match
/// the baseline because the guard rejects every chord-chain
/// candidate.
#[test]
fn t_136_d_no_regression_from_diff_planes() {
    for n in 4..=8 {
        let baseline = run_with_flags(n, 0.8, 0.05, true, false);
        let with_diff = run_with_flags(n, 0.8, 0.05, true, true);
        let d_baseline = worst_tess_nm_divergence(&baseline);
        let d_diff = worst_tess_nm_divergence(&with_diff);
        // Allow 1% headroom for floating-point noise.
        assert!(
            d_diff <= d_baseline * 1.01 + 1e-9,
            "[N={n}] diff-planes regression: baseline={:.3e} m, diff_on={:.3e} m",
            d_baseline,
            d_diff,
        );
    }
}
