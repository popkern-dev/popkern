//! PK-090 dispatch-candidate probe — print the KP20 polar dispatch candidate
//! records (fold counts, corner mismatches, accept/decline) for the
//! hexagon6_variable and pentagon5_mixed apex fixtures.
//!
//! Diagnostic probe only: run with
//! `cargo test -p popkern-fillet --test kp20_dispatch_candidate_probe -- --ignored --nocapture`

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionDiagnostic, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn probe(name: &str, n: usize, height: f64, radii: &[f64]) {
    let pyramid = BRep::n_gonal_pyramid(n, 1.0, height, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyramid, n).expect("apex must exist");
    let edges = incident_edges(&pyramid, apex);
    let specs: Vec<_> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&pyramid, &specs, &opts, Tolerance::default(), None)
        .expect("fillet should not panic");
    eprintln!("\n=== {name} — apex diagnostics containing \"dispatch\" ===");
    for d in &result.repair_report.apex_diagnostics {
        let ApexConstructionDiagnostic::ExactBoundaryCap(c) = d else {
            continue;
        };
        if !c.stage.contains("dispatch")
            && !c.code.contains("dispatch")
            && !c.stage.contains("boundary_solve")
            && !c.code.contains("boundary_solve")
            && !c.stage.contains("support_completion")
            && !c.code.contains("support_completion")
        {
            continue;
        }
        eprintln!(
            "  code={} stage={} detail={}",
            c.code,
            c.stage,
            c.detail_code.as_deref().unwrap_or("")
        );
    }
}

#[test]
#[ignore = "diagnostic probe only: PK-090 dispatch candidate records, run with --ignored --nocapture"]
fn pk090_dispatch_candidate_records() {
    probe(
        "hexagon6_variable",
        6,
        0.8,
        &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02],
    );
    probe("pentagon5_mixed", 5, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]);
}
