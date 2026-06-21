//! Stage C0 of `docs/decisions/2026-06-21-edge-curve-face-boundary-
//! coherence-invariants.md`: edge-creator attribution probe.
//!
//! Runs the headline blocker (hexagon6_variable fillet), iterates the
//! `CoherenceValidation` diagnostics emitted by
//! `corner_patch_support_completion`, and joins every L1 violation
//! (missing edge curve) against `BRep::edge_creator_locations` to
//! identify which operator authored the curve-less edge.
//!
//! Output: a `(creator file:line, count)` histogram printed to stderr,
//! sorted descending by count. The Stage C1...Cn implementer copies
//! the output into `docs/reports/2026-06-21-l1-creator-histogram.md`
//! and orders the per-creator fixes by the top entries.
//!
//! Tagged `#[ignore]` by default — this is a probe, not a gate.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionDiagnostic, BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};
use std::collections::BTreeMap;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn fillet_n_pyramid(n: usize, radii: &[f64]) -> BlendResult {
    let pyramid = BRep::n_gonal_pyramid(n, 1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyramid, n).expect("apex must exist");
    let edges = incident_edges(&pyramid, apex);
    let specs: Vec<_> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(
        &pyramid,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("pyramid fillet should not panic")
}

fn dump_l1_histogram(label: &str, result: &BlendResult) {
    let mut histogram: BTreeMap<String, usize> = BTreeMap::new();
    let mut diagnostic_index = 0usize;
    let mut total_l1 = 0usize;
    let mut untracked = 0usize;

    for d in &result.repair_report.apex_diagnostics {
        if let ApexConstructionDiagnostic::CoherenceValidation(c) = d {
            diagnostic_index += 1;
            for edge in &c.report.l1_missing_edge_curves {
                total_l1 += 1;
                let key = match result.brep.edge_creator_locations.get(edge) {
                    Some((file, line)) => format!("{file}:{line}"),
                    None => {
                        untracked += 1;
                        "<untracked>".to_string()
                    }
                };
                *histogram.entry(key).or_insert(0) += 1;
            }
            eprintln!(
                "  [{label}] diagnostic #{diagnostic_index} stage={} L1={} L2={} L3={}",
                c.stage,
                c.report.l1_missing_edge_curves.len(),
                c.report.l2_face_boundary_disagreements.len(),
                c.report.l3_missing_pcurves.len()
            );
        }
    }

    let mut sorted: Vec<(String, usize)> = histogram.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    eprintln!(
        "\n## {label} L1 creator histogram ({total_l1} L1 entries across {diagnostic_index} CoherenceValidation diagnostics; {untracked} untracked)\n"
    );
    eprintln!("| creator (file:line) | count |");
    eprintln!("|---|---|");
    for (key, count) in &sorted {
        eprintln!("| {key} | {count} |");
    }
}

#[test]
#[ignore = "Stage C0 probe — run with --ignored --nocapture and copy stderr into docs/reports/"]
fn hexagon6_variable_l1_creator_histogram() {
    let result = fillet_n_pyramid(6, &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02]);
    dump_l1_histogram("hexagon6_variable", &result);
}

#[test]
#[ignore = "Stage C2 probe — pyramid5_mixed (still L1>0 after C1)"]
fn pyramid5_mixed_l1_creator_histogram() {
    let result = fillet_n_pyramid(5, &[0.04, 0.06, 0.05, 0.07, 0.05]);
    dump_l1_histogram("pyramid5_mixed", &result);
}

#[test]
#[ignore = "Stage C2 probe — hexagon6_uniform (still L1>0 after C1)"]
fn hexagon6_uniform_l1_creator_histogram() {
    let result = fillet_n_pyramid(6, &[0.04; 6]);
    dump_l1_histogram("hexagon6_uniform", &result);
}
