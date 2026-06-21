//! Stage B2 of `docs/decisions/2026-06-21-edge-curve-face-boundary-
//! coherence-invariants.md`: baseline generator.
//!
//! Iterates a representative KP20 polar fixture matrix and prints a
//! markdown table to stderr listing the L1, L2, L3 violation counts
//! recorded by `validate_coherence` at the end of
//! `corner_patch_support_completion`. The implementer runs this test
//! with `--ignored --nocapture` and copies the output into the
//! Stage B3 baseline doc.
//!
//! Tagged `#[ignore]` by default — full fixture iteration is slow.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionDiagnostic, BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

struct FixtureCounts {
    label: &'static str,
    coherence_diagnostics: usize,
    l1: usize,
    l2: usize,
    l3: usize,
    worst_l2_gap: f64,
}

fn fillet_cube_corner(radii: &[f64]) -> BlendResult {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let specs: Vec<_> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("cube fillet should succeed")
}

fn fillet_n_pyramid_corner(n: usize, h: f64, radii: &[f64]) -> BlendResult {
    let pyramid = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).unwrap();
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
    .expect("pyramid fillet should succeed")
}

fn count_for(label: &'static str, result: &BlendResult) -> FixtureCounts {
    let mut diagnostics = 0;
    let mut l1 = 0;
    let mut l2 = 0;
    let mut l3 = 0;
    let mut worst_l2_gap = 0.0_f64;
    for d in &result.repair_report.apex_diagnostics {
        if let ApexConstructionDiagnostic::CoherenceValidation(c) = d {
            diagnostics += 1;
            l1 += c.report.l1_missing_edge_curves.len();
            l2 += c.report.l2_face_boundary_disagreements.len();
            l3 += c.report.l3_missing_pcurves.len();
            for d in &c.report.l2_face_boundary_disagreements {
                if d.worst_gap > worst_l2_gap {
                    worst_l2_gap = d.worst_gap;
                }
            }
        }
    }
    FixtureCounts {
        label,
        coherence_diagnostics: diagnostics,
        l1,
        l2,
        l3,
        worst_l2_gap,
    }
}

fn print_markdown_table(rows: &[FixtureCounts]) {
    eprintln!("\n## KP20 polar coherence baseline\n");
    eprintln!("| fixture | CoherenceValidation diags | L1 | L2 | L3 | worst L2 gap (m) |");
    eprintln!("|---|---|---|---|---|---|");
    for row in rows {
        eprintln!(
            "| {} | {} | {} | {} | {} | {:.3e} |",
            row.label, row.coherence_diagnostics, row.l1, row.l2, row.l3, row.worst_l2_gap
        );
    }
}

#[test]
#[ignore = "Stage B2 baseline generator — run with --ignored --nocapture and copy stderr into docs/reports/"]
fn kp20_polar_coherence_baseline() {
    let mut rows = Vec::new();

    rows.push(count_for(
        "cube_uniform",
        &fillet_cube_corner(&[0.10, 0.10, 0.10]),
    ));
    rows.push(count_for(
        "cube_mixed",
        &fillet_cube_corner(&[0.08, 0.10, 0.12]),
    ));
    rows.push(count_for(
        "cube_extreme",
        &fillet_cube_corner(&[0.05, 0.10, 0.15]),
    ));
    rows.push(count_for(
        "pyramid4_mixed",
        &fillet_n_pyramid_corner(4, 0.8, &[0.04, 0.06, 0.05, 0.07]),
    ));
    rows.push(count_for(
        "pyramid5_mixed",
        &fillet_n_pyramid_corner(5, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]),
    ));
    rows.push(count_for(
        "hexagon6_uniform",
        &fillet_n_pyramid_corner(6, 0.8, &[0.04; 6]),
    ));
    rows.push(count_for(
        "hexagon6_mixed",
        &fillet_n_pyramid_corner(6, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.04, 0.06]),
    ));
    rows.push(count_for(
        "hexagon6_variable",
        &fillet_n_pyramid_corner(6, 0.8, &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02]),
    ));

    print_markdown_table(&rows);
}
