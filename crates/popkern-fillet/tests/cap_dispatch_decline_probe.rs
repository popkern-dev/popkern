//! Path-1 (post-D1.1 reality check) Step A — read the actual
//! decline-detail strings for the rational_gregory and PnS cap
//! dispatches on hexagon6_variable. The detail strings record the
//! worst cap_outer_g0 value, the offending edge / face pair, the
//! sample parameter, and the cap/support uv coordinates.
//!
//! The blocker test only asserts the decline reason contains
//! "cap_outer_g0" + "exceeds 5.000e-3" — it does not surface the
//! actual G0 magnitude. This probe prints it.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionDiagnostic, BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn fillet_hexagon6_variable() -> BlendResult {
    let pyramid = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyramid, 6).expect("apex must exist");
    let edges = incident_edges(&pyramid, apex);
    let radii = [0.02, 0.04, 0.06, 0.06, 0.04, 0.02];
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
    .expect("fillet should not panic")
}

#[test]
#[ignore = "probe — run with --ignored --nocapture"]
fn hexagon6_variable_cap_dispatch_decline_details() {
    let result = fillet_hexagon6_variable();
    eprintln!("\n## hexagon6_variable cap dispatch decline detail strings\n");

    let mut count = 0;
    for d in &result.repair_report.apex_diagnostics {
        let ApexConstructionDiagnostic::ExactBoundaryCap(c) = d else {
            continue;
        };
        if c.code != "rational_gregory.dispatch.declined" && c.code != "pns.dispatch.declined" {
            continue;
        }
        count += 1;
        eprintln!("--- decline #{count} ---");
        eprintln!("  code:   {}", c.code);
        eprintln!("  stage:  {}", c.stage);
        eprintln!("  sector: {}", c.sector_name);
        match &c.detail_code {
            Some(s) => {
                for line in s.split(';') {
                    eprintln!("    {}", line);
                }
            }
            None => eprintln!("    <no detail_code>"),
        }
    }
    if count == 0 {
        eprintln!("No rational_gregory or PnS decline diagnostics found.");
    }

    eprintln!("\n## All ExactBoundaryCap diagnostics (codes + stages)\n");
    for d in &result.repair_report.apex_diagnostics {
        if let ApexConstructionDiagnostic::ExactBoundaryCap(c) = d {
            eprintln!(
                "  code={:40} stage={:60} detail={}",
                c.code,
                c.stage,
                c.detail_code.as_deref().unwrap_or("")
            );
        }
    }

    eprintln!("\n## CapSupportContinuity diagnostics (full)\n");
    for d in &result.repair_report.apex_diagnostics {
        if let ApexConstructionDiagnostic::CapSupportContinuity(c) = d {
            eprintln!(
                "  scope={} cap_face={:?} support_face={:?} edge={:?}",
                c.measurement_scope, c.cap_face, c.support_face, c.edge
            );
            eprintln!("    continuity = {:?}", c.continuity);
            eprintln!(
                "    sample_parameter={:?} cap_uv={:?} support_uv={:?}",
                c.sample_parameter, c.cap_uv, c.support_uv
            );
            eprintln!(
                "    cap_curve_gap={:?} support_curve_gap={:?}",
                c.cap_curve_gap, c.support_curve_gap
            );
            eprintln!(
                "    worst_g1_coplanarity={:.3e} gate={:.3e} edges_measured={}",
                c.worst_g1_coplanarity, c.gate, c.edges_measured
            );
        }
    }
}
