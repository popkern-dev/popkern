//! Compare routing options for hexagon6_variable: default (PnS →
//! rational_gregory → base_patch fallback) vs `use_direct_apex_construction=true`
//! (direct-apex first). Reports the cap_outer_g0 / G1 achieved on
//! each path so the routing decision is grounded in measurement.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionDiagnostic, BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn fillet_hexagon6_variable_with(options: FilletOptions) -> BlendResult {
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
    fillet_edges(&pyramid, &specs, &options, Tolerance::default(), None)
        .expect("hexagon6_variable fillet should not panic")
}

fn report_route_outcome(label: &str, result: &BlendResult) {
    eprintln!("\n## {label}\n");
    let mut declined_routes: Vec<(&str, &str)> = Vec::new();
    let mut continuity_records: Vec<(String, f64, f64)> = Vec::new();
    for d in &result.repair_report.apex_diagnostics {
        match d {
            ApexConstructionDiagnostic::ExactBoundaryCap(c) => {
                if c.code.ends_with(".dispatch.declined")
                    || c.code.ends_with(".dispatch.accepted")
                    || c.code.contains(".dispatch.")
                {
                    let detail_short = c
                        .detail_code
                        .as_deref()
                        .and_then(|d| d.split(';').find(|s| s.starts_with("cap_outer_g0 ")))
                        .unwrap_or("");
                    declined_routes.push((c.code.as_str(), detail_short));
                }
            }
            ApexConstructionDiagnostic::CapSupportContinuity(c) => {
                let g0 = c.continuity.worst_g0_gap.unwrap_or(f64::NAN);
                let g1 = c.continuity.worst_g1_normal_deficit.unwrap_or(f64::NAN);
                continuity_records.push((c.measurement_scope.clone(), g0, g1));
            }
            _ => {}
        }
    }
    eprintln!("Dispatch routes recorded:");
    for (code, detail) in &declined_routes {
        if !detail.is_empty() {
            eprintln!("  {code}: {detail}");
        } else {
            eprintln!("  {code}");
        }
    }
    eprintln!("Cap/support continuity scopes:");
    for (scope, g0, g1) in &continuity_records {
        eprintln!("  scope={scope:32} worst_g0={g0:.3e} m  worst_g1={g1:.3e}");
    }
    eprintln!(
        "Apex route summary: {}",
        result
            .repair_report
            .apex_routes
            .iter()
            .map(|r| format!("{:?}", r.route))
            .collect::<Vec<_>>()
            .join(", ")
    );
    eprintln!("Apex route fallback reasons:");
    for r in &result.repair_report.apex_routes {
        if let Some(reason) = &r.fallback_reason {
            eprintln!("  {:?}: {reason}", r.route);
        }
    }
    eprintln!("All ExactBoundaryCap diagnostic codes/stages:");
    let mut seen_codes = std::collections::BTreeSet::new();
    for d in &result.repair_report.apex_diagnostics {
        if let ApexConstructionDiagnostic::ExactBoundaryCap(c) = d
            && seen_codes.insert((c.code.clone(), c.stage.clone()))
        {
            let detail_first_line = c
                .detail_code
                .as_deref()
                .map(|d| {
                    d.lines()
                        .next()
                        .unwrap_or("")
                        .chars()
                        .take(140)
                        .collect::<String>()
                })
                .unwrap_or_default();
            eprintln!(
                "  code={:50} stage={:50} {}",
                c.code, c.stage, detail_first_line
            );
        }
    }
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_default_routing() {
    let result = fillet_hexagon6_variable_with(FilletOptions::default());
    report_route_outcome("hexagon6_variable / default options", &result);
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_direct_apex_routing() {
    let options = FilletOptions {
        use_direct_apex_construction: true,
        ..FilletOptions::default()
    };
    let result = fillet_hexagon6_variable_with(options);
    report_route_outcome(
        "hexagon6_variable / use_direct_apex_construction=true",
        &result,
    );
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_kp_polar_gregory_routing() {
    let options = FilletOptions {
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_hexagon6_variable_with(options);
    report_route_outcome("hexagon6_variable / use_kp_polar_gregory=true", &result);
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_kp_binonic_routing() {
    let options = FilletOptions {
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_hexagon6_variable_with(options);
    report_route_outcome("hexagon6_variable / use_kp_binonic_gregory=true", &result);
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_two_stripe_routing() {
    let options = FilletOptions {
        use_two_stripe_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_hexagon6_variable_with(options);
    report_route_outcome("hexagon6_variable / use_two_stripe_gregory=true", &result);
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_biseptic_routing() {
    let options = FilletOptions {
        use_biseptic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_hexagon6_variable_with(options);
    report_route_outcome("hexagon6_variable / use_biseptic_gregory=true", &result);
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_biquintic_full_routing() {
    let options = FilletOptions {
        use_biquintic_full_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_hexagon6_variable_with(options);
    report_route_outcome(
        "hexagon6_variable / use_biquintic_full_gregory=true",
        &result,
    );
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_pns_disabled_routing() {
    // SAFETY: tests run sequentially under PROFILE_LOCK in the alloc
    // profile, but this test stands alone; setting an env var only
    // affects this process. Restoring at end keeps the global state
    // tidy if other tests in this binary run after.
    let prev = std::env::var("POPKERN_PNS").ok();
    // SAFETY: single-threaded test; env writes are visible.
    unsafe {
        std::env::set_var("POPKERN_PNS", "0");
    }
    let result = fillet_hexagon6_variable_with(FilletOptions::default());
    report_route_outcome("hexagon6_variable / POPKERN_PNS=0 (PnS disabled)", &result);
    // SAFETY: ditto; restoring previous value.
    unsafe {
        if let Some(v) = prev {
            std::env::set_var("POPKERN_PNS", v);
        } else {
            std::env::remove_var("POPKERN_PNS");
        }
    }
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn hexagon6_variable_pns_disabled_kp_polar_routing() {
    let prev = std::env::var("POPKERN_PNS").ok();
    unsafe {
        std::env::set_var("POPKERN_PNS", "0");
    }
    let options = FilletOptions {
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_hexagon6_variable_with(options);
    report_route_outcome(
        "hexagon6_variable / POPKERN_PNS=0 + use_kp_polar_gregory=true",
        &result,
    );
    unsafe {
        if let Some(v) = prev {
            std::env::set_var("POPKERN_PNS", v);
        } else {
            std::env::remove_var("POPKERN_PNS");
        }
    }
}
