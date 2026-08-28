//! PK-237 boundary-ring decline-class fixtures: probes that reach
//! the ring carve's constructor classes beyond the analytic
//! cylinder saddle (swept NURBS strips, generalized rail-spur
//! recut). Split from `saddle_dispatch_tests.rs` for the structural
//! file budget; shares its fixture builders via `super`.

use super::*;

/// PK-237: variable-radius fillets on the saddle edges. Straight
/// edges with linearly varying radius emit swept NURBS strips,
/// reaching the ring carve's non-analytic constructor classes: the
/// pushed station arcs come from the exact interior iso, and the
/// mitre-chord junction (measured 8.1e-2 off the carved boundary on
/// this fixture) goes through the generalized rail-spur recut whose
/// rail is the strip chart's exact border iso curve. The serve must
/// be validate_all CLEAN with the ring fact recorded — before the
/// recut generalization this fixture declined with "strip surface
/// has no rail direction" and served the mitred fallback.
#[test]
fn variable_radius_saddle_serves_boundary_ring() {
    let (brep, vertex) = build_notched_cube()
        .expect("notched-cube boolean fixture must build (PK-231: no more silent skips)");
    let incident_edges: Vec<_> = {
        let mut seen = std::collections::BTreeSet::new();
        brep.half_edges
            .values()
            .filter(|he| he.origin == vertex)
            .map(|he| he.edge)
            .filter(|e| seen.insert(*e))
            .collect()
    };
    assert_eq!(incident_edges.len(), 3);
    let radii = [(0.05_f64, 0.09_f64), (0.07, 0.11), (0.08, 0.06)];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &(start, end))| FilletSpec {
            edge: e,
            radius: FilletRadius::Variable { start, end },
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let result = crate::fillet_edges(&brep, &specs, &opts, Tolerance::default(), None)
        .expect("variable-radius saddle fillet pipeline must complete");
    let facts: Vec<String> = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(f)
                if f.code == "boundary_ring.fact" =>
            {
                Some(f.stage.clone())
            }
            _ => None,
        })
        .collect();
    let mut structural = 0usize;
    let mut worst = 0.0f64;
    let mut cs = 0usize;
    for e in result.brep.validate_all() {
        match e {
            popkern_core::error::ValidationError::CurveSurfaceDeviation {
                max_deviation, ..
            } => {
                cs += 1;
                worst = worst.max(max_deviation);
            }
            _ => structural += 1,
        }
    }
    eprintln!(
        "PK237 variable saddle: facts={facts:?} structural={structural} check9={cs} worst={worst:.3e}"
    );
    assert!(
        facts.iter().any(|f| f == "boundary_ring.applied"),
        "variable-radius saddle must serve through the ring carve; facts: {facts:?}"
    );
    assert_eq!(structural, 0, "ring serve must be structurally clean");
    assert_eq!(
        cs, 0,
        "ring serve must be check-9 clean (worst {worst:.3e})"
    );
}
