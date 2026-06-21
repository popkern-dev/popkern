//! Stage B1 coverage of
//! `docs/decisions/2026-06-21-edge-curve-face-boundary-coherence-
//! invariants.md`: corner_patch_support_completion must emit one
//! `ApexConstructionDiagnostic::CoherenceValidation` per fillet
//! operation that exercises the cap-face completion path.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionDiagnostic, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::incident_edges;

#[test]
fn corner_patch_support_completion_emits_coherence_validation_diagnostic() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let specs: Vec<_> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("uniform cube fillet should succeed");
    let coherence_count = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter(|d| matches!(d, ApexConstructionDiagnostic::CoherenceValidation(_)))
        .count();
    assert!(
        coherence_count >= 1,
        "expected at least one CoherenceValidation diagnostic, found {coherence_count} in {} apex_diagnostics",
        result.repair_report.apex_diagnostics.len()
    );
}
