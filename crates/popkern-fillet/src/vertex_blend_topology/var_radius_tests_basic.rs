#![cfg(test)]
use super::var_radius_tests_support::*;
use super::*;

#[test]
fn classify_empty_cap_list() {
    let brep = BRep::new(Tolerance::default());
    let result = classify_cap_loops(&brep, &[], &[]);
    assert!(result.is_empty());
}

/// Per_fillet_caps shorter than cap_faces — function falls back
/// to a zero centre/radius for the missing entries (no panic,
/// twins classified as NeighbourCylinder or Other).
#[test]
fn classify_handles_missing_per_fillet_caps() {
    let brep = BRep::new(Tolerance::default());
    let result = classify_cap_loops(&brep, &[], &[(Point3::new(0.0, 0.0, 0.0), 0.1)]);
    assert!(result.is_empty());
}
