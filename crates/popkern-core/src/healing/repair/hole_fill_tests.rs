use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_has_no_holes() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let holes = detect_holes(&brep);
    assert!(holes.is_empty(), "valid cube should have no boundary holes");
}

#[test]
fn valid_cube_has_no_free_edges() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let free = find_free_edges(&brep);
    assert!(
        free.is_empty(),
        "valid cube should have no free boundary edges"
    );
}

#[test]
fn fill_holes_on_valid_cube_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let result = fill_holes(&mut brep, HoleFillStrategy::Extension, &tol);
    assert_eq!(result.holes_detected, 0);
    assert_eq!(result.holes_filled, 0);
    assert_eq!(result.strategy, HoleFillStrategy::Extension);
}

#[test]
fn empty_brep_has_no_holes() {
    let brep = BRep::new(Tolerance::default());
    let holes = detect_holes(&brep);
    assert!(holes.is_empty());
}

#[test]
fn fill_by_extension_on_valid_cube_returns_none() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    // A valid cube has no holes, so detect_holes returns empty.
    let holes = detect_holes(&brep);
    assert!(holes.is_empty());
    // Manually create a fake hole and verify extension returns None
    // when there are no actual boundary edges.
    let fake_hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 0,
    };
    let result = fill_by_extension(&brep, &fake_hole, &tol);
    assert!(
        result.is_none(),
        "extension should return None for empty hole boundary"
    );
}

#[test]
fn auto_select_small_hole_prefers_extension() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    // 3-sided hole: should prefer extension.
    let hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 3,
    };
    assert_eq!(
        auto_select_strategy(&brep, &hole, &tol),
        HoleFillStrategy::Extension
    );
}

#[test]
fn auto_select_large_hole_prefers_energy_minimizing() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    // 5-sided hole: should prefer energy-minimizing.
    let hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 5,
    };
    assert_eq!(
        auto_select_strategy(&brep, &hole, &tol),
        HoleFillStrategy::EnergyMinimizing
    );
}

#[test]
fn auto_select_4sided_hole_without_extension() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    // 4-sided hole with no real boundary edges: extension will fail,
    // should fall back to energy-minimizing.
    let hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 4,
    };
    assert_eq!(
        auto_select_strategy(&brep, &hole, &tol),
        HoleFillStrategy::EnergyMinimizing
    );
}

#[test]
fn energy_minimizing_reports_detection_only() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let result = fill_holes(&mut brep, HoleFillStrategy::EnergyMinimizing, &tol);
    assert_eq!(result.holes_detected, 0);
    assert_eq!(result.holes_filled, 0);
    assert_eq!(result.strategy, HoleFillStrategy::EnergyMinimizing);
}
