use super::*;
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn enforce_g0_on_valid_cube_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let edge = brep.edges.keys().next().unwrap();
    let result = enforce_g0(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(
        result.max_displacement < f64::EPSILON,
        "valid cube should need no displacement, got {}",
        result.max_displacement
    );
    assert_eq!(result.achieved_level, ContinuityLevel::G0);
}

#[test]
fn enforce_g0_all_on_valid_cube_returns_empty() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let results = enforce_g0_all(&mut brep, &tol, 1.0);
    assert!(
        results.is_empty(),
        "valid cube should have no edges needing enforcement"
    );
}

#[test]
fn enforce_g0_on_invalid_edge_returns_error() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = EdgeHandle::default();
    let tol = Tolerance::default();
    let result = enforce_g0(&mut brep, bogus, &tol, 1.0);
    assert!(result.is_err());
}

#[test]
fn enforce_g0_no_curve_returns_zero_displacement() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let edge = brep.edges.keys().next().unwrap();
    // Remove curve geometry.
    brep.edges[edge].curve = None;
    let result = enforce_g0(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(
        result.max_displacement < f64::EPSILON,
        "edge without curve should report zero displacement"
    );
    assert_eq!(result.achieved_level, ContinuityLevel::G0);
}

#[test]
fn enforce_g0_displacement_limit_exceeded() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let edge = brep.edges.keys().next().unwrap();

    // Deliberately move a vertex to create a gap.
    let he_handle = brep.edges[edge].half_edge;
    let vh = brep.half_edges[he_handle].origin;
    let original = brep.vertices[vh].position;
    brep.vertices[vh].position = Point3::new(original.x + 0.5, original.y + 0.5, original.z + 0.5);

    // Use a very small max_displacement limit.
    let tol = Tolerance::default();
    let result = enforce_g0(&mut brep, edge, &tol, 1e-12);
    assert!(
        result.is_err(),
        "should fail when displacement exceeds limit"
    );
}

#[test]
fn enforce_g0_closes_small_gap() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let edge = brep.edges.keys().next().unwrap();

    // Deliberately move a vertex to create a small gap.
    let he_handle = brep.edges[edge].half_edge;
    let vh = brep.half_edges[he_handle].origin;
    let original = brep.vertices[vh].position;
    let offset = 0.01;
    brep.vertices[vh].position = Point3::new(original.x + offset, original.y, original.z);

    let tol = Tolerance::default();
    let result = enforce_g0(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(
        result.max_displacement > 0.0,
        "should have non-zero displacement"
    );
    assert!(
        result.max_displacement <= offset,
        "displacement should not exceed original gap"
    );
    assert_eq!(result.achieved_level, ContinuityLevel::G0);
}

// ---- G1 enforcement tests ----

#[test]
fn enforce_g1_on_cube_edge_reports_analytical() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let edge = brep.edges.keys().next().unwrap();
    let result = enforce_g1(&mut brep, edge, &tol, 1.0).unwrap();
    // Cube faces are planes (analytical) — G1 enforcement falls back to G0.
    assert_eq!(result.achieved_level, ContinuityLevel::G0);
}

#[test]
fn enforce_g1_on_invalid_edge_returns_error() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = EdgeHandle::default();
    let tol = Tolerance::default();
    let result = enforce_g1(&mut brep, bogus, &tol, 1.0);
    assert!(result.is_err());
}

// ---- G2 enforcement tests ----

#[test]
fn enforce_g2_on_cube_edge_reports_g0() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let edge = brep.edges.keys().next().unwrap();
    let result = enforce_g2(&mut brep, edge, &tol, 1.0).unwrap();
    // Cube faces are planes (analytical) — G2 enforcement can't proceed, gets G0.
    assert_eq!(result.achieved_level, ContinuityLevel::G0);
}

#[test]
fn enforce_g2_on_invalid_edge_returns_error() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = EdgeHandle::default();
    let tol = Tolerance::default();
    let result = enforce_g2(&mut brep, bogus, &tol, 1.0);
    assert!(result.is_err());
}

// ---- G3 enforcement tests ----

#[test]
fn enforce_g3_on_cube_edge_graceful_degradation() {
    // Cube faces are analytical planes — G3 enforcement cannot proceed
    // past G1 (since G1 returns G0 for analytical surfaces). Should
    // gracefully degrade.
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let edge = brep.edges.keys().next().unwrap();
    let result = super::super::g3_enforcement::enforce_g3(&mut brep, edge, &tol, 1.0).unwrap();
    // Analytical surfaces cannot achieve G1, so G3 degrades to G0.
    assert_eq!(result.achieved_level, ContinuityLevel::G0);
}

#[test]
fn enforce_g3_on_invalid_edge_returns_error() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = EdgeHandle::default();
    let tol = Tolerance::default();
    let result = super::super::g3_enforcement::enforce_g3(&mut brep, bogus, &tol, 1.0);
    assert!(result.is_err());
}

// ---- enforce_continuity_all tests ----

#[test]
fn enforce_continuity_all_g0_on_valid_cube() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let results = enforce_continuity_all(&mut brep, ContinuityLevel::G0, &tol, 1.0);
    // Valid cube needs no enforcement — all results should have zero displacement.
    for r in &results {
        assert!(
            r.max_displacement < f64::EPSILON,
            "valid cube edge {:?} should need no displacement",
            r.edge
        );
    }
}

#[test]
fn enforce_continuity_all_g1_on_valid_cube() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let results = enforce_continuity_all(&mut brep, ContinuityLevel::G1, &tol, 1.0);
    // All faces are analytical planes — G1 enforcement returns G0 for each.
    for r in &results {
        assert_eq!(r.achieved_level, ContinuityLevel::G0);
    }
}

#[test]
fn enforce_continuity_all_g2_on_valid_cube() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let results = enforce_continuity_all(&mut brep, ContinuityLevel::G2, &tol, 1.0);
    // All faces are analytical planes — G2 enforcement returns G0 for each.
    for r in &results {
        assert_eq!(r.achieved_level, ContinuityLevel::G0);
    }
}

#[test]
fn enforce_continuity_all_g3_on_valid_cube() {
    let mut brep = BRep::unit_cube(Tolerance::default()).ok().unwrap();
    let tol = Tolerance::default();
    let results = enforce_continuity_all(&mut brep, ContinuityLevel::G3, &tol, 1.0);
    // All faces are analytical planes — G3 enforcement degrades to G0 for each.
    for r in &results {
        assert_eq!(r.achieved_level, ContinuityLevel::G0);
    }
}
