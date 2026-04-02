use super::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn thicken_cube_produces_result() {
    // A cube is a closed solid, not a sheet, but thicken should still
    // work (producing outer + inner skins with no side walls).
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    // Outer (6) + inner (6) = 12 faces.
    assert_eq!(result.brep.faces.len(), 12);
}

#[test]
fn thicken_empty_brep_errors() {
    let brep = BRep::new(tol());
    let result = thicken(&brep, 0.1, &ShellOptions::default(), tol());
    assert!(matches!(result, Err(ShellError::NoFaces)));
}

// --- New tests for zero/extreme thickness values ---

#[test]
fn thicken_cube_very_small_thickness() {
    // A very small thickness (1e-6) should still produce a valid result
    // with the expected topology: outer (6) + inner (6) = 12 faces.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 1e-6, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 12);
}

#[test]
fn thicken_cube_zero_thickness_produces_degenerate_or_error() {
    // Zero thickness means half = 0.0, so offset distance is zero.
    // This should either produce a degenerate result (outer and inner
    // coincide) or succeed trivially. Either way, no panic.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 0.0, &ShellOptions::default(), tol());
    // Accept either success (degenerate but topologically consistent)
    // or a well-typed error.
    match result {
        Ok(r) => {
            // If it succeeds, outer + inner skins should be present.
            assert!(r.brep.faces.len() >= 6);
        }
        Err(e) => {
            // Any ShellError variant is acceptable; it should not panic.
            let _ = format!("{e}");
        }
    }
}

#[test]
fn thicken_cube_negative_thickness() {
    // Negative thickness reverses the offset direction.
    // half = -0.1, so the "outward" skin goes inward and vice versa.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, -0.2, &ShellOptions::default(), tol()).unwrap();
    // Topology should be the same as positive thickness: 12 faces.
    assert_eq!(result.brep.faces.len(), 12);
}

#[test]
fn thicken_cube_negative_vs_positive_symmetry() {
    // Thickening with positive and negative thickness should produce
    // topologically equivalent results (same number of faces, edges,
    // vertices).
    let brep = BRep::unit_cube(tol()).unwrap();
    let pos = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    let neg = thicken(&brep, -0.2, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(pos.brep.faces.len(), neg.brep.faces.len());
    assert_eq!(pos.brep.edges.len(), neg.brep.edges.len());
    assert_eq!(pos.brep.vertices.len(), neg.brep.vertices.len());
}

#[test]
fn thicken_cube_large_thickness() {
    // A very large thickness should still succeed for a cube since
    // planar faces can be offset arbitrarily far.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 100.0, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 12);
}

#[test]
fn thicken_cube_topology_counts() {
    // Verify detailed topology counts for a thickened cube.
    // Input cube: V=8, E=12, F=6, Euler = 8 - 12 + 6 = 2.
    // Thickened (outer + inner, no side walls for closed solid):
    //   F = 12 (6 outer + 6 inner)
    //   V = 16 (8 outer + 8 inner)
    //   E = 24 (12 outer + 12 inner)
    //   Euler per shell: V - E + F = 16 - 24 + 12 = 4 (two disjoint shells)
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    let f = result.brep.faces.len();
    let e = result.brep.edges.len();
    let v = result.brep.vertices.len();
    assert_eq!(f, 12, "expected 12 faces (6 outer + 6 inner)");
    assert_eq!(v, 16, "expected 16 vertices (8 outer + 8 inner)");
    assert_eq!(e, 24, "expected 24 edges (12 outer + 12 inner)");
    // Euler characteristic for two disjoint closed shells.
    let euler = (v as i64) - (e as i64) + (f as i64);
    assert_eq!(euler, 4, "expected Euler characteristic 4 for two shells");
}

#[test]
fn thicken_cylinder_produces_result() {
    // Thicken a cylinder: 3 faces (lateral + 2 caps).
    // Result should have 6 faces (3 outer + 3 inner).
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn thicken_cylinder_topology_counts() {
    // Input cylinder: V=2, E=3, F=3, Euler = 2 - 3 + 3 = 2.
    // Thickened (outer + inner):
    //   F = 6 (3 outer + 3 inner)
    //   V = 4 (2 outer + 2 inner)
    //   E = 6 (3 outer + 3 inner)
    //   Euler = 4 - 6 + 6 = 4 (two disjoint shells)
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    let f = result.brep.faces.len();
    let e = result.brep.edges.len();
    let v = result.brep.vertices.len();
    assert_eq!(f, 6, "expected 6 faces (3 outer + 3 inner)");
    assert_eq!(v, 4, "expected 4 vertices (2 outer + 2 inner)");
    assert_eq!(e, 6, "expected 6 edges (3 outer + 3 inner)");
    let euler = (v as i64) - (e as i64) + (f as i64);
    assert_eq!(euler, 4, "expected Euler characteristic 4 for two shells");
}

#[test]
fn thicken_cylinder_small_thickness() {
    // A very small thickness on a cylinder should still work.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, 1e-6, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn thicken_cylinder_negative_thickness() {
    // Negative thickness on a cylinder should still produce valid topology.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, -0.2, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(result.brep.faces.len(), 6);
}

#[test]
fn thicken_cylinder_large_thickness() {
    // Large thickness on a cylinder. Planar caps can offset arbitrarily,
    // and the cylindrical face offset should also succeed.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, 50.0, &ShellOptions::default(), tol());
    // Large offset on a cylinder may or may not succeed depending on
    // whether the inner cylinder radius goes negative. Accept either.
    match result {
        Ok(r) => {
            assert!(r.brep.faces.len() >= 4);
        }
        Err(e) => {
            let _ = format!("{e}");
        }
    }
}

#[test]
fn thicken_preserves_no_warnings_for_simple_cube() {
    // A simple cube thicken should produce no warnings.
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    assert!(
        result.warnings.is_empty(),
        "expected no warnings for simple cube thicken, got: {:?}",
        result.warnings
    );
}
