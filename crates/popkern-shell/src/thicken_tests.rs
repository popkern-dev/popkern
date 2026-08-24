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
fn thicken_cube_large_thickness_declines() {
    // A 100.0 wall in a unit cube cannot exist: the inward half-offset
    // (50.0) inverts the body. The old pipeline "succeeded" by serving
    // an inside-out inner skin; the collapse gate declines with the
    // caller's thickness units (max total wall just under 1.0).
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 100.0, &ShellOptions::default(), tol());
    match result {
        Err(ShellError::OffsetTooLarge {
            requested,
            max_safe,
            ..
        }) => {
            assert!((requested - 100.0).abs() < 1e-9, "requested {requested}");
            assert!(
                max_safe > 0.9 && max_safe < 1.0,
                "unit cube max wall must sit just under 1.0, got {max_safe}"
            );
        }
        other => panic!("expected OffsetTooLarge decline, got {other:?}"),
    }
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
fn thicken_cylinder_declines_until_curved_recompute() {
    // Curved-body thicken is geometrically wrong today: edge recompute
    // keeps the un-offset clone curves for non-planar surface pairs
    // (the cap-rim circles stay at the original radius while both skins
    // moved), so the assembled body fails structural validation. The
    // serve gate turns silent wrong geometry into a diagnosed decline.
    // PK-324 un-narrows this — when curved recompute lands, these pins
    // flip to validate_all-clean serving tests with the old expected
    // counts (6 faces, 4 vertices, 6 edges, Euler 4).
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();

    // A wall thinner than the validator tolerance band serves honestly:
    // the stale rim curves (clone curves, PK-324 scope) deviate from
    // the offset skins by half the thickness, which at 1e-6 sits inside
    // tolerance — the geometry is RIGHT to within the contract.
    let tiny = thicken(&brep, 1e-6, &ShellOptions::default(), tol()).unwrap();
    assert_eq!(tiny.brep.faces.len(), 6);

    for thickness in [0.2, -0.2] {
        let result = thicken(&brep, thickness, &ShellOptions::default(), tol());
        match result {
            Err(ShellError::ValidationFailed { errors }) => {
                assert!(!errors.is_empty(), "thickness {thickness}");
            }
            Ok(_) => panic!(
                "cylinder thicken served at {thickness} — if curved recompute \
                 landed (PK-324), upgrade these pins to serving tests"
            ),
            Err(other) => panic!("thickness {thickness}: expected ValidationFailed, got {other:?}"),
        }
    }
}

#[test]
fn thicken_cylinder_past_radius_declines_offset_too_large() {
    // A 50.0 wall on a radius-1.0 cylinder: the inward half-offset (25.0)
    // exceeds the lateral face's curvature radius, so the offset surface
    // does not exist. Honest decline names the binding face; the old
    // dispatch silently NURBS-approximated inverted geometry.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, 50.0, &ShellOptions::default(), tol());
    assert!(
        matches!(result, Err(ShellError::OffsetTooLarge { .. })),
        "expected OffsetTooLarge decline, got {result:?}"
    );
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
