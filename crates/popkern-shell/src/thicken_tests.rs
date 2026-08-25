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
    // Minkowski-faithful default (PK-325): rounded outer skin
    // (6 planes + 12 cylinder strips + 8 sphere caps) + sharp inner
    // skin (erosion keeps convex edges sharp) = 32 faces.
    assert_eq!(result.brep.faces.len(), 32);
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
    // A very small thickness (1e-6) should still produce a valid result.
    // The half-offset (5e-7) is below the arc-join resolution floor
    // (sub-tolerance rounding is indistinguishable from the sharp
    // join), so the skins stay sharp: outer (6) + inner (6) = 12 faces.
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
    // Topology should be the same as positive thickness.
    assert_eq!(result.brep.faces.len(), 32);
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
    // Thickened under the Minkowski-faithful default (PK-325):
    //   outer skin = rounded box: F 26, E 48, V 24
    //   inner skin = sharp erosion: F 6, E 12, V 8
    //   Euler per shell: 24-48+26 = 2 and 8-12+6 = 2 (sum 4).
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    let f = result.brep.faces.len();
    let e = result.brep.edges.len();
    let v = result.brep.vertices.len();
    assert_eq!(f, 32, "expected 32 faces (26 rounded outer + 6 inner)");
    assert_eq!(v, 32, "expected 32 vertices (24 outer + 8 inner)");
    assert_eq!(e, 60, "expected 60 edges (48 outer + 12 inner)");
    // Euler characteristic for two disjoint closed shells.
    let euler = (v as i64) - (e as i64) + (f as i64);
    assert_eq!(euler, 4, "expected Euler characteristic 4 for two shells");
}

#[test]
fn thicken_cylinder_serves_with_exact_curved_recompute() {
    // PK-322: the analytic recompute arms (coaxial circles by affine
    // map, straight-ruling seams, curved vertex projection) make the
    // cylinder thicken geometrically exact — the pin that documented
    // the stale-clone decline flips to a measured serving test.
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();

    // Rounded outer skin (PK-325 arc joins): 3 offset faces + 2 torus
    // half-windows per rim = 7 faces; sharp inner skin = 3.
    assert_eq!(result.brep.faces.len(), 10, "7 rounded outer + 3 inner");
    assert_eq!(result.brep.vertices.len(), 10, "8 outer + 2 inner");
    assert_eq!(result.brep.edges.len(), 16, "13 outer + 3 inner");
    assert!(
        result.brep.validate_all().is_empty(),
        "thickened cylinder must be validate_all-clean: {:?}",
        result.brep.validate_all()
    );

    // Measured skins: lateral radii 1.1 (outward half) and 0.9 (inward
    // half); every vertex sits exactly on its skin's lateral radius.
    let mut radii: Vec<f64> = result
        .brep
        .faces
        .values()
        .filter_map(|f| f.surface.as_ref())
        .filter_map(|s| match &s.geometry {
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    radii.sort_by(f64::total_cmp);
    assert_eq!(radii.len(), 2, "one lateral face per skin");
    assert!((radii[0] - 0.9).abs() < 1e-12, "inner radius {}", radii[0]);
    assert!((radii[1] - 1.1).abs() < 1e-12, "outer radius {}", radii[1]);
    for v in result.brep.vertices.values() {
        let r = (v.position.x * v.position.x + v.position.y * v.position.y).sqrt();
        // Skin laterals sit at 0.9/1.1; the rim-join corner vertices on
        // the cap side stay at the original radius 1.0.
        assert!(
            (r - 0.9).abs() < 1e-9 || (r - 1.0).abs() < 1e-9 || (r - 1.1).abs() < 1e-9,
            "vertex at radius {r}"
        );
    }

    // Negative thickness mirrors the walls; still serves clean.
    let neg = thicken(&brep, -0.2, &ShellOptions::default(), tol()).unwrap();
    assert!(neg.brep.validate_all().is_empty());
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
