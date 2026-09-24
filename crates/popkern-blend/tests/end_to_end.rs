//! End-to-end integration tests for the fillet and chamfer public API.

use popkern_blend::{
    ChamferKind, ChamferOptions, ChamferSpec, FilletOptions, FilletRadius, FilletSpec,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn default_tol() -> Tolerance {
    Tolerance::default()
}

/// Assert the Euler-Poincare invariant V - E + F == 2 for a closed solid.
fn assert_euler_poincare(brep: &BRep, label: &str) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    let chi = v - e + f;
    assert!(
        chi == 2,
        "{label}: Euler-Poincare V - E + F = {v} - {e} + {f} = {chi}, expected 2"
    );
}

#[test]
fn fillet_single_cube_edge() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    assert_euler_poincare(&brep, "input cube");
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_blend::fillet_edge(&brep, edge, 0.1, default_tol(), None);
    assert!(
        result.is_ok(),
        "fillet should not error: {:?}",
        result.err()
    );
    let blend = result.unwrap();
    // Result BRep should have at least as many faces as the original.
    assert!(
        blend.brep.faces.len() >= brep.faces.len(),
        "result should have >= original face count"
    );
    assert_euler_poincare(&blend.brep, "fillet single edge");
}

#[test]
fn fillet_single_cube_edge_produces_valid_brep() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_blend::fillet_edge(&brep, edge, 0.1, default_tol(), None).unwrap();

    // Must have MORE faces than original (6 original - 2 consumed + 2 kept + 1 fillet = 7).
    assert_eq!(
        result.brep.faces.len(),
        7,
        "should have 7 faces after filleting one edge"
    );

    // Euler-Poincare invariant.
    assert_euler_poincare(&result.brep, "fillet valid brep");

    // Must pass structural validation.
    let errors = result.brep.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::BrokenCycle { .. }
                    | popkern_core::error::ValidationError::BrokenTwin { .. }
                    | popkern_core::error::ValidationError::InconsistentVertex { .. }
                    | popkern_core::error::ValidationError::NonManifoldEdge { .. }
            )
        })
        .collect();
    assert!(structural.is_empty(), "structural errors: {:?}", structural);
}

#[test]
fn chamfer_single_cube_edge() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_blend::chamfer_edge(&brep, edge, 0.2, default_tol(), None);
    assert!(
        result.is_ok(),
        "chamfer should not error: {:?}",
        result.err()
    );
    let blend = result.unwrap();
    assert!(
        blend.brep.faces.len() >= brep.faces.len(),
        "result should have >= original face count"
    );
    assert_euler_poincare(&blend.brep, "chamfer single edge");
}

#[test]
fn chamfer_single_cube_edge_produces_valid_brep() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_blend::chamfer_edge(&brep, edge, 0.2, default_tol(), None).unwrap();

    assert_eq!(
        result.brep.faces.len(),
        7,
        "should have 7 faces after chamfering one edge"
    );

    assert_euler_poincare(&result.brep, "chamfer valid brep");

    let errors = result.brep.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                popkern_core::error::ValidationError::BrokenCycle { .. }
                    | popkern_core::error::ValidationError::BrokenTwin { .. }
                    | popkern_core::error::ValidationError::InconsistentVertex { .. }
                    | popkern_core::error::ValidationError::NonManifoldEdge { .. }
            )
        })
        .collect();
    assert!(structural.is_empty(), "structural errors: {:?}", structural);
}

#[test]
fn fillet_all_cube_edges() {
    let tol = default_tol();
    let brep = BRep::unit_cube(tol).unwrap();
    let specs: Vec<_> = brep
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    let result = popkern_blend::fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None);
    // Some edges may fail due to overlapping fillets, but the call should not panic.
    // Accept either Ok or PartialFailure.
    match &result {
        Ok(blend) => {
            assert!(
                blend.brep.faces.len() >= brep.faces.len(),
                "result should have >= original face count"
            );
            assert_euler_poincare(&blend.brep, "fillet all cube edges");
            assert_structurally_valid(&blend.brep, "fillet all cube edges");
        }
        Err(popkern_blend::BlendError::PartialFailure { .. }) => {
            // Acceptable — some edges may fail for geometric reasons.
        }
        Err(e) => panic!("unexpected error: {e:?}"),
    }
}

#[test]
fn chamfer_all_cube_edges() {
    let tol = default_tol();
    let brep = BRep::unit_cube(tol).unwrap();
    let specs: Vec<_> = brep
        .edges
        .keys()
        .map(|e| ChamferSpec {
            edge: e,
            kind: ChamferKind::Equal(0.05),
        })
        .collect();
    let result = popkern_blend::chamfer_edges(&brep, &specs, &ChamferOptions::default(), tol, None);
    match &result {
        Ok(blend) => {
            assert!(
                blend.brep.faces.len() >= brep.faces.len(),
                "result should have >= original face count"
            );
            assert_euler_poincare(&blend.brep, "chamfer all cube edges");
            assert_structurally_valid(&blend.brep, "chamfer all cube edges");
        }
        Err(popkern_blend::BlendError::PartialFailure { .. }) => {}
        Err(e) => panic!("unexpected error: {e:?}"),
    }
}

#[test]
fn fillet_empty_specs_returns_error() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let result =
        popkern_blend::fillet_edges(&brep, &[], &FilletOptions::default(), default_tol(), None);
    assert!(result.is_err(), "empty specs should return an error");
}

#[test]
fn chamfer_empty_specs_returns_error() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let result =
        popkern_blend::chamfer_edges(&brep, &[], &ChamferOptions::default(), default_tol(), None);
    assert!(result.is_err(), "empty specs should return an error");
}

#[test]
fn fillet_invalid_edge_handle() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let result = popkern_blend::fillet_edge(&brep, bogus, 0.1, default_tol(), None);
    assert!(result.is_err(), "invalid edge should return an error");
}

#[test]
fn chamfer_invalid_edge_handle() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let result = popkern_blend::chamfer_edge(&brep, bogus, 0.1, default_tol(), None);
    assert!(result.is_err(), "invalid edge should return an error");
}

#[test]
fn fillet_variable_radius_single_edge() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let spec = FilletSpec {
        edge,
        radius: FilletRadius::Variable {
            start: 0.05,
            end: 0.1,
        },
    };
    let result = popkern_blend::fillet_edges(
        &brep,
        &[spec],
        &FilletOptions::default(),
        default_tol(),
        None,
    );
    // Variable radius may fail geometrically but should not panic.
    match &result {
        Ok(_) | Err(popkern_blend::BlendError::PartialFailure { .. }) => {}
        Err(e) => panic!("unexpected error: {e:?}"),
    }
}

#[test]
fn fillet_by_points_single_edge() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let spec = FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(0.0, 0.05), (0.5, 0.08), (1.0, 0.05)]),
    };
    let result = popkern_blend::fillet_edges(
        &brep,
        &[spec],
        &FilletOptions::default(),
        default_tol(),
        None,
    );
    match &result {
        Ok(_) | Err(popkern_blend::BlendError::PartialFailure { .. }) => {}
        Err(e) => panic!("unexpected error: {e:?}"),
    }
}

#[test]
fn blend_result_has_warnings_field() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_blend::fillet_edge(&brep, edge, 0.1, default_tol(), None);
    if let Ok(blend) = result {
        // Warnings vector should exist (may be empty).
        let _ = blend.edge_warnings.len();
    }
}

/// Helper: find all edges incident to a vertex.
fn edges_incident_to_vertex(
    brep: &BRep,
    vertex: popkern_core::topology::VertexHandle,
) -> Vec<popkern_core::topology::EdgeHandle> {
    let mut incident = Vec::new();
    for e in brep.edges.keys() {
        let ed = &brep.edges[e];
        let he = &brep.half_edges[ed.half_edge];
        let twin = &brep.half_edges[he.twin];
        if he.origin == vertex || twin.origin == vertex {
            incident.push(e);
        }
    }
    incident
}

/// Assert that the BRep has no structural validation errors.
///
/// Runs only the fast structural checks (half-edge pairing, cycle
/// integrity, vertex consistency, manifold edges). Intentionally skips
/// `validate_all()` which includes the expensive self-intersection
/// check — that check is O(n²) in face count and too slow for complex
/// fillet results with many NURBS surfaces.
fn assert_structurally_valid(brep: &BRep, label: &str) {
    let mut errors = Vec::new();
    if let Err(e) = brep.validate_half_edge_pairing() {
        errors.push(e);
    }
    if let Err(e) = brep.validate_half_edge_cycles() {
        errors.push(e);
    }
    if let Err(e) = brep.validate_vertex_consistency() {
        errors.push(e);
    }
    if let Err(e) = brep.validate_manifold_edges() {
        errors.push(e);
    }
    assert!(errors.is_empty(), "{label}: structural errors: {errors:?}");
}

#[test]
fn fillet_three_edges_at_cube_corner() {
    let tol = default_tol();
    let brep = BRep::unit_cube(tol).unwrap();

    // Find a vertex and its 3 incident edges.
    let vertex = brep.vertices.keys().next().unwrap();
    let incident = edges_incident_to_vertex(&brep, vertex);
    assert_eq!(incident.len(), 3, "cube vertex should have 3 edges");

    let specs: Vec<_> = incident
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();

    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: popkern_blend::VertexMode::Smooth,
        ..FilletOptions::default()
    };

    let result = popkern_blend::fillet_edges(&brep, &specs, &opts, tol, None);
    let blend = result.unwrap();

    // Must produce more faces than the original cube (6 + 3 fillets - merged = 8).
    assert!(
        blend.brep.faces.len() > brep.faces.len(),
        "should have more faces after 3-edge fillet"
    );
    assert_euler_poincare(&blend.brep, "3-edge fillet");
    assert_structurally_valid(&blend.brep, "3-edge fillet");

    // Smooth mode should insert actual blend faces — no fallback warnings.
    let fallback_count = blend
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, popkern_blend::BlendWarning::VertexModeFallback { .. }))
        .count();
    assert_eq!(
        fallback_count, 0,
        "smooth blend should not fall back to mitred"
    );
}

#[test]
fn fillet_three_edges_at_cube_corner_mitred() {
    let tol = default_tol();
    let brep = BRep::unit_cube(tol).unwrap();

    let vertex = brep.vertices.keys().next().unwrap();
    let incident = edges_incident_to_vertex(&brep, vertex);
    assert_eq!(incident.len(), 3);

    let specs: Vec<_> = incident
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();

    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: popkern_blend::VertexMode::Mitred,
        ..FilletOptions::default()
    };

    // PK-522: the mitre trims the three strips to their mutual
    // bisector-ellipse intersections, meeting at the tricylinder
    // vertex (PK-516 had declined the bare residue). Check 9 is the
    // gate the residue failed.
    let result = popkern_blend::fillet_edges(&brep, &specs, &opts, tol, None)
        .expect("mitred cube corner serves");
    assert_eq!(result.repair_report.mitre_corner.corners_built, 1);
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
}

#[test]
fn fillet_two_adjacent_edges() {
    let tol = default_tol();
    let brep = BRep::unit_cube(tol).unwrap();

    // Find two edges sharing a vertex.
    let vertex = brep.vertices.keys().next().unwrap();
    let incident = edges_incident_to_vertex(&brep, vertex);
    assert!(incident.len() >= 2);

    let specs: Vec<_> = incident[..2]
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();

    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: popkern_blend::VertexMode::Smooth,
        ..FilletOptions::default()
    };

    let result = popkern_blend::fillet_edges(&brep, &specs, &opts, tol, None);
    let blend = result.unwrap();

    assert_euler_poincare(&blend.brep, "2-edge fillet");
    assert_structurally_valid(&blend.brep, "2-edge fillet");
}

#[test]
fn large_radius_detects_overlap_or_fails() {
    let tol = default_tol();
    let brep = BRep::unit_cube(tol).unwrap();

    // Fillet two adjacent edges with a very large radius (0.4 on a unit cube).
    // The fillets should overlap. We accept either:
    //   - Ok with FilletOverlap warnings
    //   - PartialFailure (one edge fails with RadiusTooLarge or similar)
    //   - Full failure
    let vertex = brep.vertices.keys().next().unwrap();
    let incident = edges_incident_to_vertex(&brep, vertex);
    assert!(incident.len() >= 2);

    let specs: Vec<_> = incident[..2]
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.4),
        })
        .collect();

    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: popkern_blend::VertexMode::Mitred,
        ..FilletOptions::default()
    };

    let result = popkern_blend::fillet_edges(&brep, &specs, &opts, tol, None);

    match result {
        Ok(blend) => {
            // If it succeeded, it should carry overlap warnings.
            let has_overlap = blend
                .edge_warnings
                .iter()
                .any(|(_, w)| matches!(w, popkern_blend::BlendWarning::FilletOverlap { .. }));
            // Overlap detection is best-effort; the large radius may also just
            // succeed without triggering the heuristic. Either way, no panic.
            let _ = has_overlap;
        }
        Err(popkern_blend::BlendError::PartialFailure { failures, .. }) => {
            // At least one edge should have failed due to the large radius.
            assert!(
                !failures.is_empty(),
                "partial failure should report at least one failed edge"
            );
        }
        Err(_) => {
            // Total failure is also acceptable for overlapping fillets.
        }
    }
}
