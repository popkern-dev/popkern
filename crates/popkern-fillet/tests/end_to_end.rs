//! End-to-end integration tests for the fillet and chamfer public API.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ChamferKind, ChamferOptions, ChamferSpec, FilletOptions, FilletRadius, FilletSpec,
};

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
    let result = popkern_fillet::fillet_edge(&brep, edge, 0.1, default_tol());
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
    let result = popkern_fillet::fillet_edge(&brep, edge, 0.1, default_tol()).unwrap();

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
    let result = popkern_fillet::chamfer_edge(&brep, edge, 0.2, default_tol());
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
    let result = popkern_fillet::chamfer_edge(&brep, edge, 0.2, default_tol()).unwrap();

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
    let result = popkern_fillet::fillet_edges(&brep, &specs, &FilletOptions::default(), tol);
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
        Err(popkern_fillet::BlendError::PartialFailure { .. }) => {
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
    let result = popkern_fillet::chamfer_edges(&brep, &specs, &ChamferOptions::default(), tol);
    match &result {
        Ok(blend) => {
            assert!(
                blend.brep.faces.len() >= brep.faces.len(),
                "result should have >= original face count"
            );
            assert_euler_poincare(&blend.brep, "chamfer all cube edges");
            assert_structurally_valid(&blend.brep, "chamfer all cube edges");
        }
        Err(popkern_fillet::BlendError::PartialFailure { .. }) => {}
        Err(e) => panic!("unexpected error: {e:?}"),
    }
}

#[test]
fn fillet_empty_specs_returns_error() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let result = popkern_fillet::fillet_edges(&brep, &[], &FilletOptions::default(), default_tol());
    assert!(result.is_err(), "empty specs should return an error");
}

#[test]
fn chamfer_empty_specs_returns_error() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let result =
        popkern_fillet::chamfer_edges(&brep, &[], &ChamferOptions::default(), default_tol());
    assert!(result.is_err(), "empty specs should return an error");
}

#[test]
fn fillet_invalid_edge_handle() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let result = popkern_fillet::fillet_edge(&brep, bogus, 0.1, default_tol());
    assert!(result.is_err(), "invalid edge should return an error");
}

#[test]
fn chamfer_invalid_edge_handle() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let result = popkern_fillet::chamfer_edge(&brep, bogus, 0.1, default_tol());
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
    let result =
        popkern_fillet::fillet_edges(&brep, &[spec], &FilletOptions::default(), default_tol());
    // Variable radius may fail geometrically but should not panic.
    match &result {
        Ok(_) | Err(popkern_fillet::BlendError::PartialFailure { .. }) => {}
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
    let result =
        popkern_fillet::fillet_edges(&brep, &[spec], &FilletOptions::default(), default_tol());
    match &result {
        Ok(_) | Err(popkern_fillet::BlendError::PartialFailure { .. }) => {}
        Err(e) => panic!("unexpected error: {e:?}"),
    }
}

#[test]
fn blend_result_has_warnings_field() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_fillet::fillet_edge(&brep, edge, 0.1, default_tol());
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
fn assert_structurally_valid(brep: &BRep, label: &str) {
    let errors = brep.validate_all();
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
    assert!(
        structural.is_empty(),
        "{label}: structural errors: {structural:?}"
    );
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
        vertex_mode: popkern_fillet::VertexMode::Smooth,
    };

    let result = popkern_fillet::fillet_edges(&brep, &specs, &opts, tol);
    let blend = result.unwrap();

    // Must produce more faces than the original cube (6 + 3 fillets - merged = 8).
    assert!(
        blend.brep.faces.len() > brep.faces.len(),
        "should have more faces after 3-edge fillet"
    );
    assert_euler_poincare(&blend.brep, "3-edge fillet");
    assert_structurally_valid(&blend.brep, "3-edge fillet");

    // Smooth mode currently falls back to Mitred — verify a warning is emitted.
    let has_fallback = blend
        .edge_warnings
        .iter()
        .any(|(_, w)| matches!(w, popkern_fillet::BlendWarning::VertexModeFallback { .. }));
    assert!(
        has_fallback,
        "expected VertexModeFallback warning for Smooth mode"
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
        vertex_mode: popkern_fillet::VertexMode::Mitred,
    };

    let result = popkern_fillet::fillet_edges(&brep, &specs, &opts, tol);
    let blend = result.unwrap();

    assert_euler_poincare(&blend.brep, "3-edge mitred fillet");
    assert_structurally_valid(&blend.brep, "3-edge mitred fillet");

    // Mitred mode should produce no vertex-mode warnings.
    let has_fallback = blend
        .edge_warnings
        .iter()
        .any(|(_, w)| matches!(w, popkern_fillet::BlendWarning::VertexModeFallback { .. }));
    assert!(
        !has_fallback,
        "mitred mode should not produce VertexModeFallback warnings"
    );
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
        vertex_mode: popkern_fillet::VertexMode::Smooth,
    };

    let result = popkern_fillet::fillet_edges(&brep, &specs, &opts, tol);
    let blend = result.unwrap();

    assert_euler_poincare(&blend.brep, "2-edge fillet");
    assert_structurally_valid(&blend.brep, "2-edge fillet");
}
