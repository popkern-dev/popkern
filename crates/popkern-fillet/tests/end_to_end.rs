//! End-to-end integration tests for the fillet and chamfer public API.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ChamferKind, ChamferOptions, ChamferSpec, FilletOptions, FilletRadius, FilletSpec,
};

fn default_tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn fillet_single_cube_edge() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
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
