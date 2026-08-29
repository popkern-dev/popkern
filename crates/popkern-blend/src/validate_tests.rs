use super::*;
use crate::error::BlendError;
use crate::types::{ChamferKind, ChamferSpec, FilletRadius, FilletSpec, VertexMode};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_and_edge() -> (BRep, popkern_core::topology::EdgeHandle) {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    (brep, edge)
}

fn cube_edge_and_face() -> (
    BRep,
    popkern_core::topology::EdgeHandle,
    popkern_core::topology::FaceHandle,
) {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let face = brep.faces.keys().next().unwrap();
    (brep, edge, face)
}

// ---------- fillet spec validation ----------

#[test]
fn validate_fillet_specs_rejects_empty() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let err = validate_fillet_specs(&brep, &[]).unwrap_err();
    match err {
        BlendError::InvalidParameter {
            parameter, reason, ..
        } => {
            assert_eq!(parameter, "specs");
            assert!(reason.contains("at least one"));
        }
        other => panic!("expected InvalidParameter, got {other:?}"),
    }
}

#[test]
fn validate_fillet_specs_accepts_constant_positive_radius() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    assert!(validate_fillet_specs(&brep, &specs).is_ok());
}

#[test]
fn validate_fillet_specs_rejects_nan_radius() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(f64::NAN),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { value, reason, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert_eq!(value, "NaN");
    assert!(reason.contains("NaN"));
}

#[test]
fn validate_fillet_specs_rejects_negative_radius() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(-0.1),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(reason.contains("strictly positive"));
}

#[test]
fn validate_fillet_specs_rejects_zero_radius() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.0),
    }];
    assert!(matches!(
        validate_fillet_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_fillet_specs_rejects_infinite_radius() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(f64::INFINITY),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(reason.contains("finite"));
}

#[test]
fn validate_fillet_specs_rejects_variable_with_zero_start() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Variable {
            start: 0.0,
            end: 0.1,
        },
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { parameter, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(parameter.ends_with(".start"));
}

#[test]
fn validate_fillet_specs_rejects_variable_with_negative_end() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Variable {
            start: 0.1,
            end: -0.1,
        },
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { parameter, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(parameter.ends_with(".end"));
}

#[test]
fn validate_fillet_specs_rejects_bypoints_empty() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![]),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(reason.contains("at least one"));
}

#[test]
fn validate_fillet_specs_accepts_bypoints_single_point() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(0.5, 0.1)]),
    }];
    assert!(validate_fillet_specs(&brep, &specs).is_ok());
}

#[test]
fn validate_fillet_specs_rejects_bypoints_out_of_range_parameter() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(1.5, 0.1)]),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter {
        parameter, reason, ..
    } = &err
    else {
        panic!("expected InvalidParameter");
    };
    assert!(parameter.contains(".t"));
    assert!(reason.contains("[0, 1]"));
}

#[test]
fn validate_fillet_specs_rejects_bypoints_nan_parameter() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(f64::NAN, 0.1)]),
    }];
    assert!(matches!(
        validate_fillet_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_fillet_specs_rejects_bypoints_non_monotonic_parameter() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(0.5, 0.1), (0.3, 0.1)]),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(reason.contains("strictly increasing"));
}

#[test]
fn validate_fillet_specs_rejects_bypoints_duplicate_parameter() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(0.3, 0.1), (0.3, 0.1)]),
    }];
    assert!(matches!(
        validate_fillet_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_fillet_specs_rejects_bypoints_negative_radius() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(0.5, -0.1)]),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { parameter, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(parameter.contains(".r"));
}

#[test]
fn validate_fillet_specs_rejects_unknown_edge_handle() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let bogus = popkern_core::topology::EdgeHandle::default();
    let specs = vec![FilletSpec {
        edge: bogus,
        radius: FilletRadius::Constant(0.1),
    }];
    let err = validate_fillet_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { reason, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(reason.contains("does not exist"));
}

// ---------- chamfer spec validation ----------

#[test]
fn validate_chamfer_specs_rejects_empty() {
    let brep = BRep::unit_cube(tol()).unwrap();
    assert!(matches!(
        validate_chamfer_specs(&brep, &[]),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_chamfer_specs_accepts_equal_positive() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    assert!(validate_chamfer_specs(&brep, &specs).is_ok());
}

#[test]
fn validate_chamfer_specs_rejects_equal_nan() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(f64::NAN),
    }];
    assert!(matches!(
        validate_chamfer_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_chamfer_specs_rejects_two_distance_negative() {
    let (brep, edge, face) = cube_edge_and_face();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::TwoDistance {
            face_1: face,
            dist_1: 0.1,
            dist_2: -0.1,
        },
    }];
    assert!(matches!(
        validate_chamfer_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_chamfer_specs_rejects_two_distance_unknown_face() {
    let (brep, edge) = cube_and_edge();
    let bogus = popkern_core::topology::FaceHandle::default();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::TwoDistance {
            face_1: bogus,
            dist_1: 0.1,
            dist_2: 0.1,
        },
    }];
    let err = validate_chamfer_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter { parameter, .. } = &err else {
        panic!("expected InvalidParameter");
    };
    assert!(parameter.contains(".face_1"));
}

#[test]
fn validate_chamfer_specs_accepts_distance_angle_in_range() {
    let (brep, edge, face) = cube_edge_and_face();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::DistanceAngle {
            face,
            distance: 0.1,
            angle: std::f64::consts::FRAC_PI_4,
        },
    }];
    assert!(validate_chamfer_specs(&brep, &specs).is_ok());
}

#[test]
fn validate_chamfer_specs_rejects_distance_angle_zero_angle() {
    let (brep, edge, face) = cube_edge_and_face();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::DistanceAngle {
            face,
            distance: 0.1,
            angle: 0.0,
        },
    }];
    let err = validate_chamfer_specs(&brep, &specs).unwrap_err();
    let BlendError::InvalidParameter {
        parameter, reason, ..
    } = &err
    else {
        panic!("expected InvalidParameter");
    };
    assert!(parameter.contains(".angle"));
    assert!(reason.contains("(0, π)"));
}

#[test]
fn validate_chamfer_specs_rejects_distance_angle_pi_angle() {
    let (brep, edge, face) = cube_edge_and_face();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::DistanceAngle {
            face,
            distance: 0.1,
            angle: std::f64::consts::PI,
        },
    }];
    assert!(matches!(
        validate_chamfer_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_chamfer_specs_rejects_distance_angle_nan_angle() {
    let (brep, edge, face) = cube_edge_and_face();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::DistanceAngle {
            face,
            distance: 0.1,
            angle: f64::NAN,
        },
    }];
    assert!(matches!(
        validate_chamfer_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn validate_chamfer_specs_rejects_distance_angle_unknown_face() {
    let (brep, edge) = cube_and_edge();
    let bogus = popkern_core::topology::FaceHandle::default();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::DistanceAngle {
            face: bogus,
            distance: 0.1,
            angle: std::f64::consts::FRAC_PI_4,
        },
    }];
    assert!(matches!(
        validate_chamfer_specs(&brep, &specs),
        Err(BlendError::InvalidParameter { .. })
    ));
}

// Suppress "unused import" for VertexMode — kept because callers of
// these tests may extend the suite to test vertex-mode validation.
#[allow(dead_code)]
fn _vertex_mode_is_imported(_: VertexMode) {}
