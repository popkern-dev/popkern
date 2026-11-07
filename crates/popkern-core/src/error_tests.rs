use super::*;

#[test]
fn validation_error_display() {
    let err = ValidationError::LowValence {
        vertex: Default::default(),
        valence: 1,
    };
    let msg = err.to_string();
    assert!(msg.contains("valence"));
}

#[test]
fn core_error_from_validation() {
    let val_err = ValidationError::LowValence {
        vertex: Default::default(),
        valence: 1,
    };
    let core_err: CoreError = val_err.into();
    assert!(matches!(core_err, CoreError::Validation(_)));
}

#[test]
fn core_error_invalid_handle() {
    let err = CoreError::InvalidHandle("vertex 42 not found".into());
    assert!(err.to_string().contains("vertex 42"));
}

#[test]
fn core_error_precondition_failed() {
    let err = CoreError::PreconditionFailed("vertices must be in same loop".into());
    assert!(err.to_string().contains("same loop"));
}

#[test]
fn core_error_tolerance_violation() {
    let err = CoreError::ToleranceViolation("points too far apart".into());
    assert!(err.to_string().contains("too far apart"));
}

#[test]
fn core_error_invalid_nurbs() {
    let err = CoreError::InvalidNurbs("knot vector not non-decreasing".into());
    assert!(err.to_string().contains("knot vector"));
}

#[test]
fn validation_error_euler_poincare() {
    let err = ValidationError::EulerPoincare {
        solid: Default::default(),
        expected: 0,
        actual: 2,
    };
    assert!(err.to_string().contains("expected 0"));
    assert!(err.to_string().contains("got 2"));
}

#[test]
fn self_intersection_error_display() {
    let err = ValidationError::SelfIntersection {
        face_a: Default::default(),
        face_b: Default::default(),
    };
    assert!(err.to_string().contains("Self-intersection"));
}
