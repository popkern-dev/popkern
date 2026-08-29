//! Tests for [`super::types`]. Extracted from `types.rs` to
//! keep it under the 500-line structural cap.

use super::types::*;

#[test]
fn fillet_options_default() {
    let opts = FilletOptions::default();
    assert!(opts.propagate_tangent);
    assert_eq!(opts.vertex_mode, VertexMode::Smooth);
}

#[test]
fn chamfer_options_default() {
    let opts = ChamferOptions::default();
    assert!(opts.propagate_tangent);
}

#[test]
fn fillet_radius_constant() {
    let r = FilletRadius::Constant(2.5);
    match r {
        FilletRadius::Constant(v) => assert!((v - 2.5).abs() < f64::EPSILON),
        _ => panic!("expected Constant"),
    }
}
