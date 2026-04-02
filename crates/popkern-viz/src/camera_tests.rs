//! Additional camera edge-case tests.

use super::*;

#[test]
fn view_matrix_is_finite_for_default_camera() {
    let cam = Camera::default_iso();
    let m = cam.view_matrix();
    for col in 0..4 {
        for row in 0..4 {
            assert!(
                m.col(col)[row].is_finite(),
                "view matrix element [{row},{col}] should be finite"
            );
        }
    }
}

#[test]
fn projection_matrix_is_finite_for_default_camera() {
    let cam = Camera::default_iso();
    let m = cam.projection_matrix();
    for col in 0..4 {
        for row in 0..4 {
            assert!(
                m.col(col)[row].is_finite(),
                "projection matrix element [{row},{col}] should be finite"
            );
        }
    }
}

#[test]
fn view_projection_combines_both_matrices() {
    let cam = Camera::default_iso();
    let vp = cam.view_projection();
    let expected = cam.projection_matrix() * cam.view_matrix();
    for col in 0..4 {
        for row in 0..4 {
            let diff = (vp.col(col)[row] - expected.col(col)[row]).abs();
            assert!(diff < 1e-5, "vp element [{row},{col}] mismatch: {diff}");
        }
    }
}

#[test]
fn camera_with_zero_length_direction_orbit_is_noop() {
    // Eye == target means zero-length view direction.
    let mut cam = Camera {
        eye: Vec3::ZERO,
        target: Vec3::ZERO,
        up: Vec3::Y,
        fov_y: std::f32::consts::FRAC_PI_4,
        aspect: 1.0,
        near: 0.1,
        far: 100.0,
    };
    // Orbit should be a no-op (early return due to zero distance).
    cam.orbit(0.1, 0.1);
    assert_eq!(cam.eye, Vec3::ZERO);
    assert_eq!(cam.target, Vec3::ZERO);
}

#[test]
fn camera_with_zero_length_direction_zoom_is_noop() {
    let mut cam = Camera {
        eye: Vec3::ZERO,
        target: Vec3::ZERO,
        up: Vec3::Y,
        fov_y: std::f32::consts::FRAC_PI_4,
        aspect: 1.0,
        near: 0.1,
        far: 100.0,
    };
    cam.zoom(1.0);
    assert_eq!(cam.eye, Vec3::ZERO);
}

#[test]
fn zoom_clamps_to_min_distance() {
    let mut cam = Camera::default_iso();
    // Zoom in by a huge amount — should clamp to MIN_DISTANCE.
    cam.zoom(10_000.0);
    let dist = (cam.eye - cam.target).length();
    assert!(
        dist >= MIN_DISTANCE - f32::EPSILON,
        "distance should be at least MIN_DISTANCE, got {dist}"
    );
}

#[test]
fn zoom_out_increases_distance() {
    let mut cam = Camera::default_iso();
    let initial_dist = (cam.eye - cam.target).length();
    cam.zoom(-5.0);
    let new_dist = (cam.eye - cam.target).length();
    assert!(
        new_dist > initial_dist,
        "zoom out should increase distance: {initial_dist} -> {new_dist}"
    );
}

#[test]
fn pan_moves_both_eye_and_target() {
    let mut cam = Camera::default_iso();
    let initial_eye = cam.eye;
    let initial_target = cam.target;
    cam.pan(1.0, 0.0);
    // Both should move by the same amount.
    let eye_delta = cam.eye - initial_eye;
    let target_delta = cam.target - initial_target;
    let diff = (eye_delta - target_delta).length();
    assert!(
        diff < 1e-5,
        "pan should move eye and target equally, diff={diff}"
    );
    // They should have actually moved.
    assert!(eye_delta.length() > 0.01, "pan should move the camera");
}

#[test]
fn orbit_preserves_distance() {
    let mut cam = Camera::default_iso();
    let initial_dist = (cam.eye - cam.target).length();
    cam.orbit(0.3, 0.2);
    let new_dist = (cam.eye - cam.target).length();
    assert!(
        (new_dist - initial_dist).abs() < 1e-4,
        "orbit should preserve distance: {initial_dist} vs {new_dist}"
    );
}

#[test]
fn set_view_preset_preserves_distance() {
    let mut cam = Camera::default_iso();
    let initial_dist = (cam.eye - cam.target).length();
    cam.set_view_preset(ViewPreset::Front);
    let new_dist = (cam.eye - cam.target).length();
    assert!(
        (new_dist - initial_dist).abs() < 1e-4,
        "preset should preserve distance: {initial_dist} vs {new_dist}"
    );
}

#[test]
fn set_view_preset_with_zero_distance_uses_default() {
    let mut cam = Camera {
        eye: Vec3::ZERO,
        target: Vec3::ZERO,
        up: Vec3::Y,
        fov_y: std::f32::consts::FRAC_PI_4,
        aspect: 1.0,
        near: 0.1,
        far: 100.0,
    };
    cam.set_view_preset(ViewPreset::Top);
    let dist = (cam.eye - cam.target).length();
    assert!(
        (dist - PRESET_DISTANCE).abs() < 1e-3,
        "zero-distance preset should use PRESET_DISTANCE, got {dist}"
    );
}

#[test]
fn all_view_presets_produce_finite_matrices() {
    let presets = [
        ViewPreset::Front,
        ViewPreset::Back,
        ViewPreset::Top,
        ViewPreset::Bottom,
        ViewPreset::Left,
        ViewPreset::Right,
        ViewPreset::Iso,
    ];
    for preset in &presets {
        let mut cam = Camera::default_iso();
        cam.set_view_preset(*preset);
        let vp = cam.view_projection();
        for col in 0..4 {
            for row in 0..4 {
                assert!(
                    vp.col(col)[row].is_finite(),
                    "VP element [{row},{col}] not finite for preset {preset:?}"
                );
            }
        }
    }
}

#[test]
fn fit_to_bounds_moves_camera_to_see_box() {
    let mut cam = Camera::default_iso();
    cam.fit_to_bounds(Vec3::new(-10.0, -10.0, -10.0), Vec3::new(10.0, 10.0, 10.0));
    // Target should be at the box center.
    assert!((cam.target - Vec3::ZERO).length() < 1e-3);
    // Eye should be far enough to see the box.
    let dist = (cam.eye - cam.target).length();
    assert!(
        dist > 10.0,
        "camera should be far enough to see 20-unit box"
    );
}

#[test]
fn fit_to_bounds_zero_fov_uses_fallback() {
    let mut cam = Camera {
        eye: Vec3::new(0.0, 0.0, 5.0),
        target: Vec3::ZERO,
        up: Vec3::Y,
        fov_y: 0.0, // degenerate zero FOV
        aspect: 1.0,
        near: 0.1,
        far: 100.0,
    };
    cam.fit_to_bounds(Vec3::new(-1.0, -1.0, -1.0), Vec3::new(1.0, 1.0, 1.0));
    // Should not panic; distance should be a reasonable fallback.
    let dist = (cam.eye - cam.target).length();
    assert!(dist.is_finite(), "distance should be finite after fit");
}

#[test]
fn camera_with_extreme_aspect_ratio_produces_finite_projection() {
    let cam = Camera {
        eye: Vec3::new(0.0, 0.0, 5.0),
        target: Vec3::ZERO,
        up: Vec3::Y,
        fov_y: std::f32::consts::FRAC_PI_4,
        aspect: 100.0, // very wide
        near: 0.1,
        far: 100.0,
    };
    let m = cam.projection_matrix();
    for col in 0..4 {
        for row in 0..4 {
            assert!(
                m.col(col)[row].is_finite(),
                "projection [{row},{col}] should be finite with extreme aspect"
            );
        }
    }
}
