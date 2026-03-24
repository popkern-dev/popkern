use glam::Vec3;
use popkern_viz::camera::{Camera, ViewPreset};

#[test]
fn default_iso_eye_is_not_at_origin() {
    let cam = Camera::default_iso();
    assert_ne!(cam.eye, Vec3::ZERO);
    assert_eq!(cam.eye, Vec3::new(3.0, 3.0, 3.0));
}

#[test]
fn view_matrix_is_invertible() {
    let cam = Camera::default_iso();
    let view = cam.view_matrix();
    let det = view.determinant();
    assert!(det.abs() > f32::EPSILON, "view matrix should be invertible");
    let inv = view.inverse();
    let product = view * inv;
    for i in 0..4 {
        for j in 0..4 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(
                (product.col(j)[i] - expected).abs() < 1e-5,
                "product[{i}][{j}] = {} (expected {expected})",
                product.col(j)[i]
            );
        }
    }
}

#[test]
fn orbit_changes_eye_position() {
    let mut cam = Camera::default_iso();
    let before = cam.eye;
    cam.orbit(0.3, 0.1);
    let after = cam.eye;
    assert_ne!(before, after, "orbit should change the eye position");

    // Distance to target should be preserved
    let dist_before = (before - cam.target).length();
    let dist_after = (after - cam.target).length();
    assert!(
        (dist_before - dist_after).abs() < 1e-4,
        "orbit should preserve distance"
    );
}

#[test]
fn pan_moves_both_eye_and_target() {
    let mut cam = Camera::default_iso();
    let eye_before = cam.eye;
    let target_before = cam.target;
    cam.pan(1.0, 0.5);
    assert_ne!(cam.eye, eye_before, "pan should move the eye");
    assert_ne!(cam.target, target_before, "pan should move the target");

    // The eye-target vector should be unchanged
    let offset_before = eye_before - target_before;
    let offset_after = cam.eye - cam.target;
    assert!(
        (offset_before - offset_after).length() < 1e-5,
        "pan should preserve eye-target offset"
    );
}

#[test]
fn zoom_changes_distance() {
    let mut cam = Camera::default_iso();
    let dist_before = (cam.eye - cam.target).length();
    cam.zoom(1.0);
    let dist_after = (cam.eye - cam.target).length();
    assert!(
        dist_after < dist_before,
        "positive zoom should decrease distance"
    );

    // Zoom out
    cam.zoom(-5.0);
    let dist_out = (cam.eye - cam.target).length();
    assert!(
        dist_out > dist_after,
        "negative zoom should increase distance"
    );
}

#[test]
fn zoom_clamps_minimum_distance() {
    let mut cam = Camera::default_iso();
    cam.zoom(10000.0);
    let dist = (cam.eye - cam.target).length();
    assert!(dist >= 0.01, "zoom should clamp to minimum distance");
}

#[test]
fn fit_to_bounds_centers_on_box() {
    let mut cam = Camera::default_iso();
    let min = Vec3::new(10.0, 10.0, 10.0);
    let max = Vec3::new(20.0, 20.0, 20.0);
    cam.fit_to_bounds(min, max);

    let center = Vec3::new(15.0, 15.0, 15.0);
    assert_eq!(cam.target, center, "target should be at box center");

    // Eye should be along the original look direction from center
    let dir = (cam.eye - cam.target).normalize();
    // Original direction was (3,3,3) - (0,0,0) normalized
    let expected_dir = Vec3::new(1.0, 1.0, 1.0).normalize();
    assert!(
        (dir - expected_dir).length() < 1e-4,
        "look direction should be preserved"
    );
}

#[test]
fn view_preset_front_looks_along_z() {
    let mut cam = Camera::default_iso();
    cam.set_view_preset(ViewPreset::Front);

    // For Front, eye should be on +Z side of target
    let dir = cam.eye - cam.target;
    assert!(
        dir.z > 0.0,
        "front preset: eye should be on +Z side of target"
    );
    assert!(dir.x.abs() < 1e-5, "front preset: no X offset");
    assert!(dir.y.abs() < 1e-5, "front preset: no Y offset");
}

#[test]
fn view_preset_preserves_distance() {
    let mut cam = Camera::default_iso();
    let dist_before = (cam.eye - cam.target).length();
    cam.set_view_preset(ViewPreset::Top);
    let dist_after = (cam.eye - cam.target).length();
    assert!(
        (dist_before - dist_after).abs() < 1e-4,
        "preset should preserve eye-target distance"
    );
}

#[test]
fn all_presets_produce_valid_view_matrix() {
    let presets = [
        ViewPreset::Front,
        ViewPreset::Back,
        ViewPreset::Top,
        ViewPreset::Bottom,
        ViewPreset::Left,
        ViewPreset::Right,
        ViewPreset::Iso,
    ];
    for preset in presets {
        let mut cam = Camera::default_iso();
        cam.set_view_preset(preset);
        let det = cam.view_matrix().determinant();
        assert!(
            det.abs() > f32::EPSILON,
            "preset {preset:?} should produce invertible view matrix"
        );
    }
}

#[test]
fn view_projection_combines_view_and_projection() {
    let cam = Camera::default_iso();
    let vp = cam.view_projection();
    let expected = cam.projection_matrix() * cam.view_matrix();
    for i in 0..4 {
        for j in 0..4 {
            assert!(
                (vp.col(j)[i] - expected.col(j)[i]).abs() < 1e-5,
                "view_projection should equal projection * view"
            );
        }
    }
}
