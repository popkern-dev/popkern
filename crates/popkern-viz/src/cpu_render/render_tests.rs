//! Additional CPU render edge-case tests.

use super::*;
use crate::camera::{Camera, ViewPreset};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate};

#[test]
fn render_empty_brep_with_wireframe_is_all_background() {
    let brep = BRep::new(Tolerance::default());
    let mesh = popkern_mesh::TriMesh::default();
    let camera = Camera::default_iso();
    let pixels = render(&brep, &mesh, &camera, 16, 16, true);
    for chunk in pixels.chunks(4) {
        assert_eq!(chunk, &[38, 38, 38, 255]);
    }
}

#[test]
fn render_cube_from_different_views() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");

    let presets = [
        ViewPreset::Front,
        ViewPreset::Back,
        ViewPreset::Top,
        ViewPreset::Left,
    ];
    for preset in &presets {
        let mut camera = Camera::default_iso();
        camera.set_view_preset(*preset);
        let pixels = render(&brep, &mesh, &camera, 32, 32, false);
        assert_eq!(pixels.len(), 32 * 32 * 4);
        // Should have some non-background pixels (model visible).
        let has_non_bg = pixels
            .chunks(4)
            .any(|px| px[0] != 38 || px[1] != 38 || px[2] != 38);
        assert!(has_non_bg, "cube should be visible from preset {preset:?}");
    }
}

#[test]
fn render_1x1_does_not_panic() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
    let camera = Camera::default_iso();
    let pixels = render(&brep, &mesh, &camera, 1, 1, false);
    assert_eq!(pixels.len(), 4);
}

#[test]
fn render_large_resolution_is_correct_size() {
    let brep = BRep::new(Tolerance::default());
    let mesh = popkern_mesh::TriMesh::default();
    let camera = Camera::default_iso();
    let pixels = render(&brep, &mesh, &camera, 200, 100, false);
    assert_eq!(pixels.len(), 200 * 100 * 4);
}

#[test]
fn render_sphere_produces_visible_pixels() {
    let brep = BRep::sphere(1.0, Tolerance::default()).expect("sphere builds");
    let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
    let mut camera = Camera::default_iso();
    camera.fit_to_bounds(
        glam::Vec3::new(-1.0, -1.0, -1.0),
        glam::Vec3::new(1.0, 1.0, 1.0),
    );
    let pixels = render(&brep, &mesh, &camera, 64, 64, false);
    let non_bg_count = pixels
        .chunks(4)
        .filter(|px| px[0] != 38 || px[1] != 38 || px[2] != 38)
        .count();
    assert!(
        non_bg_count > 0,
        "sphere render should produce visible pixels"
    );
}
