//! Additional headless rendering edge-case tests.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn render_to_png_software_empty_brep_produces_background() {
    let brep = BRep::new(Tolerance::default());
    let camera = Camera::default_iso();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_empty_brep_sw.png");

    // Empty BRep has no faces -> tessellation produces empty mesh -> render
    // produces all-background image (no error because the software renderer
    // handles empty meshes gracefully).
    let result = render_to_png_software(&brep, &camera, 32, 32, &path, false);
    assert!(result.is_ok(), "empty brep should render: {result:?}");

    // Verify the image was written.
    assert!(path.exists(), "PNG file should be created");
    let meta = std::fs::metadata(&path).expect("metadata");
    assert!(meta.len() > 0, "PNG file should not be empty");

    std::fs::remove_file(&path).ok();
}

#[test]
fn render_to_png_software_empty_brep_with_wireframe() {
    let brep = BRep::new(Tolerance::default());
    let camera = Camera::default_iso();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_empty_brep_wf.png");

    let result = render_to_png_software(&brep, &camera, 16, 16, &path, true);
    assert!(
        result.is_ok(),
        "empty brep + wireframe should render: {result:?}"
    );

    std::fs::remove_file(&path).ok();
}

#[test]
fn render_to_png_software_unit_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let camera = Camera::default_iso();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_cube_sw.png");

    let result = render_to_png_software(&brep, &camera, 64, 64, &path, true);
    assert!(result.is_ok(), "cube should render: {result:?}");
    assert!(path.exists());

    std::fs::remove_file(&path).ok();
}

#[test]
fn render_to_png_fallback_produces_file() {
    // render_to_png tries GPU first, then falls back to software.
    // In CI (no GPU), this exercises the fallback path.
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let camera = Camera::default_iso();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_cube_fallback.png");

    let result = render_to_png(&brep, &camera, 32, 32, &path, false);
    assert!(result.is_ok(), "fallback render should succeed: {result:?}");
    assert!(path.exists());

    std::fs::remove_file(&path).ok();
}

#[test]
fn write_png_zero_dimensions_is_error() {
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_zero_dim.png");
    // Zero width/height with empty pixels.
    let result = write_png(&path, 0, 0, &[]);
    assert!(result.is_err(), "zero dimensions should error");
}

#[test]
fn write_png_1x1_pixel() {
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_1x1.png");
    let pixels = vec![128, 64, 32, 255];
    let result = write_png(&path, 1, 1, &pixels);
    assert!(result.is_ok(), "1x1 PNG should succeed");
    assert!(path.exists());
    std::fs::remove_file(&path).ok();
}

#[test]
fn render_software_tiny_resolution() {
    // Render at minimum resolution (1x1).
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let camera = Camera::default_iso();
    let dir = std::env::temp_dir();
    let path = dir.join("popkern_test_1x1_render.png");

    let result = render_to_png_software(&brep, &camera, 1, 1, &path, false);
    assert!(result.is_ok(), "1x1 render should succeed: {result:?}");

    std::fs::remove_file(&path).ok();
}
