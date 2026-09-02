//! PNG rendering of the active body through the CPU software renderer,
//! with the camera fitted to the body's bounds.

use std::path::Path;

use glam::Vec3;
use popkern_core::topology::BRep;
use popkern_viz::camera::Camera;

/// Render `brep` to `path`, isometric view fitted to its bounds.
pub fn render_png(
    brep: &BRep,
    path: &Path,
    width: u32,
    height: u32,
    edges: bool,
) -> Result<(), String> {
    let mut min = Vec3::splat(f32::INFINITY);
    let mut max = Vec3::splat(f32::NEG_INFINITY);
    for vertex in brep.vertices.values() {
        let p = Vec3::new(
            vertex.position.x as f32,
            vertex.position.y as f32,
            vertex.position.z as f32,
        );
        min = min.min(p);
        max = max.max(p);
    }
    let mut camera = Camera::default_iso();
    camera.up = Vec3::Z;
    camera.eye = Vec3::new(3.0, -3.0, 2.5);
    if min.x.is_finite() && max.x.is_finite() {
        camera.fit_to_bounds(min, max);
    }
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    popkern_viz::headless::render_to_png(brep, &camera, width, height, path, edges)
        .map_err(|e| e.to_string())
}
