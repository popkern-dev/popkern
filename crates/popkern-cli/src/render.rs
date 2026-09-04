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
    let mut grow = |x: f64, y: f64, z: f64| {
        let p = Vec3::new(x as f32, y as f32, z as f32);
        min = min.min(p);
        max = max.max(p);
    };
    for vertex in brep.vertices.values() {
        grow(vertex.position.x, vertex.position.y, vertex.position.z);
    }
    // Vertices alone under-frame curved bodies — a helical tube has
    // only its seam endpoints as vertices, so the fitted camera
    // cropped the spring to a corner (PK-440). Boundary-curve samples
    // trace the real extent.
    // A helical wall's boundary winds many turns; sparse samples
    // under-frame it and the fitted camera ends up inside the model.
    for face in brep.faces.keys() {
        for p in brep.face_outer_loop_sample_positions(face, 96) {
            grow(p.x, p.y, p.z);
        }
    }
    let mut camera = Camera::default_iso();
    camera.up = Vec3::Z;
    camera.eye = Vec3::new(3.0, -3.0, 2.5);
    // The projection must match the output frame, or the picture lands
    // stretched and off-centre (default_iso carries a 16:9 aspect).
    camera.aspect = width as f32 / height.max(1) as f32;
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
