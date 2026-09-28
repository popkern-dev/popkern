//! PNG rendering of the active body through the CPU software renderer,
//! with the camera fitted to the body's bounds.

use std::path::Path;

use glam::Vec3;
use popkern_core::topology::BRep;
use popkern_viz::camera::Camera;
use popkern_viz::zebra::ZebraNormals;
use serde::Serialize;

/// The isometric camera fitted to `brep`'s bounds for a `width` ×
/// `height` frame.
pub fn fitted_camera(brep: &BRep, width: u32, height: u32) -> Camera {
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
    camera
}

fn create_parent(path: &Path) -> Result<(), String> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Render `brep` to `path`, isometric view fitted to its bounds.
pub fn render_png(
    brep: &BRep,
    path: &Path,
    width: u32,
    height: u32,
    edges: bool,
) -> Result<(), String> {
    let camera = fitted_camera(brep, width, height);
    create_parent(path)?;
    popkern_viz::headless::render_to_png(brep, &camera, width, height, path, edges)
        .map_err(|e| e.to_string())
}

/// The mesh-vs-surface normal disagreement of a zebra render (PK-503),
/// in the units a reader compares against the display facet bound.
#[derive(Debug, Clone, Serialize)]
pub struct ZebraSummary {
    /// Which normals the written image reflects in.
    pub normals: &'static str,
    /// Stripe count.
    pub stripes: u32,
    /// Angular facet bound of the rendered mesh (degrees).
    pub facet_angle_deg: f64,
    /// Pixels covered by the part.
    pub covered_pixels: usize,
    /// Covered pixels with an evaluated surface normal on the face the
    /// mesh shows.
    pub evaluated_pixels: usize,
    /// Covered pixels whose surface view shows a neighbouring face
    /// (between the mesh seam and the true seam; not in the angles).
    pub reassigned_pixels: usize,
    /// Covered pixels whose surface point lies outside the visible
    /// facet's footprint (silhouette slivers; not in the angles).
    pub silhouette_pixels: usize,
    /// Tessellation accuracy: largest mesh-vs-surface normal angle at
    /// corresponding points (degrees).
    pub facet_max_angle_deg: f64,
    /// Its 99th percentile (degrees).
    pub facet_p99_angle_deg: f64,
    /// Evaluated pixels whose surface point is the view ray's hit.
    pub ray_hit_pixels: usize,
    /// Largest mesh-vs-surface normal angle (degrees).
    pub max_angle_deg: f64,
    /// Mean angle (degrees).
    pub mean_angle_deg: f64,
    /// 99th-percentile angle (degrees).
    pub p99_angle_deg: f64,
    /// Pixel `[x, y]` of the largest angle.
    pub worst_pixel: Option<[u32; 2]>,
    /// Surface kind of the face owning the worst pixel.
    pub worst_surface: Option<&'static str>,
    /// Pixels whose stripe band differs between the two views.
    pub band_mismatch_pixels: usize,
}

/// Render `brep` to `path` as a zebra view and summarize how far its
/// mesh and surface normals disagree.
pub fn render_zebra_png(
    brep: &BRep,
    path: &Path,
    width: u32,
    height: u32,
    stripes: u32,
    normals: ZebraNormals,
    facet_angle_deg: f64,
) -> Result<ZebraSummary, String> {
    let camera = fitted_camera(brep, width, height);
    create_parent(path)?;
    let tessellation = popkern_mesh::TessellationOptions {
        max_facet_angle: facet_angle_deg.to_radians(),
        ..popkern_mesh::TessellationOptions::display()
    };
    let report = popkern_viz::surface_analysis::render_zebra_to_png(
        brep,
        &camera,
        width,
        height,
        stripes,
        normals,
        &tessellation,
        path,
    )
    .map_err(|e| e.to_string())?;
    Ok(ZebraSummary {
        normals: match normals {
            ZebraNormals::Mesh => "mesh",
            ZebraNormals::Surface => "surface",
        },
        stripes,
        facet_angle_deg,
        covered_pixels: report.covered_pixels,
        evaluated_pixels: report.evaluated_pixels,
        reassigned_pixels: report.reassigned_pixels,
        silhouette_pixels: report.silhouette_pixels,
        facet_max_angle_deg: report.facet_max_angle_rad.to_degrees(),
        facet_p99_angle_deg: report.facet_p99_angle_rad.to_degrees(),
        ray_hit_pixels: report.ray_hit_pixels,
        max_angle_deg: report.max_angle_rad.to_degrees(),
        mean_angle_deg: report.mean_angle_rad.to_degrees(),
        p99_angle_deg: report.p99_angle_rad.to_degrees(),
        worst_pixel: report.worst.map(|w| [w.x, w.y]),
        worst_surface: report.worst.and_then(|w| {
            brep.faces
                .get(w.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.kind_name())
        }),
        band_mismatch_pixels: report.band_mismatch_pixels,
    })
}
