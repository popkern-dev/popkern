//! CPU software rasterizer — produces RGBA pixel buffers without GPU.

pub mod rasterize;
pub mod shading;
pub mod wireframe;

use popkern_core::topology::BRep;
use popkern_mesh::TriMesh;

use crate::camera::Camera;
use crate::mesh_data;
use rasterize::{Framebuffer, ScreenVertex, rasterize_triangle};
use shading::shade_vertex;

/// Render a BRep to an RGBA pixel buffer using CPU software rasterization.
///
/// Returns a `Vec<u8>` of length `width * height * 4` in RGBA order.
pub fn render(
    brep: &BRep,
    mesh: &TriMesh,
    camera: &Camera,
    width: u32,
    height: u32,
    show_wireframe: bool,
) -> Vec<u8> {
    let (vertices, indices) = mesh_data::build_mesh_vertices(brep, mesh);
    let vp = camera.view_projection();

    let mut fb = Framebuffer::new(width, height);
    fb.clear([38, 38, 38, 255]);

    let w = width as f32;
    let h = height as f32;

    // Project and shade all vertices.
    let projected: Vec<Option<ScreenVertex>> = vertices
        .iter()
        .map(|v| project_and_shade(v, &vp, w, h))
        .collect();

    // Rasterize triangles.
    for tri in indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        let i0 = tri[0] as usize;
        let i1 = tri[1] as usize;
        let i2 = tri[2] as usize;

        if let (Some(sv0), Some(sv1), Some(sv2)) = (&projected[i0], &projected[i1], &projected[i2])
        {
            rasterize_triangle(&mut fb, sv0, sv1, sv2);
        }
    }

    // Wireframe overlay.
    if show_wireframe {
        let edges = mesh_data::extract_edge_pairs(brep);
        for (pa, pb) in &edges {
            if let (Some((x0, y0)), Some((x1, y1))) =
                (project_point(*pa, &vp, w, h), project_point(*pb, &vp, w, h))
            {
                wireframe::draw_line(&mut fb, x0, y0, x1, y1, wireframe::WIREFRAME_COLOR);
            }
        }
    }

    fb.to_rgba()
}

/// Project a world-space vertex through VP to screen space and shade it.
fn project_and_shade(
    v: &mesh_data::MeshVertex,
    vp: &glam::Mat4,
    width: f32,
    height: f32,
) -> Option<ScreenVertex> {
    let pos = glam::Vec4::new(v.position[0], v.position[1], v.position[2], 1.0);
    let clip = *vp * pos;

    if clip.w <= 0.0 {
        return None;
    }

    let ndc = clip / clip.w;
    let sx = (ndc.x + 1.0) * 0.5 * width;
    let sy = (1.0 - ndc.y) * 0.5 * height;
    let sz = (ndc.z + 1.0) * 0.5;

    let shaded = shade_vertex(v.normal, v.color);
    let color = [
        (shaded[0].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
        (shaded[1].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
        (shaded[2].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
        (shaded[3].clamp(0.0, 1.0) * 255.0 + 0.5) as u8,
    ];

    Some(ScreenVertex {
        x: sx,
        y: sy,
        z: sz,
        color,
    })
}

/// Project a world-space point to screen coordinates.
fn project_point(pos: [f32; 3], vp: &glam::Mat4, width: f32, height: f32) -> Option<(f32, f32)> {
    let p = glam::Vec4::new(pos[0], pos[1], pos[2], 1.0);
    let clip = *vp * p;

    if clip.w <= 0.0 {
        return None;
    }

    let ndc = clip / clip.w;
    let sx = (ndc.x + 1.0) * 0.5 * width;
    let sy = (1.0 - ndc.y) * 0.5 * height;

    Some((sx, sy))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::Camera;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    use popkern_mesh::{TessellationOptions, tessellate};

    #[test]
    fn render_unit_cube_produces_image() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
        let camera = Camera::default_iso();
        let pixels = render(&brep, &mesh, &camera, 64, 64, false);
        assert_eq!(pixels.len(), 64 * 64 * 4);
        // Should not be all-background
        let has_non_bg = pixels
            .chunks(4)
            .any(|px| px[0] != 38 || px[1] != 38 || px[2] != 38);
        assert!(has_non_bg, "image should contain rendered pixels");
    }

    #[test]
    fn render_with_wireframe() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
        let camera = Camera::default_iso();
        let pixels = render(&brep, &mesh, &camera, 64, 64, true);
        assert_eq!(pixels.len(), 64 * 64 * 4);
    }

    #[test]
    fn render_empty_brep_is_all_background() {
        let brep = BRep::new(Tolerance::default());
        let mesh = popkern_mesh::TriMesh::default();
        let camera = Camera::default_iso();
        let pixels = render(&brep, &mesh, &camera, 16, 16, false);
        for chunk in pixels.chunks(4) {
            assert_eq!(chunk, &[38, 38, 38, 255]);
        }
    }
}
