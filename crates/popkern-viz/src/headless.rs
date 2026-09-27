//! Headless (offscreen) rendering: render a BRep model to a PNG file
//! without opening a window.
//!
//! When the `gpu` feature is enabled, uses wgpu with no surface (headless
//! adapter), renders to an offscreen texture, reads back the pixels, and
//! encodes them as a PNG via the `image` crate.

use std::path::Path;

use crate::camera::Camera;
use crate::error::VizError;
use popkern_core::topology::BRep;

#[cfg(feature = "gpu")]
use crate::mesh_data::build_mesh_vertices;
#[cfg(feature = "gpu")]
use crate::renderer::Renderer;
#[cfg(feature = "gpu")]
use crate::wireframe;
#[cfg(feature = "gpu")]
use wgpu::util::DeviceExt;

/// The texture format used for headless rendering.
#[cfg(feature = "gpu")]
const HEADLESS_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;

/// Render a BRep model to a PNG file using the GPU without opening a window.
///
/// Creates a headless wgpu device, tessellates the model, renders it
/// with the given camera parameters, and writes the result as a PNG.
///
/// # Arguments
///
/// * `brep` - The BRep model to render.
/// * `camera` - Camera defining the viewpoint.
/// * `width` - Output image width in pixels.
/// * `height` - Output image height in pixels.
/// * `output_path` - Where to write the PNG file.
/// * `show_wireframe` - Whether to overlay wireframe edges.
#[cfg(feature = "gpu")]
pub fn render_to_png_gpu(
    brep: &BRep,
    camera: &Camera,
    width: u32,
    height: u32,
    output_path: &Path,
    show_wireframe: bool,
) -> Result<(), VizError> {
    let (device, queue) = create_headless_device()?;

    // Tessellate the model at display quality (angular facet bound).
    let options = popkern_mesh::TessellationOptions::display();
    let mesh =
        popkern_mesh::tessellate(brep, &options).map_err(|e| VizError::Render(e.to_string()))?;

    if mesh.triangles.is_empty() {
        return Err(VizError::Render("model has no triangles".into()));
    }

    // Build GPU vertex data with per-face colors
    let (vertices, indices) = build_mesh_vertices(brep, &mesh);

    // Create offscreen texture (render target)
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("headless_color"),
        size: wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: HEADLESS_FORMAT,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let color_view = texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Create renderer (owns pipeline, uniform buffer, etc.)
    let renderer = Renderer::new(device, queue, HEADLESS_FORMAT);

    // Upload camera uniforms
    let mut cam = camera.clone();
    cam.aspect = width as f32 / height as f32;
    renderer.update_uniforms(&cam);

    // Upload vertex + index buffers
    let vertex_buffer = renderer
        .device
        .create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("headless_vb"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });

    let _index_buffer = renderer
        .device
        .create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("headless_ib"),
            contents: bytemuck::cast_slice(&indices),
            usage: wgpu::BufferUsages::INDEX,
        });

    // Depth texture
    let depth_tex = crate::renderer::DepthTexture::new(&renderer.device, width, height);

    // Wireframe data (optional)
    let wireframe_data = if show_wireframe {
        wireframe::GpuWireframe::from_brep(&renderer.device, brep)
    } else {
        None
    };

    let line_pipeline = if wireframe_data.is_some() {
        Some(wireframe::create_line_pipeline(
            &renderer.device,
            HEADLESS_FORMAT,
            renderer.uniform_bind_group_layout(),
        ))
    } else {
        None
    };

    // Record and submit render commands
    let mut encoder = renderer
        .device
        .create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("headless_encoder"),
        });

    {
        let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("headless_pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: &color_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: 0.15,
                        g: 0.15,
                        b: 0.15,
                        a: 1.0,
                    }),
                    store: wgpu::StoreOp::Store,
                },
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: &depth_tex.view,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(1.0),
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            ..Default::default()
        });

        renderer.render_mesh(&mut pass, &vertex_buffer, vertices.len() as u32);

        if let (Some(wf), Some(lp)) = (&wireframe_data, &line_pipeline) {
            wireframe::render_wireframe(&mut pass, lp, renderer.uniform_bind_group(), wf);
        }
    }

    // Copy texture to a buffer for readback
    let bytes_per_pixel = 4u32;
    let unpadded_bytes_per_row = width * bytes_per_pixel;
    let align = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    let padded_bytes_per_row = unpadded_bytes_per_row.div_ceil(align) * align;

    let output_buffer = renderer.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("headless_readback"),
        size: (padded_bytes_per_row * height) as u64,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    encoder.copy_texture_to_buffer(
        wgpu::TexelCopyTextureInfo {
            texture: &texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        wgpu::TexelCopyBufferInfo {
            buffer: &output_buffer,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(padded_bytes_per_row),
                rows_per_image: Some(height),
            },
        },
        wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
    );

    renderer.queue.submit(std::iter::once(encoder.finish()));

    // Map buffer and read pixels
    let buffer_slice = output_buffer.slice(..);
    let (sender, receiver) = std::sync::mpsc::channel();
    buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
        let _ = sender.send(result);
    });
    renderer.device.poll(wgpu::Maintain::Wait);

    receiver
        .recv()
        .map_err(|e| VizError::Render(format!("buffer map recv failed: {e}")))?
        .map_err(|e| VizError::Render(format!("buffer map failed: {e}")))?;

    let data = buffer_slice.get_mapped_range();

    // Remove row padding and write PNG
    let mut pixels = Vec::with_capacity((width * height * bytes_per_pixel) as usize);
    for row in 0..height {
        let start = (row * padded_bytes_per_row) as usize;
        let end = start + (width * bytes_per_pixel) as usize;
        pixels.extend_from_slice(&data[start..end]);
    }
    drop(data);
    output_buffer.unmap();

    write_png(output_path, width, height, &pixels)?;

    Ok(())
}

/// Render a BRep to PNG using the CPU software rasterizer.
///
/// Works on any machine — no GPU required.
pub fn render_to_png_software(
    brep: &BRep,
    camera: &Camera,
    width: u32,
    height: u32,
    output_path: &Path,
    show_wireframe: bool,
) -> Result<(), VizError> {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::display())
        .map_err(|e| VizError::SoftwareRender(format!("tessellation failed: {e}")))?;
    let pixels = crate::cpu_render::render(brep, &mesh, camera, width, height, show_wireframe);
    write_png(output_path, width, height, &pixels)
}

/// Render a BRep to PNG, using GPU if available, software otherwise.
pub fn render_to_png(
    brep: &BRep,
    camera: &Camera,
    width: u32,
    height: u32,
    output_path: &Path,
    show_wireframe: bool,
) -> Result<(), VizError> {
    #[cfg(feature = "gpu")]
    {
        match render_to_png_gpu(brep, camera, width, height, output_path, show_wireframe) {
            Ok(()) => return Ok(()),
            Err(e) => {
                log::warn!("GPU render failed ({e}), falling back to software renderer");
            }
        }
    }
    render_to_png_software(brep, camera, width, height, output_path, show_wireframe)
}

/// Create a headless wgpu device and queue (no window surface).
#[cfg(feature = "gpu")]
fn create_headless_device() -> Result<(wgpu::Device, wgpu::Queue), VizError> {
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        force_fallback_adapter: false,
        compatible_surface: None,
    }))
    .ok_or_else(|| VizError::Render("no GPU adapter available".into()))?;

    let (device, queue) = pollster::block_on(adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("headless_device"),
            ..Default::default()
        },
        None,
    ))
    .map_err(|e| VizError::Render(format!("device creation failed: {e}")))?;

    Ok((device, queue))
}

/// Encode RGBA pixels (`width * height * 4`) as a PNG and write to
/// disk.
pub fn write_png(path: &Path, width: u32, height: u32, pixels: &[u8]) -> Result<(), VizError> {
    let img = image::RgbaImage::from_raw(width, height, pixels.to_vec())
        .ok_or_else(|| VizError::Render("pixel buffer size mismatch".into()))?;
    img.save(path)
        .map_err(|e| VizError::Render(format!("PNG write failed: {e}")))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_data::build_mesh_vertices;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn build_mesh_vertices_from_unit_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let options = popkern_mesh::TessellationOptions::default();
        let mesh = popkern_mesh::tessellate(&brep, &options).expect("tessellate succeeds");

        let (vertices, indices) = build_mesh_vertices(&brep, &mesh);
        assert!(!vertices.is_empty());
        assert!(!indices.is_empty());
        // Each triangle contributes 3 vertices
        assert_eq!(vertices.len(), mesh.triangles.len() * 3);
        assert_eq!(indices.len(), mesh.triangles.len() * 3);
    }

    #[test]
    fn build_mesh_vertices_empty_mesh() {
        let brep = BRep::new(Tolerance::default());
        let mesh = popkern_mesh::TriMesh::default();

        let (vertices, indices) = build_mesh_vertices(&brep, &mesh);
        assert!(vertices.is_empty());
        assert!(indices.is_empty());
    }

    #[test]
    fn write_png_creates_file() {
        let dir = std::env::temp_dir();
        let path = dir.join("popkern_test_headless.png");

        // 2x2 red image
        let pixels: Vec<u8> = [255, 0, 0, 255].repeat(4);
        write_png(&path, 2, 2, &pixels).expect("png write should succeed");

        assert!(path.exists());
        let meta = std::fs::metadata(&path).expect("metadata");
        assert!(meta.len() > 0);

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn write_png_bad_size_is_error() {
        let dir = std::env::temp_dir();
        let path = dir.join("popkern_test_bad_size.png");

        // Mismatched pixel count
        let pixels: Vec<u8> = vec![255, 0, 0, 255];
        let result = write_png(&path, 10, 10, &pixels);
        assert!(result.is_err());
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn headless_format_is_rgba8() {
        assert_eq!(HEADLESS_FORMAT, wgpu::TextureFormat::Rgba8UnormSrgb);
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn uniforms_size_check() {
        // Verify we can still construct uniforms for headless
        let u = crate::renderer::Uniforms::default();
        assert_eq!(std::mem::size_of_val(&u), 80);
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn depth_format_is_depth32float() {
        assert_eq!(
            crate::renderer::DEPTH_FORMAT,
            wgpu::TextureFormat::Depth32Float
        );
    }
}

#[cfg(test)]
#[path = "headless_tests.rs"]
mod headless_tests;
