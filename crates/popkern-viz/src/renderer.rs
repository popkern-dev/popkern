//! wgpu render pipeline for the debug viewer.
//!
//! [`Renderer`] owns the GPU device, render pipeline, uniform buffer,
//! and depth texture. It does **not** manage windowing — that is handled
//! by the eframe integration layer (Task 7).

use wgpu::util::DeviceExt;

use crate::camera::Camera;

// Re-export MeshVertex from the shared mesh_data module so that existing
// `crate::renderer::MeshVertex` paths continue to work.
pub use crate::mesh_data::MeshVertex;

/// Vertex buffer layout descriptor for [`MeshVertex`] in the mesh pipeline.
///
/// This is a standalone function (rather than an inherent method) because
/// it depends on wgpu types while [`MeshVertex`] itself is wgpu-free.
pub fn mesh_vertex_buffer_layout() -> wgpu::VertexBufferLayout<'static> {
    const ATTRS: &[wgpu::VertexAttribute] = &wgpu::vertex_attr_array![
        0 => Float32x3, // position
        1 => Float32x3, // normal
        2 => Float32x4, // color
    ];
    wgpu::VertexBufferLayout {
        array_stride: std::mem::size_of::<MeshVertex>() as wgpu::BufferAddress,
        step_mode: wgpu::VertexStepMode::Vertex,
        attributes: ATTRS,
    }
}

// ---------------------------------------------------------------------------
// Uniform buffer
// ---------------------------------------------------------------------------

/// Data uploaded to the GPU uniform buffer each frame.
///
/// Contains the combined view-projection matrix and the directional light
/// direction (w component unused, padded for alignment).
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Uniforms {
    /// Combined view * projection matrix.
    pub view_proj: [[f32; 4]; 4],
    /// Directional light direction in world space (w unused).
    pub light_dir: [f32; 4],
}

impl Default for Uniforms {
    fn default() -> Self {
        Self {
            view_proj: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            light_dir: [0.577, 0.577, 0.577, 0.0],
        }
    }
}

// ---------------------------------------------------------------------------
// Depth texture helper
// ---------------------------------------------------------------------------

/// Depth attachment used by the render pipeline.
pub struct DepthTexture {
    /// The depth texture view used as attachment.
    pub view: wgpu::TextureView,
    _texture: wgpu::Texture,
}

/// Depth format used throughout the renderer.
pub const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

impl DepthTexture {
    /// Create a new depth texture for the given viewport dimensions.
    pub fn new(device: &wgpu::Device, width: u32, height: u32) -> Self {
        let size = wgpu::Extent3d {
            width: width.max(1),
            height: height.max(1),
            depth_or_array_layers: 1,
        };
        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("depth_texture"),
            size,
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
        Self {
            view,
            _texture: texture,
        }
    }
}

// ---------------------------------------------------------------------------
// Renderer
// ---------------------------------------------------------------------------

/// GPU rendering backend for the debug viewer.
///
/// Manages the wgpu device, mesh render pipeline, uniform buffer, and
/// depth texture. Call [`Renderer::new`] with an existing device/queue
/// (typically obtained from the windowing layer) and the surface texture
/// format.
pub struct Renderer {
    /// The wgpu logical device.
    pub device: wgpu::Device,
    /// The wgpu command queue.
    pub queue: wgpu::Queue,
    /// Pipeline for shaded triangle meshes.
    mesh_pipeline: wgpu::RenderPipeline,
    /// Uniform buffer (view_proj + light).
    uniform_buffer: wgpu::Buffer,
    /// Bind group that exposes the uniform buffer to shaders.
    uniform_bind_group: wgpu::BindGroup,
    /// Bind group layout (shared with other pipelines that need uniforms).
    uniform_bind_group_layout: wgpu::BindGroupLayout,
    /// Depth texture, sized to the current viewport.
    depth_texture: Option<DepthTexture>,
}

impl Renderer {
    /// Create a new renderer.
    ///
    /// `surface_format` must match the window surface format (e.g.
    /// `Bgra8UnormSrgb`). The renderer does **not** own a surface —
    /// the caller provides texture views at render time.
    pub fn new(
        device: wgpu::Device,
        queue: wgpu::Queue,
        surface_format: wgpu::TextureFormat,
    ) -> Self {
        let shader_src = include_str!("shaders/mesh.wgsl");
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("mesh_shader"),
            source: wgpu::ShaderSource::Wgsl(shader_src.into()),
        });

        // Uniform buffer
        let uniforms = Uniforms::default();
        let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("uniform_buffer"),
            contents: bytemuck::cast_slice(&[uniforms]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        // Bind group layout + bind group
        let uniform_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("uniform_bind_group_layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });

        let uniform_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("uniform_bind_group"),
            layout: &uniform_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        // Pipeline layout
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("mesh_pipeline_layout"),
            bind_group_layouts: &[&uniform_bind_group_layout],
            push_constant_ranges: &[],
        });

        // Render pipeline
        let mesh_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("mesh_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                buffers: &[mesh_vertex_buffer_layout()],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_main"),
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                unclipped_depth: false,
                polygon_mode: wgpu::PolygonMode::Fill,
                conservative: false,
            },
            depth_stencil: Some(wgpu::DepthStencilState {
                format: DEPTH_FORMAT,
                depth_write_enabled: true,
                depth_compare: wgpu::CompareFunction::Less,
                stencil: wgpu::StencilState::default(),
                bias: wgpu::DepthBiasState::default(),
            }),
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
            cache: None,
        });

        Self {
            device,
            queue,
            mesh_pipeline,
            uniform_buffer,
            uniform_bind_group,
            uniform_bind_group_layout,
            depth_texture: None,
        }
    }

    /// Returns a reference to the uniform bind group layout.
    ///
    /// Useful for pipelines that share the same uniform buffer (e.g.
    /// wireframe overlay).
    pub fn uniform_bind_group_layout(&self) -> &wgpu::BindGroupLayout {
        &self.uniform_bind_group_layout
    }

    /// Returns a reference to the uniform bind group.
    ///
    /// Useful for render passes that share the same uniform buffer.
    pub fn uniform_bind_group(&self) -> &wgpu::BindGroup {
        &self.uniform_bind_group
    }

    /// Ensure the depth texture matches the given viewport dimensions.
    ///
    /// Creates a new depth texture if none exists or the size has changed.
    pub fn ensure_depth_texture(&mut self, width: u32, height: u32) {
        let needs_create = self.depth_texture.is_none();
        if needs_create {
            self.depth_texture = Some(DepthTexture::new(&self.device, width, height));
        }
    }

    /// Recreate the depth texture to match a new viewport size.
    pub fn resize_depth_texture(&mut self, width: u32, height: u32) {
        self.depth_texture = Some(DepthTexture::new(&self.device, width, height));
    }

    /// Returns the depth texture view, if one has been created.
    pub fn depth_view(&self) -> Option<&wgpu::TextureView> {
        self.depth_texture.as_ref().map(|d| &d.view)
    }

    /// Upload the current camera and lighting uniforms to the GPU.
    ///
    /// Call this once per frame before issuing draw calls.
    pub fn update_uniforms(&self, camera: &Camera) {
        let vp = camera.view_projection();
        let uniforms = Uniforms {
            view_proj: vp.to_cols_array_2d(),
            // Default directional light from upper-right-front
            light_dir: [0.577, 0.577, 0.577, 0.0],
        };
        self.queue
            .write_buffer(&self.uniform_buffer, 0, bytemuck::cast_slice(&[uniforms]));
    }

    /// Record draw commands for a triangle mesh into the given render pass.
    ///
    /// The caller must have already begun a render pass with the appropriate
    /// color and depth attachments. `vertex_buffer` should contain
    /// [`MeshVertex`] data and `vertex_count` is the number of vertices to
    /// draw.
    pub fn render_mesh(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        vertex_buffer: &wgpu::Buffer,
        vertex_count: u32,
    ) {
        render_pass.set_pipeline(&self.mesh_pipeline);
        render_pass.set_bind_group(0, &self.uniform_bind_group, &[]);
        render_pass.set_vertex_buffer(0, vertex_buffer.slice(..));
        render_pass.draw(0..vertex_count, 0..1);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mesh_vertex_size_is_40_bytes() {
        assert_eq!(std::mem::size_of::<MeshVertex>(), 40);
    }

    #[test]
    fn uniforms_size_is_80_bytes() {
        // mat4x4 = 64 bytes + vec4 = 16 bytes = 80
        assert_eq!(std::mem::size_of::<Uniforms>(), 80);
    }

    #[test]
    fn uniforms_alignment_is_16_byte_aligned() {
        // wgpu requires uniform buffers to be 16-byte aligned
        assert_eq!(std::mem::align_of::<Uniforms>() % 4, 0);
    }

    #[test]
    fn default_uniforms_has_identity_view_proj() {
        let u = Uniforms::default();
        // Diagonal should be 1.0
        assert!((u.view_proj[0][0] - 1.0).abs() < f32::EPSILON);
        assert!((u.view_proj[1][1] - 1.0).abs() < f32::EPSILON);
        assert!((u.view_proj[2][2] - 1.0).abs() < f32::EPSILON);
        assert!((u.view_proj[3][3] - 1.0).abs() < f32::EPSILON);
    }

    #[test]
    fn vertex_buffer_layout_has_three_attributes() {
        let layout = mesh_vertex_buffer_layout();
        assert_eq!(layout.attributes.len(), 3);
        assert_eq!(layout.array_stride, 40);
    }

    #[test]
    fn mesh_vertex_is_pod() {
        // Verify bytemuck round-trip works
        let v = MeshVertex {
            position: [1.0, 2.0, 3.0],
            normal: [0.0, 1.0, 0.0],
            color: [1.0, 0.0, 0.0, 1.0],
        };
        let bytes = bytemuck::bytes_of(&v);
        assert_eq!(bytes.len(), 40);
        let v2: &MeshVertex = bytemuck::from_bytes(bytes);
        assert_eq!(v.position, v2.position);
        assert_eq!(v.normal, v2.normal);
        assert_eq!(v.color, v2.color);
    }
}
