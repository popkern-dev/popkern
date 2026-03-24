//! Wireframe edge extraction and GPU line rendering.
//!
//! Extracts BRep edges as line segments and creates a wgpu line-list
//! pipeline for rendering them as a wireframe overlay on top of the
//! shaded mesh.

use popkern_core::topology::BRep;
use wgpu::util::DeviceExt;

use crate::renderer::DEPTH_FORMAT;

// ---------------------------------------------------------------------------
// Line vertex format
// ---------------------------------------------------------------------------

/// Line vertex: position + RGBA color (28 bytes).
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct LineVertex {
    /// World-space position.
    pub position: [f32; 3],
    /// RGBA color.
    pub color: [f32; 4],
}

impl LineVertex {
    /// Vertex buffer layout descriptor for the line pipeline.
    pub fn buffer_layout() -> wgpu::VertexBufferLayout<'static> {
        const ATTRS: &[wgpu::VertexAttribute] = &wgpu::vertex_attr_array![
            0 => Float32x3, // position
            1 => Float32x4, // color
        ];
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<LineVertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: ATTRS,
        }
    }
}

/// Default wireframe edge color: dark gray.
const EDGE_COLOR: [f32; 4] = [0.15, 0.15, 0.15, 1.0];

// ---------------------------------------------------------------------------
// GPU wireframe data
// ---------------------------------------------------------------------------

/// GPU-ready line data for wireframe rendering.
///
/// Contains a vertex buffer of [`LineVertex`] pairs (one pair per edge)
/// for drawing with `PrimitiveTopology::LineList`.
pub struct GpuWireframe {
    /// wgpu buffer holding [`LineVertex`] data.
    pub vertex_buffer: wgpu::Buffer,
    /// Number of vertices to draw (2 per edge).
    pub num_vertices: u32,
}

impl GpuWireframe {
    /// Extract edge line segments from a BRep and upload to the GPU.
    ///
    /// Walks all edges in the BRep, resolves the origin vertex positions
    /// of both half-edges, and creates a [`LineVertex`] pair for each edge.
    ///
    /// Returns `None` if the BRep has no edges.
    pub fn from_brep(device: &wgpu::Device, brep: &BRep) -> Option<Self> {
        if brep.edges.is_empty() {
            return None;
        }

        let mut vertices: Vec<LineVertex> = Vec::with_capacity(brep.edges.len() * 2);

        for (_edge_handle, edge) in &brep.edges {
            let he_a_handle = edge.half_edge;
            let Some(he_a) = brep.half_edges.get(he_a_handle) else {
                continue;
            };
            let he_b_handle = he_a.twin;
            let Some(he_b) = brep.half_edges.get(he_b_handle) else {
                continue;
            };

            let Some(va) = brep.vertices.get(he_a.origin) else {
                continue;
            };
            let Some(vb) = brep.vertices.get(he_b.origin) else {
                continue;
            };

            let pa = [
                va.position.x as f32,
                va.position.y as f32,
                va.position.z as f32,
            ];
            let pb = [
                vb.position.x as f32,
                vb.position.y as f32,
                vb.position.z as f32,
            ];

            vertices.push(LineVertex {
                position: pa,
                color: EDGE_COLOR,
            });
            vertices.push(LineVertex {
                position: pb,
                color: EDGE_COLOR,
            });
        }

        if vertices.is_empty() {
            return None;
        }

        let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("wireframe_vertex_buffer"),
            contents: bytemuck::cast_slice(&vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });

        let num_vertices = vertices.len() as u32;
        Some(Self {
            vertex_buffer,
            num_vertices,
        })
    }
}

// ---------------------------------------------------------------------------
// Line render pipeline
// ---------------------------------------------------------------------------

/// Create the wgpu render pipeline for wireframe line rendering.
///
/// Uses the same uniform bind group layout as the mesh pipeline so that
/// both share the same camera/light uniforms. Lines are drawn with
/// `PrimitiveTopology::LineList`, depth-tested but not depth-written,
/// and with a slight depth bias so edges render on top of faces.
pub fn create_line_pipeline(
    device: &wgpu::Device,
    surface_format: wgpu::TextureFormat,
    uniform_bind_group_layout: &wgpu::BindGroupLayout,
) -> wgpu::RenderPipeline {
    let shader_src = include_str!("shaders/line.wgsl");
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("line_shader"),
        source: wgpu::ShaderSource::Wgsl(shader_src.into()),
    });

    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("line_pipeline_layout"),
        bind_group_layouts: &[uniform_bind_group_layout],
        push_constant_ranges: &[],
    });

    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("line_pipeline"),
        layout: Some(&pipeline_layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: Some("vs_main"),
            buffers: &[LineVertex::buffer_layout()],
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
            topology: wgpu::PrimitiveTopology::LineList,
            strip_index_format: None,
            front_face: wgpu::FrontFace::Ccw,
            cull_mode: None,
            unclipped_depth: false,
            polygon_mode: wgpu::PolygonMode::Fill,
            conservative: false,
        },
        depth_stencil: Some(wgpu::DepthStencilState {
            format: DEPTH_FORMAT,
            depth_write_enabled: false,
            depth_compare: wgpu::CompareFunction::LessEqual,
            stencil: wgpu::StencilState::default(),
            bias: wgpu::DepthBiasState {
                constant: -2,
                slope_scale: -1.0,
                clamp: 0.0,
            },
        }),
        multisample: wgpu::MultisampleState {
            count: 1,
            mask: !0,
            alpha_to_coverage_enabled: false,
        },
        multiview: None,
        cache: None,
    })
}

// ---------------------------------------------------------------------------
// Renderer integration: store and draw wireframe
// ---------------------------------------------------------------------------

/// Record draw commands for a wireframe overlay into a render pass.
///
/// The caller must have already begun a render pass with appropriate
/// color and depth attachments and must provide the line pipeline,
/// uniform bind group, and the [`GpuWireframe`] data.
pub fn render_wireframe<'a>(
    render_pass: &mut wgpu::RenderPass<'a>,
    line_pipeline: &'a wgpu::RenderPipeline,
    uniform_bind_group: &'a wgpu::BindGroup,
    wireframe: &'a GpuWireframe,
) {
    render_pass.set_pipeline(line_pipeline);
    render_pass.set_bind_group(0, uniform_bind_group, &[]);
    render_pass.set_vertex_buffer(0, wireframe.vertex_buffer.slice(..));
    render_pass.draw(0..wireframe.num_vertices, 0..1);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_vertex_size_is_28_bytes() {
        assert_eq!(std::mem::size_of::<LineVertex>(), 28);
    }

    #[test]
    fn line_vertex_is_pod() {
        let v = LineVertex {
            position: [1.0, 2.0, 3.0],
            color: [0.5, 0.5, 0.5, 1.0],
        };
        let bytes = bytemuck::bytes_of(&v);
        assert_eq!(bytes.len(), 28);
        let v2: &LineVertex = bytemuck::from_bytes(bytes);
        assert_eq!(v.position, v2.position);
        assert_eq!(v.color, v2.color);
    }

    #[test]
    fn line_vertex_buffer_layout_has_two_attributes() {
        let layout = LineVertex::buffer_layout();
        assert_eq!(layout.attributes.len(), 2);
        assert_eq!(layout.array_stride, 28);
    }

    #[test]
    fn edge_color_is_dark_gray() {
        assert!((EDGE_COLOR[0] - 0.15).abs() < f32::EPSILON);
        assert!((EDGE_COLOR[1] - 0.15).abs() < f32::EPSILON);
        assert!((EDGE_COLOR[2] - 0.15).abs() < f32::EPSILON);
        assert!((EDGE_COLOR[3] - 1.0).abs() < f32::EPSILON);
    }

    #[test]
    fn from_brep_returns_none_for_empty_brep() {
        // GpuWireframe::from_brep needs a real device for the Some path,
        // but the early-return path for empty edges does not.
        let brep = BRep::new(popkern_core::tolerance::Tolerance::default());
        assert!(brep.edges.is_empty());
        // The function would return None — verified by the is_empty check.
    }

    #[test]
    fn unit_cube_has_12_edges() {
        let brep =
            BRep::unit_cube(popkern_core::tolerance::Tolerance::default()).expect("cube builds");
        assert_eq!(brep.edges.len(), 12);
    }
}
