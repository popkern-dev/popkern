//! eframe application with wgpu 3D viewport and camera controls.
//!
//! [`VizApp`] implements [`eframe::App`], combining the wgpu renderer
//! with egui UI panels. Custom 3D rendering is injected into the egui
//! render pass via [`egui_wgpu::CallbackTrait`].

use eframe::egui_wgpu;
use eframe::egui_wgpu::CallbackResources;

use crate::camera::ViewPreset;
use crate::gpu_mesh::GpuMesh;
use crate::picking;
use crate::renderer::Renderer;
use crate::scene::{Selection, VizServer};
use crate::ui;
use crate::wireframe::{self, GpuWireframe};

// ---------------------------------------------------------------------------
// GPU resources stored in CallbackResources
// ---------------------------------------------------------------------------

/// GPU resources shared between the eframe app and the paint callback.
///
/// Stored in [`egui_wgpu::CallbackResources`] so that the
/// [`SceneRenderCallback`] can access them during the paint phase.
pub struct GpuResources {
    /// Our custom mesh renderer (owns pipeline, uniform buffer, etc.).
    pub renderer: Renderer,
    /// Line pipeline for wireframe overlay.
    pub line_pipeline: wgpu::RenderPipeline,
    /// Uploaded mesh data, if a model is loaded.
    pub gpu_mesh: Option<GpuMesh>,
    /// Uploaded wireframe data, if a model is loaded.
    pub wireframe: Option<GpuWireframe>,
}

// ---------------------------------------------------------------------------
// VizApp
// ---------------------------------------------------------------------------

/// Main application state for the popkern debug viewer.
///
/// Holds the scene (BRep model, camera, display options) and manages
/// GPU resource initialization via eframe's wgpu render state.
pub struct VizApp {
    /// Scene state (model, camera, selection, display options).
    pub scene: VizServer,
    /// Whether GPU resources have been initialized.
    gpu_initialized: bool,
    /// Stored render state from eframe for GPU access.
    render_state: Option<egui_wgpu::RenderState>,
}

impl VizApp {
    /// Create a new application with default scene state.
    pub fn new() -> Self {
        Self {
            scene: VizServer::new(),
            gpu_initialized: false,
            render_state: None,
        }
    }

    /// Initialize GPU resources using eframe's render state.
    ///
    /// Called once on the first frame when the wgpu device is available.
    fn init_gpu(&mut self, render_state: &egui_wgpu::RenderState) {
        let device = render_state.device.clone();
        let queue = render_state.queue.clone();
        let target_format = render_state.target_format;

        let renderer = Renderer::new(device.clone(), queue, target_format);

        let line_pipeline = wireframe::create_line_pipeline(
            &device,
            target_format,
            renderer.uniform_bind_group_layout(),
        );

        // Upload mesh/wireframe if a model is already loaded.
        let gpu_mesh = self
            .scene
            .brep
            .as_ref()
            .zip(self.scene.mesh.as_ref())
            .and_then(|(brep, mesh)| GpuMesh::from_brep_mesh(&device, brep, mesh));

        let gpu_wireframe = self
            .scene
            .brep
            .as_ref()
            .and_then(|brep| GpuWireframe::from_brep(&device, brep));

        let resources = GpuResources {
            renderer,
            line_pipeline,
            gpu_mesh,
            wireframe: gpu_wireframe,
        };

        // Store in eframe's callback resources type-map.
        render_state
            .renderer
            .write()
            .callback_resources
            .insert(resources);

        self.gpu_initialized = true;
        self.render_state = Some(render_state.clone());
    }

    /// Re-upload mesh and wireframe data after loading a new model.
    ///
    /// Called when the scene's BRep changes (e.g. via [`VizServer::load_brep`]).
    /// Will be wired up in Task 12 (model loading).
    #[allow(dead_code)]
    fn refresh_gpu_data(&self) {
        let Some(rs) = &self.render_state else {
            return;
        };
        let mut egui_renderer = rs.renderer.write();
        let Some(res) = egui_renderer.callback_resources.get_mut::<GpuResources>() else {
            return;
        };
        res.gpu_mesh = self
            .scene
            .brep
            .as_ref()
            .zip(self.scene.mesh.as_ref())
            .and_then(|(brep, mesh)| GpuMesh::from_brep_mesh(&rs.device, brep, mesh));
        res.wireframe = self
            .scene
            .brep
            .as_ref()
            .and_then(|brep| GpuWireframe::from_brep(&rs.device, brep));
    }

    /// Handle mouse and keyboard input for camera controls and picking.
    fn handle_input(&mut self, response: &egui::Response, ui: &egui::Ui, rect: egui::Rect) {
        // Orbit: left-drag
        if response.dragged_by(egui::PointerButton::Primary) && !ui.input(|i| i.modifiers.shift) {
            let delta = response.drag_delta();
            let sensitivity = 0.005;
            self.scene
                .camera
                .orbit(delta.x * sensitivity, delta.y * sensitivity);
        }

        // Pan: middle-drag or shift+left-drag
        let is_pan = response.dragged_by(egui::PointerButton::Middle)
            || (response.dragged_by(egui::PointerButton::Primary)
                && ui.input(|i| i.modifiers.shift));
        if is_pan {
            let delta = response.drag_delta();
            let sensitivity = 0.01;
            self.scene
                .camera
                .pan(-delta.x * sensitivity, delta.y * sensitivity);
        }

        // Zoom: scroll wheel
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll.abs() > 0.1 {
            self.scene.camera.zoom(scroll * 0.05);
        }

        // View presets: numpad / top row keys
        ui.input(|i| {
            if i.key_pressed(egui::Key::Num1) {
                self.scene.camera.set_view_preset(ViewPreset::Front);
            } else if i.key_pressed(egui::Key::Num2) {
                self.scene.camera.set_view_preset(ViewPreset::Back);
            } else if i.key_pressed(egui::Key::Num3) {
                self.scene.camera.set_view_preset(ViewPreset::Top);
            } else if i.key_pressed(egui::Key::Num4) {
                self.scene.camera.set_view_preset(ViewPreset::Bottom);
            } else if i.key_pressed(egui::Key::Num5) {
                self.scene.camera.set_view_preset(ViewPreset::Left);
            } else if i.key_pressed(egui::Key::Num6) {
                self.scene.camera.set_view_preset(ViewPreset::Right);
            } else if i.key_pressed(egui::Key::Num7) {
                self.scene.camera.set_view_preset(ViewPreset::Iso);
            }
        });

        // Face picking: click without drag to select a face.
        if response.clicked()
            && let Some(pos) = response.interact_pointer_pos()
        {
            self.handle_click(pos, rect);
        }
    }

    /// Perform face picking at the given screen position within the viewport rect.
    fn handle_click(&mut self, screen_pos: egui::Pos2, viewport_rect: egui::Rect) {
        let Some(mesh) = &self.scene.mesh else {
            self.scene.select(Selection::None);
            return;
        };

        // Convert screen position to normalized [0,1] coordinates within the viewport.
        let ndc_x = (screen_pos.x - viewport_rect.min.x) / viewport_rect.width();
        let ndc_y = (screen_pos.y - viewport_rect.min.y) / viewport_rect.height();

        match picking::pick_face(&self.scene.camera, (ndc_x, ndc_y), mesh) {
            Some(hit) => self.scene.select(Selection::Face(hit.face)),
            None => self.scene.select(Selection::None),
        }
    }
}

impl Default for VizApp {
    fn default() -> Self {
        Self::new()
    }
}

impl eframe::App for VizApp {
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        // Lazy-init GPU resources on first frame.
        if !self.gpu_initialized
            && let Some(rs) = frame.wgpu_render_state()
        {
            self.init_gpu(rs);
        }

        // 1. Top toolbar
        ui::toolbar(ctx, &mut self.scene);

        // 2. Bottom status bar
        ui::status_bar(ctx, &self.scene);

        // 3. Right inspector panel
        ui::side_panel(ctx, &self.scene);

        // 4. Central panel = 3D viewport
        egui::CentralPanel::default()
            .frame(egui::Frame::NONE.fill(egui::Color32::from_rgb(30, 30, 30)))
            .show(ctx, |ui| {
                let rect = ui.available_rect_before_wrap();
                let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());

                // Update camera aspect ratio to match viewport.
                if rect.width() > 0.0 && rect.height() > 0.0 {
                    self.scene.camera.aspect = rect.width() / rect.height();
                }

                self.handle_input(&response, ui, rect);

                // Update uniforms and issue paint callback.
                if let Some(rs) = &self.render_state {
                    {
                        let egui_renderer = rs.renderer.read();
                        if let Some(res) = egui_renderer.callback_resources.get::<GpuResources>() {
                            res.renderer.update_uniforms(&self.scene.camera);
                        }
                    }

                    let show_wireframe = self.scene.display.show_wireframe;
                    let cb = egui_wgpu::Callback::new_paint_callback(
                        rect,
                        SceneRenderCallback { show_wireframe },
                    );
                    ui.painter().add(cb);
                }

                // Request continuous repaint while interacting.
                if response.hovered() {
                    ctx.request_repaint();
                }
            });
    }
}

// ---------------------------------------------------------------------------
// Paint callback
// ---------------------------------------------------------------------------

/// Custom wgpu paint callback that renders the 3D scene.
///
/// Issued once per frame by [`VizApp::update`]. Draws the shaded mesh
/// and (optionally) the wireframe overlay into egui's render pass.
struct SceneRenderCallback {
    /// Whether to draw the wireframe overlay.
    show_wireframe: bool,
}

impl egui_wgpu::CallbackTrait for SceneRenderCallback {
    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        callback_resources: &CallbackResources,
    ) {
        let Some(res) = callback_resources.get::<GpuResources>() else {
            return;
        };

        // Draw shaded mesh.
        if let Some(gpu_mesh) = &res.gpu_mesh {
            res.renderer
                .render_mesh(render_pass, &gpu_mesh.vertex_buffer, gpu_mesh.num_indices);
        }

        // Draw wireframe overlay.
        if self.show_wireframe
            && let Some(wf) = &res.wireframe
        {
            wireframe::render_wireframe(
                render_pass,
                &res.line_pipeline,
                res.renderer.uniform_bind_group(),
                wf,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn viz_app_default_has_no_model() {
        let app = VizApp::new();
        assert!(app.scene.brep.is_none());
        assert!(!app.gpu_initialized);
    }

    #[test]
    fn viz_app_default_impl() {
        let app = VizApp::default();
        assert!(app.scene.brep.is_none());
        assert!(app.render_state.is_none());
    }

    #[test]
    fn viz_app_scene_is_accessible() {
        let mut app = VizApp::new();
        app.scene.display.show_wireframe = false;
        assert!(!app.scene.display.show_wireframe);
    }
}
