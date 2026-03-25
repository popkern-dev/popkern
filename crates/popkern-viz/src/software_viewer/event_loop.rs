//! winit event loop with camera input handling for the software viewer.
//!
//! Opens a window using winit and presents software-rendered frames via
//! softbuffer. Camera controls: left-drag = orbit, middle-drag or
//! shift+left-drag = pan, scroll = zoom, numpad = view presets.

use std::num::NonZeroU32;
use std::sync::Arc;

use winit::application::ApplicationHandler;
use winit::dpi::PhysicalSize;
use winit::event::{ElementState, MouseButton, MouseScrollDelta, WindowEvent};
use winit::event_loop::{ActiveEventLoop, EventLoop};
use winit::keyboard::{KeyCode, PhysicalKey};
use winit::window::{Window, WindowId};

use crate::camera::ViewPreset;
use crate::scene::VizServer;

use super::present::rgba_buf_to_u32_buf;

/// Application state for the software (non-GPU) interactive viewer.
struct SoftwareViewerApp {
    /// Scene state: model, mesh, camera, display options.
    scene: VizServer,
    /// The winit window (created on `resumed`).
    window: Option<Arc<Window>>,
    /// Softbuffer context (created on `resumed`).
    context: Option<softbuffer::Context<Arc<Window>>>,
    /// Softbuffer surface for pixel presentation.
    surface: Option<softbuffer::Surface<Arc<Window>, Arc<Window>>>,
    /// Whether the left mouse button is currently pressed.
    left_pressed: bool,
    /// Whether the middle mouse button is currently pressed.
    middle_pressed: bool,
    /// Whether a shift key is currently held.
    shift_held: bool,
    /// Last known cursor position for computing deltas.
    last_cursor: Option<(f64, f64)>,
    /// Whether the framebuffer needs to be re-rendered.
    needs_redraw: bool,
}

impl SoftwareViewerApp {
    /// Create a new viewer app wrapping the given scene.
    fn new(scene: VizServer) -> Self {
        Self {
            scene,
            window: None,
            context: None,
            surface: None,
            left_pressed: false,
            middle_pressed: false,
            shift_held: false,
            last_cursor: None,
            needs_redraw: true,
        }
    }

    /// Redraw the window contents using the CPU software renderer.
    fn redraw(&mut self) {
        let (Some(window), Some(surface)) = (self.window.as_ref(), self.surface.as_mut()) else {
            return;
        };

        let size = window.inner_size();
        let width = size.width;
        let height = size.height;

        if width == 0 || height == 0 {
            return;
        }

        // Render RGBA pixels via the CPU rasterizer (or a blank background).
        let rgba = match (self.scene.brep.as_ref(), self.scene.mesh.as_ref()) {
            (Some(brep), Some(mesh)) => crate::cpu_render::render(
                brep,
                mesh,
                &self.scene.camera,
                width,
                height,
                self.scene.display.show_wireframe,
            ),
            _ => {
                // No model loaded — fill with dark-gray background.
                [38u8, 38, 38, 255].repeat((width * height) as usize)
            }
        };

        // Convert RGBA to softbuffer's native u32 pixel format.
        let native = rgba_buf_to_u32_buf(&rgba);

        // Resize the surface buffer if needed, then copy pixels and present.
        if let (Some(w), Some(h)) = (NonZeroU32::new(width), NonZeroU32::new(height)) {
            if surface.resize(w, h).is_err() {
                return;
            }
            if let Ok(mut buffer) = surface.buffer_mut() {
                let len = buffer.len().min(native.len());
                buffer[..len].copy_from_slice(&native[..len]);
                let _ = buffer.present();
            }
        }

        self.needs_redraw = false;
    }

    /// Handle keyboard view-preset shortcuts (numpad keys).
    fn handle_key(&mut self, key_code: KeyCode) {
        match key_code {
            KeyCode::Numpad1 => self.scene.camera.set_view_preset(ViewPreset::Front),
            KeyCode::Numpad2 => self.scene.camera.set_view_preset(ViewPreset::Back),
            KeyCode::Numpad3 => self.scene.camera.set_view_preset(ViewPreset::Right),
            KeyCode::Numpad4 => self.scene.camera.set_view_preset(ViewPreset::Left),
            KeyCode::Numpad5 => self.scene.camera.set_view_preset(ViewPreset::Iso),
            KeyCode::Numpad7 => self.scene.camera.set_view_preset(ViewPreset::Top),
            KeyCode::Numpad8 => self.scene.camera.set_view_preset(ViewPreset::Bottom),
            _ => return,
        }
        self.needs_redraw = true;
    }
}

impl ApplicationHandler for SoftwareViewerApp {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        let attrs = Window::default_attributes()
            .with_title("popkern-viz (software)")
            .with_inner_size(PhysicalSize::new(1280u32, 720u32));

        let window = match event_loop.create_window(attrs) {
            Ok(w) => Arc::new(w),
            Err(e) => {
                log::error!("failed to create window: {e}");
                event_loop.exit();
                return;
            }
        };

        let context = match softbuffer::Context::new(window.clone()) {
            Ok(c) => c,
            Err(e) => {
                log::error!("failed to create softbuffer context: {e}");
                event_loop.exit();
                return;
            }
        };

        let surface = match softbuffer::Surface::new(&context, window.clone()) {
            Ok(s) => s,
            Err(e) => {
                log::error!("failed to create softbuffer surface: {e}");
                event_loop.exit();
                return;
            }
        };

        self.window = Some(window);
        self.context = Some(context);
        self.surface = Some(surface);
        self.needs_redraw = true;
    }

    fn window_event(
        &mut self,
        event_loop: &ActiveEventLoop,
        _window_id: WindowId,
        event: WindowEvent,
    ) {
        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }

            WindowEvent::Resized(size) => {
                if size.width > 0 && size.height > 0 {
                    self.scene.camera.aspect = size.width as f32 / size.height as f32;
                    self.needs_redraw = true;
                    if let Some(window) = self.window.as_ref() {
                        window.request_redraw();
                    }
                }
            }

            WindowEvent::MouseInput { state, button, .. } => {
                let pressed = state == ElementState::Pressed;
                match button {
                    MouseButton::Left => self.left_pressed = pressed,
                    MouseButton::Middle => self.middle_pressed = pressed,
                    _ => {}
                }
            }

            WindowEvent::CursorMoved { position, .. } => {
                let (cx, cy) = (position.x, position.y);
                if let Some((lx, ly)) = self.last_cursor {
                    let dx = (cx - lx) as f32;
                    let dy = (cy - ly) as f32;

                    let is_orbit = self.left_pressed && !self.shift_held;
                    let is_pan = self.middle_pressed || (self.left_pressed && self.shift_held);

                    if is_orbit {
                        self.scene.camera.orbit(dx * 0.01, dy * 0.01);
                        self.needs_redraw = true;
                    } else if is_pan {
                        self.scene.camera.pan(dx * 0.005, dy * 0.005);
                        self.needs_redraw = true;
                    }
                }
                self.last_cursor = Some((cx, cy));

                if self.needs_redraw
                    && let Some(window) = self.window.as_ref()
                {
                    window.request_redraw();
                }
            }

            WindowEvent::MouseWheel { delta, .. } => {
                let y = match delta {
                    MouseScrollDelta::LineDelta(_, y) => y,
                    MouseScrollDelta::PixelDelta(pos) => pos.y as f32 * 0.1,
                };
                if y.abs() > f32::EPSILON {
                    self.scene.camera.zoom(y * 0.5);
                    self.needs_redraw = true;
                    if let Some(window) = self.window.as_ref() {
                        window.request_redraw();
                    }
                }
            }

            WindowEvent::ModifiersChanged(modifiers) => {
                self.shift_held = modifiers.state().shift_key();
            }

            WindowEvent::KeyboardInput { event, .. } => {
                if event.state == ElementState::Pressed
                    && let PhysicalKey::Code(code) = event.physical_key
                {
                    self.handle_key(code);
                    if self.needs_redraw
                        && let Some(window) = self.window.as_ref()
                    {
                        window.request_redraw();
                    }
                }
            }

            WindowEvent::RedrawRequested => {
                if self.needs_redraw {
                    self.redraw();
                }
            }

            _ => {}
        }
    }
}

/// Run the software interactive viewer.
///
/// Opens a window and enters the event loop. Returns when the window is closed.
pub fn run(scene: VizServer) -> Result<(), crate::error::VizError> {
    let event_loop =
        EventLoop::new().map_err(|e| crate::error::VizError::Render(format!("event loop: {e}")))?;
    event_loop.set_control_flow(winit::event_loop::ControlFlow::Wait);
    let mut app = SoftwareViewerApp::new(scene);
    event_loop
        .run_app(&mut app)
        .map_err(|e| crate::error::VizError::Render(format!("event loop: {e}")))?;
    Ok(())
}
