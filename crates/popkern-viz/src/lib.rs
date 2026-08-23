//! popkern-viz: interactive 3D debug viewer for BRep models.
//!
//! Provides both windowed (human) and headless (agent) modes for
//! inspecting BRep topology, geometry, and boolean operation results.

// The eframe/wgpu callback-resources type-map recurses deeply through
// auto-trait resolution for GpuResources; nightly rustc now lints the
// default depth (recursion_depth_exceeding_limit, future hard error).
#![recursion_limit = "256"]
// Always available
pub mod backend;
pub mod camera;
pub mod colors;
pub mod commands;
pub mod cpu_render;
pub mod error;
pub mod headless;
pub mod imprecision;
pub mod loader;
pub mod mesh_data;
pub mod picking;
pub mod scene;
pub mod software_viewer;

// GPU-only
#[cfg(feature = "gpu")]
pub mod app;
#[cfg(feature = "gpu")]
pub mod gpu_mesh;
#[cfg(feature = "gpu")]
pub mod renderer;
#[cfg(feature = "gpu")]
pub mod ui;
#[cfg(feature = "gpu")]
pub mod wireframe;
