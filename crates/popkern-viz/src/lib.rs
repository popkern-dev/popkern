//! popkern-viz: interactive 3D debug viewer for BRep models.
//!
//! Provides both windowed (human) and headless (agent) modes for
//! inspecting BRep topology, geometry, and boolean operation results.

// Always available
pub mod backend;
pub mod camera;
pub mod colors;
pub mod commands;
pub mod error;
pub mod headless;
pub mod loader;
pub mod mesh_data;
pub mod picking;
pub mod scene;

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
