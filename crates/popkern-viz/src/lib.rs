//! popkern-viz: interactive 3D debug viewer for BRep models.
//!
//! Provides both windowed (human) and headless (agent) modes for
//! inspecting BRep topology, geometry, and boolean operation results.

pub mod app;
pub mod camera;
pub mod colors;
pub mod commands;
pub mod error;
pub mod gpu_mesh;
pub mod headless;
pub mod loader;
pub mod mesh_data;
pub mod picking;
pub mod renderer;
pub mod scene;
pub mod ui;
pub mod wireframe;
