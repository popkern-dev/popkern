//! The popkern live engine (design §6.5, phase 4a): the same
//! record/replay engine every other host runs, wrapped for a page.
//!
//! A [`Session`] runs a Rhai model script under host budgets and
//! answers with an [`Envelope`]: how the run ended (the one class
//! partition every host shares — [`popkern_rhai::Outcome`]), what the
//! script printed, the replay report, the declared parameters, and
//! one triangle buffer set per alive body. Parameter edits replay in
//! place; STEP and `.pkn` leave as bytes for a download; the emitted
//! Rhai is the "show as code" form.
//!
//! Only the JS boundary ([`bindings`]) is target-specific. Everything
//! else compiles and is tested natively, so the wasm module carries no
//! logic of its own — the node smoke test only has to prove that the
//! same code answers the same way through the boundary.

#![forbid(unsafe_code)]

/// The worker protocol's version (`web/live/PROTOCOL.md`). Bumped when
/// a verb's shape changes; the page checks it before speaking.
pub const PROTOCOL_VERSION: &str = "1";

pub mod bake;
pub mod envelope;
pub mod mesh_buffers;
pub mod protocol;
pub mod session;
pub mod session_assembly_timeline;
mod session_document;
pub mod session_project;
pub mod session_timeline;

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
pub mod bindings;

pub use envelope::{AssemblyView, BodyMesh, DefinitionView, Envelope, InstanceView, WallBudget};
pub use popkern_rhai::{HostOptions, Outcome};
pub use popkern_script::ContinuityCertification;
pub use session::{CertifyStep, RunOptions, Session};
