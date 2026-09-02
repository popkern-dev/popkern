//! The browser playground engine (design §6.5, phase 4a): the same
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

pub mod envelope;
pub mod mesh_buffers;
pub mod session;

#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
pub mod bindings;

pub use envelope::{BodyMesh, Envelope, WallBudget};
pub use popkern_rhai::{HostOptions, Outcome};
pub use session::{RunOptions, Session};
