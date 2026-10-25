//! Thermal and structural analysis solver for popkern parts.
//!
//! The solver's only geometric input is an `AnalysisModel` (closed,
//! face-tagged triangle meshes plus face metadata). It depends on
//! `popkern-math` only, so it stays wasm-clean and never sees popkern
//! topology types. See `docs/architecture/analysis-library-boundary.md`.

#![warn(missing_docs)]
#![forbid(unsafe_code)]

pub mod aggregate;
pub mod assess;
pub mod clip;
pub mod element;
pub mod error;
pub mod face_ref;
pub mod grid;
pub mod linalg;
pub mod material;
pub mod model;
pub mod moments;
pub mod options;
pub mod parallel;
pub mod recover;
pub mod report;
pub mod setup;
pub mod solve;
pub mod structural;
pub mod thermal;
pub mod thermal_setup;

#[cfg(test)]
#[path = "mesh_test_helpers.rs"]
mod mesh_test_helpers;
#[cfg(test)]
#[path = "model_tests.rs"]
mod model_tests;
#[cfg(test)]
#[path = "setup_tests.rs"]
mod setup_tests;

pub use error::{AnalysisError, Degeneracy};
pub use face_ref::{Expect, FaceRef, Selector, Tolerance};
pub use material::{Material, MaterialError, Property};
pub use model::{
    AnalysisModel, ConcaveEdge, FaceId, FaceInfo, FaceRole, ModelError, Provenance, RoleSource,
    SolidId, SolidMesh, Surface, Units,
};
pub use options::AnalysisOptions;
pub use report::SolveReport;
pub use setup::{
    Axis, FaceDef, LoadDirection, LoadSpec, SetupError, StructuralSetup, SupportKind, SupportSpec,
};
