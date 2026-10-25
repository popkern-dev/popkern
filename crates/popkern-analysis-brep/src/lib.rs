//! Adapter from a popkern BRep to the solver's `AnalysisModel`, and
//! resolution of geometric face selectors.
//!
//! This is the only crate that sees popkern topology or geometry types.

#![warn(missing_docs)]
#![forbid(unsafe_code)]

pub mod adapter;
pub mod features;
mod fold_check;
pub mod geometry;
pub mod orient;
pub mod resolve;
mod roles;
pub mod weld;

pub use adapter::{AdapterError, NamedModel, model_from_breps, model_from_breps_named};
pub use features::{FeatureTolerances, SolidFeatures, detect_features};
pub use resolve::{Resolution, resolve_faces};

#[cfg(test)]
#[path = "adapter_tests.rs"]
mod adapter_tests;
#[cfg(test)]
#[path = "orient_tests.rs"]
mod orient_tests;
