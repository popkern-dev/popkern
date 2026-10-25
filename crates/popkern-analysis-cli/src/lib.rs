//! Command-line runner: a STEP part and a setup file in, a report out.
//!
//! This crate touches files, STEP and output formats; the setup file
//! format itself is `popkern-analysis-setup`.

#![warn(missing_docs)]
#![forbid(unsafe_code)]

pub mod png;
pub mod run;
pub mod samples;
pub mod step_model;
pub mod vtu;
pub mod vtu_design;

// The setup file format lives in popkern-analysis-setup (shared with
// popkern live); these keep the CLI's module paths.
pub use popkern_analysis_setup::{setup_file, thermal_file};
