//! Analysis setup files: the TOML format a structural or thermal setup
//! is written in, parsed to the solver's typed setups and written back.
//!
//! The `popkern-analysis` CLI reads these files beside a STEP part;
//! popkern live sends the same text with its `simulate` request and
//! saves it with the project, so a setup made in the panel runs
//! unchanged from the command line (docs/architecture/analysis-integration.md).

#![warn(missing_docs)]
#![forbid(unsafe_code)]

pub mod setup_file;
pub mod thermal_file;

pub use setup_file::{SetupFileError, parse_setup, write_setup};
pub use thermal_file::{ThermalFileError, parse_thermal_setup, write_thermal_setup};

#[cfg(test)]
#[path = "setup_file_tests.rs"]
mod setup_file_tests;
