//! Cell operators on the uniform grid, built from fragment moments.

pub mod operators;
pub mod poly;
pub mod scalar;
pub mod thermal_strain;
pub mod vector;
pub mod vector_boundary;

pub use scalar::{Matrix8, ScalarTables, SurfaceMoments, surface_moments};
pub use vector::{ETA_EAS, Enhancement, Matrix24, VectorCell, VectorTables, lame};

#[cfg(test)]
#[path = "cut_eas_tests.rs"]
mod cut_eas_tests;
#[cfg(test)]
#[path = "scalar_tests.rs"]
mod scalar_tests;
#[cfg(test)]
#[path = "vector_tests.rs"]
mod vector_tests;
