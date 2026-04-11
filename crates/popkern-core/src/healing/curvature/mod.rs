//! Curvature analysis for surfaces and curves.
//!
//! Provides differential-geometry tools built on the fundamental forms:
//! principal curvatures, Gaussian/mean curvature, minimum radius of curvature,
//! and curve curvature/torsion.

pub mod curve_curvature;
pub mod fundamental_forms;
pub mod surface_curvature;

pub use curve_curvature::{curve_curvature_at, curve_torsion_at};
pub use fundamental_forms::{FirstFundamentalForm, SecondFundamentalForm, ShapeOperator};
pub use surface_curvature::{
    SurfaceCurvature, gaussian_curvature, mean_curvature, min_radius_of_curvature,
    principal_curvatures, surface_curvature_at,
};
