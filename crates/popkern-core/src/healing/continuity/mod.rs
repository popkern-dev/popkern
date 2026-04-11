//! Continuity analysis for surface boundaries.
//!
//! Measures cross-boundary continuity (G0/G1/G2) between adjacent faces
//! and classifies edges by their achieved continuity level. Smooth chains
//! of edges are detected for downstream fillet and blending operations.

pub mod classification;
pub mod measurement;

pub use classification::{
    EdgeClassification, SmoothChain, classify_all_edges, classify_edge_smoothness,
    detect_smooth_chains,
};
pub use measurement::{BoundarySample, ContinuityMeasurement, cross_boundary_continuity};
