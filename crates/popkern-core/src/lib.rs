#![warn(missing_docs)]

//! popkern-core: Core topology, geometry, and validation for the popkern BRep kernel.
//!
//! This crate provides the foundational data structures for boundary representation
//! (BRep) solid modeling: topological entities (Vertex, Edge, Face, Shell, Solid),
//! geometric entities (Point, Curve, Surface), and validation utilities.

pub mod adaptive_precision;
pub mod cancellation;
pub mod check_finite;
pub mod constructors;
pub mod diagnostics;
pub mod error;
pub mod euler;
pub mod geometry;
pub mod healing;
pub mod history;
pub mod naming;
pub mod naming_criteria;
pub mod naming_resolve;
pub mod precision;
pub mod progress;
pub mod tolerance;
pub mod topology;
pub mod trace;
pub mod validation;
