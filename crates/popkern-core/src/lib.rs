#![warn(missing_docs)]

//! popkern-core: Core topology, geometry, and validation for the popkern BRep kernel.
//!
//! This crate provides the foundational data structures for boundary representation
//! (BRep) solid modeling: topological entities (Vertex, Edge, Face, Shell, Solid),
//! geometric entities (Point, Curve, Surface), and validation utilities.

pub mod constructors;
pub mod error;
pub mod euler;
pub mod geometry;
pub mod tolerance;
pub mod topology;
pub mod validation;
