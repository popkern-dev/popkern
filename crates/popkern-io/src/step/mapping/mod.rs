//! Mapping layers that convert resolved STEP entities into popkern-core types.
//!
//! Each submodule handles a different domain:
//! - [`geometry`] — points, directions, placements, curves, surfaces
//! - [`topology`] — BRep construction from STEP topology entities
//! - [`product`] — product structure and assembly tree
//! - [`presentation`] — colors, layers, and visual metadata

pub mod geometry;
pub mod presentation;
pub mod product;
pub mod topology;
