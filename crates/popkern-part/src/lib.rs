//! Code-first part modeling over the popkern BRep kernel.
//!
//! A part is a short, readable program: sketch profiles, extrude
//! them, combine solids, dress edges up, gate the result, export it.
//! The API is a thin facade — every operation delegates to the
//! kernel crates (`popkern-sweep`, `popkern-bool`, `popkern-fillet`,
//! `popkern-mesh`, `popkern-io`) and preserves their diagnostics.
//!
//! ```no_run
//! use popkern_part::prelude::*;
//!
//! fn boss_plate() -> Result<Part, PartError> {
//!     let plate = Profile::on(Plane::XY).rect(2.4, 1.6).extrude(0.4)?;
//!     let boss = Profile::on(Plane::xy_at(0.2)).circle(0.35).extrude(0.8)?;
//!     plate
//!         .union(&boss)?
//!         .fillet(0.08, |e| {
//!             e.circular()
//!                 .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
//!         })
//! }
//! ```
//!
//! Design: `docs/architecture/part-modeling-api-design.md`. The API
//! follows the repo's UI-forward rules — `Result`-first (no deferred
//! errors), selector queries that treat EMPTY selections as errors,
//! feature logging for replay, and error payloads that carry the
//! values and limits a GUI needs to explain itself.

pub mod dressup;
mod engrave;
mod error;
mod feature;
mod part;
mod part_polyhedron;
pub mod poly;
mod profile;
mod select;
mod shaping;
#[cfg(feature = "text")]
pub mod text;

#[cfg(test)]
#[path = "part_tests.rs"]
mod part_tests;

pub use engrave::{
    EngravePlacement, FaceFrame, HoledPolygon, engrave_polygons, face_frame, polygons_bbox,
};
pub use error::PartError;
pub use feature::{Feature, FeatureLog};
pub use part::Part;
pub use profile::{Plane, Profile};
pub use select::{Dir, EdgeInfo, EdgeQuery, FaceInfo, FaceQuery, SurfaceKind};

/// Everything a modeling script usually needs.
pub mod prelude {
    pub use crate::dressup::{BlendReport, CapKind, Continuity};
    pub use crate::error::PartError;
    pub use crate::part::Part;
    pub use crate::profile::{Plane, Profile};
    pub use crate::select::{Dir, EdgeQuery, FaceQuery, SurfaceKind};
    pub use popkern_core::tolerance::Tolerance;
}
