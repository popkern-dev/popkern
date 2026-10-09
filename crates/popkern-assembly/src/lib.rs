//! `popkern-assembly` — product structure above the scripting recorder.
//!
//! An [`Assembly`] is a document of **part definitions** (each a
//! recorded [`popkern_script::Model`], built once) and **instances**:
//! named, rigidly [`Placement`]-ed occurrences of a definition,
//! nestable through sub-assembly definitions. Geometry lives only in
//! definitions; [`Assembly::flatten`] moves bodies to world
//! coordinates on demand (export as one solid, interference,
//! cross-part booleans). Exports keep the structure: STEP writes a
//! product per definition and an occurrence per instance
//! ([`Assembly::to_step_bytes`]); `.pkn` version 5 carries the tree
//! ([`Assembly::to_pkn_bytes`], [`Document::from_pkn_bytes`]).
//!
//! The design (invariants, wire forms, the popkern live tree, the
//! stage 2–3 families, interference, mates and STEP import) is
//! `docs/architecture/assemblies.md`.
//!
//! ```no_run
//! # use popkern_assembly::{Assembly, Parent, Placement};
//! # fn ball() -> popkern_script::Model { unimplemented!() }
//! let mut asm = Assembly::new("bearing");
//! let ball = asm.add_part("ball", ball())?;          // built once
//! for i in 0..8 {
//!     let a = i as f64 * std::f64::consts::TAU / 8.0;
//!     let at = Placement::translation(0.7 * a.cos(), 0.7 * a.sin(), 0.0);
//!     asm.place(Parent::Root, ball, Some(&format!("ball_{i}")), at)?;
//! }
//! let step = asm.to_step_bytes("bearing")?;          // one solid, eight occurrences
//! # Ok::<(), popkern_assembly::AssemblyError>(())
//! ```

#![forbid(unsafe_code)]

pub mod assembly;
pub mod document;
pub mod document_ops;
pub mod error;
pub mod family;
pub mod flatten;
pub mod hinge;
pub mod hinge_seam;
pub mod history;
pub mod interference;
mod interference_certify;
mod interference_occurrences;
mod interference_refine;
pub mod joints;
mod mate_clusters;
pub mod mate_linalg;
mod mate_predict;
pub mod mate_solve;
mod mate_steps;
pub mod mates;
pub mod mesh_bvh;
mod mesh_cache;
pub mod names;
mod parallel_map;
pub mod params;
mod params_mates;
pub mod part_edit;
mod parts_cache;
pub mod placement;
pub mod step;
pub mod step_import;
pub mod swing;
pub mod tree;
pub mod tri_geometry;
mod trim_side;

pub use assembly::{
    Assembly, Definition, DefinitionContent, DefinitionId, Instance, InstanceRef, Parent, TreeRow,
};
pub use document::{Document, FIRST_ASSEMBLY_VERSION};
pub use error::AssemblyError;
pub use family::{FamilyVariant, ParamValue};
pub use flatten::FlatBody;
pub use hinge::HingeParts;
pub use hinge_seam::SeamHingeParts;
pub use interference::{InterferenceOptions, InterferenceReport, PairReport, PairStatus};
pub use joints::{Joint, JointDrive, RevoluteOptions, RevoluteReport};
pub use mate_solve::{ClusterReport, MATE_TOLERANCE, MateReport};
pub use mates::{Datum, Mate, MateDrive, MateEntity, MateKind, MateRef};
pub use names::{InstancePath, entity_path, validate_name};
pub use params::{AssemblyParam, Drive};
pub use placement::{Placement, RIGID_TOLERANCE};
pub use step_import::StepImportReport;
pub use swing::{SwingOptions, SwingReport, SwingSample};
