#![warn(missing_docs)]
//! popkern-sweep: Sweep and loft operations for the popkern BRep kernel.

pub mod arc_length;
pub mod error;
pub mod frame;
pub mod types;

pub use error::{SweepError, SweepWarning};
pub use types::*;
