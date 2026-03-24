#![warn(missing_docs)]
//! popkern-fillet: Fillets and chamfers for the popkern BRep kernel.

pub mod center_curve;
pub mod error;
pub mod fillet_surface;
pub mod offset;
pub mod trim;
pub mod types;

pub use error::{BlendError, BlendFailure, BlendWarning};
pub use types::*;
