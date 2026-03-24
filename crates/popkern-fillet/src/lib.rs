#![warn(missing_docs)]
//! popkern-fillet: Fillets and chamfers for the popkern BRep kernel.

pub mod error;
pub mod offset;
pub mod types;

pub use error::{BlendError, BlendFailure, BlendWarning};
pub use types::*;
