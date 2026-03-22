//! Euler operators — the only way to modify BRep topology.
//!
//! Each operator maintains the Euler-Poincaré invariant by construction.

mod kemr;
mod kfmrh;
mod mef;
mod mev;
mod mvfs;
mod sweep;
pub use kemr::{KemrResult, MekrResult};
pub use kfmrh::MfkrhResult;
pub use mef::MefResult;
pub use mev::MevResult;
pub use mvfs::MvfsResult;
pub use sweep::ExtrudeResult;
pub use sweep::RevolveResult;
