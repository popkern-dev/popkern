//! Euler operators — the only way to modify BRep topology.
//!
//! Each operator maintains the Euler-Poincaré invariant by construction.

mod mev;
mod mvfs;
pub use mev::MevResult;
pub use mvfs::MvfsResult;
