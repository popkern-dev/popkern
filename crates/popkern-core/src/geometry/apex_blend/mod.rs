//! Apex-blend termination geometry.
//!
//! For each fillet meeting at a vertex blend, computes the curve
//! where the cap surface naturally terminates the fillet — the
//! "fillet-termination" curve. BRep vertices for the apex region
//! are placed at the termination endpoints, and BRep edges between
//! cap and fillet trace the termination curves.
//!
//! See ADR `docs/decisions/2026-04-27-apex-blend-topology-rebuild.md`
//! for the full design.
//!
//! # Cap-type-specific termination math
//!
//! - **Sphere-fit** (this module's `sphere_fit` submodule):
//!   closed-form algebraic intersection of the cap sphere with
//!   the fillet surface. Implemented for the analytic surface
//!   pairs popkern's cap pipeline produces.
//!
//! - **KP-binonic** (`kp_binonic_term` submodule, future):
//!   iterative spoke-G1 termination for N ∈ {3, 4} Gregory caps.
//!
//! - **KP-polar** (`kp_polar_term` submodule, future): iterative
//!   termination for N ≥ 5 polar K-P caps.
//!
//! # Stage 6c-B0 (spike)
//!
//! The spike-prototype lives in `kp_binonic_term` as a single
//! hand-implemented test case for the symmetric cube corner. Its
//! goal is to prove the spoke-G1 termination math is tractable
//! before committing to the full the outer-fillet-G1 work implementation.

pub mod sphere_fit;

// kp_binonic_term submodule lives in popkern-blend (Stage 6c-B
// requires the cap construction's pub(crate) machinery, which is
// only accessible from inside that crate). The decisive-negative
// apex-termination measurement that motivated this split is
// recorded in docs/decisions/2026-04-27-apex-blend-topology-rebuild.md
// §15.1/§16.2 (the probe module itself was retired with PK-181).
