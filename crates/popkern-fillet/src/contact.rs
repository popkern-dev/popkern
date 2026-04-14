//! Rolling-ball contact solver (§11.2).
//!
//! Solves the rolling-ball contact problem at one spine point: given two
//! surfaces, a desired ball radius, and a spine plane, find the surface
//! parameters `(u0, v0, u1, v1)` such that the rolling-ball center lies on
//! the inward side of both surfaces at the requested radius and inside the
//! spine plane. The system is closed by the spine-plane constraint into a
//! square 4×4 nonlinear system, solved by damped Newton iteration with
//! analytical Jacobian propagated through the Weingarten map.
//!
//! This module is the foundation of phase 1a of the requirements doc's
//! §11.2 contact problem (`docs/requirements/fillets_blends.md`). It
//! coexists with the existing offset-intersect short-cut in
//! `center_curve.rs` rather than replacing it: analytic-pair fast paths
//! stay; cones, NURBS, and ill-conditioned analytic pairs route through
//! this module.
//!
//! See `docs/superpowers/plans/2026-04-14-fillets-rolling-ball-contact-solver.md`
//! for the full design rationale.
//!
//! ## Module layout
//!
//! Implementation is split across sibling files for the 500-line
//! structural cap. The four pieces are wired here via `#[path]` and
//! their public API is re-exported.
//!
//! - [`types`] — public types: `BallSide`, `ContactRequest`,
//!   `ContactOptions`, `ContactSolution`, `ContactError`, etc.
//! - [`math`] — math primitives: `signed_unit_normal`, Weingarten map,
//!   4×4 Gaussian elimination, condition estimate, Jacobian builder,
//!   Levenberg-Marquardt step.
//! - [`solver`] — the damped-Newton entry point
//!   [`solve_rolling_ball_contact`] and its private helpers.
//! - [`feasible`] — [`pre_check_contact`] and [`solve_with_max_feasible`].

#[path = "contact_feasible.rs"]
pub(crate) mod feasible;
#[path = "contact_math.rs"]
pub(crate) mod math;
#[path = "contact_solver.rs"]
pub(crate) mod solver;
#[path = "contact_types.rs"]
pub(crate) mod types;

pub use feasible::{pre_check_contact, solve_with_max_feasible};
pub use solver::solve_rolling_ball_contact;
pub use types::{
    BallSide, ContactError, ContactFeasibility, ContactOptions, ContactQuality, ContactRequest,
    ContactSolution, ContactTrace, ContactTraceStep, SingularityKind,
};

#[cfg(test)]
#[path = "contact_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "contact_newton_tests.rs"]
mod newton_tests;

#[cfg(test)]
#[path = "contact_jacobian_fd_tests.rs"]
mod jacobian_fd_tests;

#[cfg(test)]
#[path = "contact_robustness_tests.rs"]
mod robustness_tests;

#[cfg(test)]
#[path = "contact_feasible_tests.rs"]
mod feasible_tests;

#[cfg(test)]
#[path = "contact_cone_demo_tests.rs"]
mod cone_demo_tests;
