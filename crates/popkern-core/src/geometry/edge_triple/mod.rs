//! EdgeTriple factory — "construct, don't fit" entry point.
//!
//! Owns the construction of all three pieces of edge geometry — a 3-D
//! curve plus a pcurve on each adjacent surface — from one solve, so
//! the consistency invariant
//!
//! ```text
//! surf_a(pcurve_a(t)) ≡ curve_3d(t) ≡ surf_b(pcurve_b(t))
//! ```
//!
//! holds by construction. See ADR
//! `docs/decisions/2026-04-25-edge-triple-factory.md` for the full
//! design (hard guarantees G1–G6, dispatch table, multi-branch
//! disambiguation, performance budget, UX affordances).
//!
//! This module ships the public API surface. Stage 1 of the
//! migration: every dispatch path returns
//! [`EdgeTripleError::Unimplemented`]. Subsequent stages fill in
//! analytic, Gregory-boundary, and marching specialisations.

mod analytic;
mod analytic_helpers;
mod analytic_trig;
mod analytic_trig_cone;
mod analytic_trig_cone_cone;
mod analytic_trig_cone_torus_axial;
mod analytic_trig_cyl_cone;
mod analytic_trig_cyl_cone_general;
mod analytic_trig_cyl_cyl;
mod analytic_trig_cyl_cyl_perp;
mod analytic_trig_cyl_sphere;
mod analytic_trig_cyl_torus_axial;
mod analytic_trig_helpers;
mod analytic_trig_plane_cone_hyperbola;
mod analytic_trig_plane_cone_oblique;
mod analytic_trig_plane_cone_parabola;
mod analytic_trig_plane_cyl;
mod analytic_trig_plane_cyl_oblique;
mod analytic_trig_plane_sphere;
mod analytic_trig_plane_sphere_general;
mod analytic_trig_plane_torus;
mod analytic_trig_plane_torus_oblique;
mod analytic_trig_plane_torus_villarceau;
mod analytic_trig_sphere_cone;
mod analytic_trig_sphere_cone_general;
mod analytic_trig_sphere_sphere;
mod analytic_trig_sphere_sphere_general;
mod analytic_trig_sphere_torus_axial;
mod analytic_trig_torus_torus_coaxial;
mod dispatch;
mod error;
mod inputs;
mod levin;
mod options;
mod types;
mod verify;

pub use dispatch::{
    EdgeRequest, EdgeTripleBatchOptions, EdgeTriplePreview, compute_edge_geometry,
    compute_edge_geometry_all_branches, compute_edge_geometry_batch, edge_triple_curvature_comb,
    preview_edge_geometry, revalidate_edge_triple,
};
pub use error::{EdgeTripleError, EdgeTripleVerifyError};
pub use inputs::{BoundarySide, EdgeAnchors, EdgeHint, SeamAxis, SurfaceSide};
pub use options::{BranchPolicy, CancellationToken, EdgeTripleOptions};
pub use types::{
    AnalyticPair, BranchSelector, CurvatureSample, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, GregorySide, ParameterMap,
};

#[cfg(test)]
mod tests;
