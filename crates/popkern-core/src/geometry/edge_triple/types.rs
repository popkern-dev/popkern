//! Core data types: [`EdgeTriple`], provenance, quality, topology.

use serde::{Deserialize, Serialize};

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};

/// Whether the edge is open (start ≠ end) or closed (loop).
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum EdgeTopology {
    /// `curve_3d.point_at(t_start)` is the start anchor;
    /// `curve_3d.point_at(t_end)` is the end anchor.
    Open,
    /// Closed loop (cylinder seam, sphere equator). The 3-D curve's
    /// endpoints coincide; pcurves are periodic in the matching
    /// surface dimension.
    Closed {
        /// `t_end - t_start`.
        period: f64,
    },
}

/// How the factory built the triple. Diagnostics and UI read this
/// directly to label edges by class without re-inferring from the
/// surface pair.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum EdgeTripleProvenance {
    /// Closed-form intersection of two analytic surfaces.
    AnalyticIntersection {
        /// Discrete tag identifying the analytic surface-pair class.
        kind: AnalyticPair,
        /// Which intersection branch was picked.
        branch: BranchSelector,
    },
    /// Built from a known canonical surface boundary (Gregory cap
    /// boundary line in `[0,1]²`).
    GregoryBoundary {
        /// Which surface(s) the Gregory boundary is anchored to.
        side: GregorySide,
        /// Which of the N boundaries: 0..N for an N-sided patch.
        boundary_index: u8,
    },
    /// Solved by the marching solver.
    Marching {
        /// Achieved per-sample positional residual.
        residual: f64,
        /// Iterations consumed before convergence or abort.
        iterations: u32,
        /// True if the solver converged within tolerance.
        converged: bool,
    },
    /// Self-loop edge: surf_a == surf_b on a degenerate seam.
    DegenerateSeam {
        /// The fixed UV coordinate value of the seam (e.g. `u=0`
        /// for a cylinder seam in `(u, v)` coordinates).
        seam_param: f64,
    },
}

/// Discrete tag for analytic surface-pair classes. Members are
/// alphabetised so the canonical ordering in the dispatch is stable.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalyticPair {
    /// Plane × plane — line, or empty when parallel.
    PlanePlane,
    /// Plane × cylinder — ellipse, two lines, or one line (tangent).
    PlaneCylinder,
    /// Plane × sphere — circle, point (tangent), or empty.
    PlaneSphere,
    /// Plane × cone — conic section.
    PlaneCone,
    /// Plane × torus — Villarceau circles or general quartic.
    PlaneTorus,
    /// Cylinder × cylinder — ellipses (parallel/coplanar) or quartic.
    CylinderCylinder,
    /// Cylinder × sphere — Viviani's curve or two circles.
    CylinderSphere,
    /// Cylinder × cone — closed-form when axes coplanar.
    CylinderCone,
    /// Cylinder × torus — closed-form when cylinder axis aligned.
    CylinderTorus,
    /// Sphere × sphere — circle.
    SphereSphere,
    /// Sphere × cone — closed-form when sphere centre on cone axis.
    SphereCone,
    /// Sphere × torus — closed-form when sphere centre on axis.
    SphereTorus,
    /// Cone × cone — closed-form when axes coplanar.
    ConeCone,
}

/// Which intersection branch the disambiguation logic picked. The
/// branch space is documented per pair class in the dispatch
/// implementation.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BranchSelector {
    /// The pair has exactly one branch (sphere×sphere → 1 circle,
    /// plane×sphere → 1 circle, …).
    OnlyBranch,
    /// One of multiple branches; index in the canonical enumeration.
    Indexed(u32),
}

/// On a Gregory-boundary triple, which surface(s) the construction
/// is anchored to.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GregorySide {
    /// Cap is `surf_a` (so `pcurve_a` is the canonical UV boundary).
    A,
    /// Cap is `surf_b`.
    B,
    /// Both surfaces are Gregory caps (corner-of-corner).
    Both,
}

/// Achieved quality of a successfully-constructed triple. Always
/// populated, even on analytic dispatch paths where every field
/// reads zero or near-ulp.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct EdgeTripleQuality {
    /// Maximum positional residual over the verification samples
    /// (G2): `‖surf(pcurve(t)) − curve_3d(t)‖`, max over both sides.
    pub max_positional_residual: f64,
    /// `‖curve_3d(t_start) − anchors.start‖` + likewise at end. Zero
    /// on analytic dispatch paths that hit endpoints exactly via CP
    /// equality (G1).
    pub max_endpoint_error: f64,
    /// Number of verification samples taken.
    pub samples: u32,
    /// `tol.linear()` at construction time. Snapshot tests pin this
    /// so a tolerance-system regression surfaces here first.
    pub tolerance_used: f64,
}

/// The triple itself. Returned by
/// [`crate::geometry::edge_triple::compute_edge_geometry`].
///
/// Hard guarantees on every successful return are listed in the
/// module-level docstring and in the ADR §1 (G1–G6).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeTriple {
    /// 3-D curve on the edge.
    pub curve_3d: NurbsCurve,
    /// Pcurve on `surf_a` — UV stored as `Point3(u, v, 0.0)`.
    pub pcurve_a: NurbsCurve,
    /// Pcurve on `surf_b`.
    pub pcurve_b: NurbsCurve,
    /// Shared parameter range start (G3).
    pub t_start: f64,
    /// Shared parameter range end (G3); `t_start ≤ t_end` (G5).
    pub t_end: f64,
    /// Open or closed-loop topology.
    pub topology: EdgeTopology,
    /// How the factory built this triple — read by diagnostics.
    pub provenance: EdgeTripleProvenance,
    /// Achieved residual / endpoint error from construction.
    pub quality: EdgeTripleQuality,
}

/// One sample on the curvature comb of a 3-D edge curve. UI renders
/// the result as a ribbon of normals scaled by curvature magnitude.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CurvatureSample {
    /// Curve parameter value at the sample.
    pub t: f64,
    /// Position on the 3-D curve at `t`.
    pub position: Point3,
    /// Unit tangent at `t`.
    pub tangent: Vec3,
    /// Magnitude of the curvature vector.
    pub curvature: f64,
    /// Direction of the principal normal (curvature vector,
    /// normalised).
    pub normal: Vec3,
}
