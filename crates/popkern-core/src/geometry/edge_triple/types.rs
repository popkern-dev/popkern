//! Core data types: [`EdgeTriple`], provenance, quality, topology.

use serde::{Deserialize, Serialize};

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::inputs::SurfaceSide;
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
        /// Cardinality report: total branches in the analytic pair
        /// and the conic kind of each. `None` when the dispatcher
        /// doesn't compute this (most specialised paths produce a
        /// single branch and leave it implicit). Levin's quadric-pair
        /// dispatcher fills it with `Some` so UI / boolean ops can
        /// see the structure of the intersection.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        cardinality: Option<BranchCardinality>,
    },
    /// Built from a known canonical surface boundary (Gregory cap
    /// boundary line in `[0,1]²`).
    GregoryBoundary {
        /// Which surface(s) the Gregory boundary is anchored to.
        side: GregorySide,
        /// Which of the N boundaries: 0..N for an N-sided patch.
        boundary_index: u8,
    },
    /// Built from a known canonical bounded-surface boundary.
    SurfaceBoundary {
        /// Which surface owns the canonical boundary curve.
        side: SurfaceSide,
        /// Which boundary was used, encoded with the same order as [`BoundarySide`].
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

/// Conic / curve type of one branch of an analytic intersection.
/// Aligned with `CurveGeometry`'s analytic variants.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BranchKind {
    /// Straight line (e.g., plane × plane).
    Line,
    /// Full circle (e.g., plane × sphere, plane × cyl perpendicular).
    Circle,
    /// Circular arc (closed-form circle that's a sub-range of `[0, 2π)`).
    Arc,
    /// Ellipse (plane × cyl oblique, plane × cone elliptic section).
    Ellipse,
    /// Parabola (plane × cone parallel-to-slant).
    Parabola,
    /// Hyperbola (plane × cone oblique past the slant, e.g., cone × cone
    /// non-coaxial Levin).
    Hyperbola,
    /// Generic NURBS curve (sample-and-fit fallback).
    Nurbs,
}

impl BranchKind {
    /// Classify a `CurveGeometry` into a `BranchKind`. The mapping
    /// is fixed per analytic curve variant.
    pub fn from_curve(c: &crate::geometry::curve::CurveGeometry) -> Self {
        use crate::geometry::curve::CurveGeometry::*;
        match c {
            Line { .. } => BranchKind::Line,
            Circle { .. } => BranchKind::Circle,
            Arc { .. } => BranchKind::Arc,
            Ellipse { .. } => BranchKind::Ellipse,
            Parabola { .. } => BranchKind::Parabola,
            Hyperbola { .. } => BranchKind::Hyperbola,
            Nurbs(_) => BranchKind::Nurbs,
        }
    }
}

/// Cardinality report for an analytic intersection: how many real
/// branches the surface pair admits, and the conic kind of each.
///
/// Levin's algorithm computes this from the rank-decomposition of
/// the pencil's degenerate fibre (rank-1 ⇒ 1 plane, rank-2 mixed-sign
/// ⇒ 2 planes); each plane × surface gives a known conic kind.
/// Specialised single-branch dispatchers leave `cardinality = None`
/// on the provenance — `total = 1, kinds = [the_branch_kind]` is
/// implicit.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BranchCardinality {
    /// Total number of real intersection branches.
    pub total: u8,
    /// Conic kind per branch, in canonical (dispatch) order.
    pub kinds: Vec<BranchKind>,
    /// True iff the surfaces meet *tangentially* along this curve —
    /// i.e., the two surfaces share a tangent plane at every point of
    /// the intersection (multiplicity-2 contact). Detected exactly via
    /// the symbolic Q rational-pencil rank-1 fibre; downstream
    /// boolean / mesh / IO consumers branch on this when handling
    /// kissing edges, double-curves, or G1 stitching across the join.
    /// Defaults to `false` for transverse intersections.
    #[serde(default)]
    pub tangent: bool,
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
    /// Cone × torus — closed-form when cone apex on torus axis and
    /// axes parallel.
    ConeTorus,
    /// Torus × torus — closed-form when axes coincident.
    TorusTorus,
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
/// module-level docstring and in the ADR §1 (G1–G6) with §13's
/// revised G3 (parameter-map correspondence rather than literal
/// shared knots).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct EdgeTriple {
    /// 3-D curve on the edge. Stored as a `CurveGeometry` enum so
    /// analytic primitives (Line, Circle, Arc) can be carried
    /// exactly — matching commercial-kernel and OpenCASCADE practice.
    pub curve_3d: CurveGeometry,
    /// Pcurve on `surf_a`.
    pub pcurve_a: crate::geometry::pcurve::Pcurve,
    /// Pcurve on `surf_b`.
    pub pcurve_b: crate::geometry::pcurve::Pcurve,
    /// 3-D curve parameter range start.
    pub t_start: f64,
    /// 3-D curve parameter range end; `t_start ≤ t_end` (G5).
    pub t_end: f64,
    /// Open or closed-loop topology.
    pub topology: EdgeTopology,
    /// Bijective map from `curve_3d`'s parameter to each pcurve's
    /// parameter (ADR §13 G3). For analytic dispatches this is
    /// `Identity`; for NURBS marching it's `SharedKnots`.
    pub parameter_map: ParameterMap,
    /// How the factory built this triple — read by diagnostics.
    pub provenance: EdgeTripleProvenance,
    /// Achieved residual / endpoint error from construction.
    pub quality: EdgeTripleQuality,
}

/// Bijective correspondence between the 3-D curve's parameter and
/// each pcurve's parameter. Used by G2 verification and by
/// downstream consumers that need to evaluate `surf(pcurve(t))` at
/// the same point as `curve_3d(t)`.
///
/// See ADR §13.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ParameterMap {
    /// `Φ_a(t) = Φ_b(t) = t`. The 3-D curve and both pcurves share
    /// the natural parameter (e.g. circle on cylinder where the
    /// cyl's u-axis IS the circle's angle).
    Identity,
    /// `Φ_x(t) = a_x + b_x · t` for each side. Linear
    /// reparameterisation per side.
    LinearAffine {
        /// `Φ_a(t) = a_a + b_a · t`.
        a_a: f64,
        /// `Φ_a(t) = a_a + b_a · t`.
        b_a: f64,
        /// `Φ_b(t) = a_b + b_b · t`.
        a_b: f64,
        /// `Φ_b(t) = a_b + b_b · t`.
        b_b: f64,
    },
    /// All three curves share knot vector and parameter range —
    /// the stage-1 invariant. Used for NURBS×NURBS marching
    /// where there's no analytic correspondence.
    SharedKnots,
}

impl ParameterMap {
    /// Apply the map to obtain `(t_pcurve_a, t_pcurve_b)` from a
    /// 3-D parameter `t`.
    pub fn apply(&self, t: f64) -> (f64, f64) {
        match self {
            ParameterMap::Identity | ParameterMap::SharedKnots => (t, t),
            ParameterMap::LinearAffine { a_a, b_a, a_b, b_b } => (a_a + b_a * t, a_b + b_b * t),
        }
    }
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
