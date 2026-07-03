//! Public types for the rolling-ball contact solver. Split from
//! `contact.rs` to keep individual files under the 500-line structural
//! cap; re-exported through `crate::contact`.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use thiserror::Error;

/// Which side of a surface the rolling-ball *center* sits on, defined
/// **relative to the surface's natural `normal_at(u, v)` direction**.
///
/// `Outside` ⇒ ball center on the `+normal_at` side.
/// `Inside`  ⇒ ball center on the `-normal_at` side.
///
/// **The caller is responsible for translating their notion of "outside
/// the solid" into the right `BallSide` value, because the natural
/// normal direction depends on the surface parameterization.** For
/// example, a `Sphere` parameterized by `(longitude, colatitude)` has
/// a natural normal pointing *inward* (`partial_u × partial_v` is
/// inward at the equator); a ball *outside* the sphere therefore
/// requires `BallSide::Inside`.
///
/// In a popkern-fillet context, callers typically already have a
/// `face.orientation: bool` flag that records "does this face's
/// natural surface normal agree with the outward-from-solid direction?"
/// — translate as `BallSide::Outside` if `orientation` agrees with the
/// rolling ball's side, `BallSide::Inside` otherwise.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BallSide {
    /// Ball center on the `+normal_at(u, v)` side of the surface.
    Outside,
    /// Ball center on the `-normal_at(u, v)` side of the surface.
    Inside,
}

/// Translate popkern face orientations + edge convexity into the
/// contact solver's [`BallSide`] convention.
///
/// popkern face normals point INTO the solid, and a face's natural
/// surface normal equals that face normal when `orientation = true`
/// (flipped when `false`). Which side the generating ball sits on
/// depends on the EDGE CLASS, not the face orientation alone:
///
/// - **Convex edge** (material wedge < 180°): the ball is inscribed
///   in the material — center on the into-solid side of both faces
///   (the classic dihedral fillet cylinder axis). There is no
///   exterior ball tangent to both material faces, only to their
///   extensions; converging to that mirror solution builds a
///   "weld bead" blend with contacts on the non-material plane
///   extensions (PK-091).
/// - **Concave edge** (reentrant): the ball rolls in the void —
///   center opposite the into-solid normal.
///
/// `Mixed` edges (tangent or near-tangent surfaces) map like
/// concave: tangent-continuation edges arise in concave grooves —
/// a secondary blend rolling between a primary fillet and its
/// support face — where the ball sits in the void. This also
/// preserves the behaviour the blend-on-blend acceptance gates were
/// built on.
pub(crate) fn ball_sides_for_edge(
    convexity: crate::analysis::Convexity,
    orientation_1: bool,
    orientation_2: bool,
) -> (BallSide, BallSide) {
    use crate::analysis::Convexity;
    let material = |orientation: bool| {
        if orientation {
            BallSide::Outside
        } else {
            BallSide::Inside
        }
    };
    let void = |orientation: bool| {
        if orientation {
            BallSide::Inside
        } else {
            BallSide::Outside
        }
    };
    match convexity {
        Convexity::Concave => (void(orientation_1), void(orientation_2)),
        Convexity::Convex | Convexity::Mixed => (material(orientation_1), material(orientation_2)),
    }
}

impl BallSide {
    /// Sign multiplier on the surface normal vector. `Outside → +1`,
    /// `Inside → -1`. The ball center is `S(u, v) + sign · r · n̂`.
    pub(crate) fn sign(self) -> f64 {
        match self {
            BallSide::Outside => 1.0,
            BallSide::Inside => -1.0,
        }
    }
}

/// Inputs to a single contact-solver call.
///
/// The contact problem at one spine point is a 3-equation contact
/// equality plus a spine-plane constraint, yielding a 4×4 nonlinear
/// system in `(u0, v0, u1, v1)` solved by damped Newton iteration
/// (§11.2). The spine-plane constraint pins the otherwise
/// under-determined system: the rolling-ball center is required to
/// lie in the plane through `spine_point` perpendicular to
/// `spine_plane_normal`.
///
/// Walkers pass the spine tangent as `spine_plane_normal`. One-shot
/// callers pass the edge tangent at the parameter of interest.
pub struct ContactRequest<'a> {
    /// First contact surface (any type implementing [`SurfaceEval`]).
    pub surface_0: &'a dyn SurfaceEval,
    /// Second contact surface.
    pub surface_1: &'a dyn SurfaceEval,
    /// Which side of `surface_0` the ball center sits on.
    pub side_0: BallSide,
    /// Which side of `surface_1` the ball center sits on.
    pub side_1: BallSide,
    /// Rolling-ball radius. Must be > 0; convexity is encoded via `side_*`.
    pub radius: f64,
    /// Initial guess for `surface_0`'s contact parameter.
    pub initial_uv0: (f64, f64),
    /// Valid range `((u_min, u_max), (v_min, v_max))` for `surface_0`.
    /// Newton steps are clamped to this box; drifting outside returns
    /// [`ContactError::OutOfBounds`].
    pub uv0_range: ((f64, f64), (f64, f64)),
    /// Initial guess for `surface_1`'s contact parameter.
    pub initial_uv1: (f64, f64),
    /// Valid range for `surface_1`'s `(u, v)`.
    pub uv1_range: ((f64, f64), (f64, f64)),
    /// 3D point on the spine. The rolling-ball center is constrained
    /// to lie in the plane through this point with normal
    /// `spine_plane_normal`.
    pub spine_point: Point3,
    /// Normal vector defining the spine plane. Walkers pass the spine
    /// tangent here; one-shot callers pass the edge tangent at the
    /// parameter of interest. Need not be unit length.
    pub spine_plane_normal: Vec3,
}

/// Tunable knobs for the contact solver. Use [`ContactOptions::default`]
/// for the commercial-grade preset.
#[derive(Debug, Clone, Copy)]
pub struct ContactOptions {
    /// Maximum Newton iterations before giving up. Default 20.
    pub max_iterations: usize,
    /// Residual norm tolerance (3D distance, model units). Default
    /// `1e-6`; setting it here overrides whatever the call-site
    /// `Tolerance` would otherwise imply.
    pub residual_tolerance: f64,
    /// Initial Levenberg-Marquardt damping factor. `0.0` is pure
    /// Newton (default); positive values regularize from the first
    /// step.
    pub initial_damping: f64,
    /// Quality preset. `Final` runs to full tolerance; `Preview`
    /// stops early (3 iterations, looser residual) for sub-100ms
    /// slider response in interactive UX.
    pub quality: ContactQuality,
    /// When true, clamp out-of-bounds initial UV guesses to the
    /// range boundary instead of rejecting with `OutOfBounds`.
    /// Used for §12.3 blend-on-blend where torus UV projection
    /// can produce wrapped parameters.
    pub clamp_initial_uv: bool,
}

impl Default for ContactOptions {
    fn default() -> Self {
        Self {
            max_iterations: 20,
            residual_tolerance: 1e-6,
            initial_damping: 0.0,
            quality: ContactQuality::Final,
            clamp_initial_uv: false,
        }
    }
}

/// Cost vs accuracy preset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactQuality {
    /// 3 iterations, 10× looser residual. Use during a UI slider drag
    /// where sub-100ms response matters more than full precision.
    Preview,
    /// Full tolerance, up to `max_iterations` Newton steps. Use for
    /// the final commit where precision is essential.
    Final,
}

/// A converged contact solution. The fields satisfy:
/// `center == surface_0.point_at(uv0) + side_0·radius·n̂0(uv0)`
/// `        == surface_1.point_at(uv1) + side_1·radius·n̂1(uv1)`
/// to within `residual`.
#[derive(Debug, Clone)]
pub struct ContactSolution {
    /// Contact parameters on surface 0.
    pub uv0: (f64, f64),
    /// Contact parameters on surface 1.
    pub uv1: (f64, f64),
    /// Rolling-ball center in 3D space.
    pub center: Point3,
    /// 3D contact point on surface 0 (`= surface_0.point_at(uv0)`).
    pub contact0: Point3,
    /// 3D contact point on surface 1.
    pub contact1: Point3,
    /// Number of Newton iterations consumed.
    pub iterations: usize,
    /// §18 telemetry: total backtracking-line-search rejections
    /// across all Newton iterations. Each rejected step is a
    /// candidate `α` value whose residual failed to improve on
    /// the current state before the line search shrunk `α`.
    /// A healthy convergent solve has very few rejections (often
    /// zero); a high rejection count per iteration indicates
    /// ill-conditioning or a near-tangent surface pair.
    pub step_rejections: usize,
    /// Final residual norm (3D, model units).
    pub residual: f64,
    /// Estimated condition number of the final 4×4 Jacobian (1-norm).
    pub condition: f64,
    /// The final 4×4 Jacobian. Reused by:
    /// - the spine walker (phase 1b) for continuation step prediction;
    /// - G2 cross-section construction (phase 1d) for curvature matching;
    /// - sensitivity analysis (§14.8).
    pub jacobian: [[f64; 4]; 4],
}

/// Optional per-iteration trace for debugging and telemetry. Pass
/// `Some(&mut ContactTrace)` to `solve_rolling_ball_contact` to record
/// every Newton step's state. Default `None` → zero overhead.
#[derive(Debug, Clone, Default)]
pub struct ContactTrace {
    /// One entry per Newton iteration (in iteration order).
    pub steps: Vec<ContactTraceStep>,
}

/// One Newton step's state, captured in [`ContactTrace`].
#[derive(Debug, Clone, Copy)]
pub struct ContactTraceStep {
    /// 0-based iteration index.
    pub iteration: usize,
    /// Residual norm at the start of this step.
    pub residual: f64,
    /// Norm of the Newton step `||Δx||`.
    pub step_length: f64,
    /// Estimated condition of the Jacobian for this step.
    pub condition: f64,
    /// Levenberg-Marquardt damping factor used (0.0 if pure Newton).
    pub damping: f64,
    /// Surface 0 parameter at the start of this step.
    pub uv0: (f64, f64),
    /// Surface 1 parameter at the start of this step.
    pub uv1: (f64, f64),
}

/// Why a [`ContactError::SingularJacobian`] occurred. Lets a GUI
/// highlight the specific surface (or pair) responsible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SingularityKind {
    /// `surface_0`'s normal vector has length below the numerical floor.
    Normal0Degenerate,
    /// `surface_1`'s normal vector has length below the numerical floor.
    Normal1Degenerate,
    /// The two surface normals are nearly parallel; no rolling ball
    /// fits between them.
    NormalsParallel,
    /// The 4×4 Jacobian is rank-deficient for an unspecified reason
    /// (typically a parameterization singularity such as a sphere
    /// pole or cone apex).
    JacobianRankDeficient,
}

/// Failure modes of [`super::solve_rolling_ball_contact`]. Each variant
/// carries enough structured data for `BlendFailure::ContactSolverFailed`
/// to be populated with a 3D location and remediation suggestion.
#[derive(Debug, Clone, Error)]
pub enum ContactError {
    /// Newton residual increased over multiple consecutive steps.
    #[error(
        "rolling-ball contact solver diverged after {iterations} iterations \
         (residual {last_residual:.3e})"
    )]
    Diverged {
        /// Iterations attempted before declaring divergence.
        iterations: usize,
        /// Residual at the last iteration.
        last_residual: f64,
        /// Surface 0 parameter at the last iteration.
        last_uv0: (f64, f64),
        /// Surface 1 parameter at the last iteration.
        last_uv1: (f64, f64),
    },
    /// Newton did not reach `residual_tolerance` within `max_iterations`.
    #[error(
        "rolling-ball contact solver did not converge \
         (residual {residual:.3e} after {iterations} iterations)"
    )]
    NotConverged {
        /// Iterations attempted (always `max_iterations`).
        iterations: usize,
        /// Final residual.
        residual: f64,
        /// Surface 0 parameter at the last iteration.
        last_uv0: (f64, f64),
        /// Surface 1 parameter at the last iteration.
        last_uv1: (f64, f64),
    },
    /// Jacobian became near-singular and Levenberg-Marquardt
    /// regularization could not recover.
    #[error(
        "rolling-ball contact solver hit a singular Jacobian \
         ({kind:?}, condition {condition:.3e}) after {iterations} iterations"
    )]
    SingularJacobian {
        /// Iterations consumed.
        iterations: usize,
        /// Estimated condition number when LM gave up.
        condition: f64,
        /// 3D location of the failure (surface midpoint at the last
        /// iteration).
        location: Point3,
        /// What kind of mathematical singularity caused the failure.
        kind: SingularityKind,
    },
    /// A surface normal evaluation produced a degenerate vector that
    /// could not be repaired by perturbation.
    #[error(
        "surface {surface_index} normal is degenerate at parameter \
         ({:.4}, {:.4})",
        uv.0,
        uv.1
    )]
    DegenerateNormal {
        /// `0` for surface 0, `1` for surface 1.
        surface_index: usize,
        /// Parameter at which the degeneracy was found.
        uv: (f64, f64),
        /// 3D location of the degenerate point.
        location: Point3,
    },
    /// The solution drifted outside the supplied UV bounds and
    /// clamping pushed it back to the boundary, producing an
    /// infeasible state.
    #[error(
        "rolling-ball contact drifted off surface {surface_index} \
         at parameter ({:.4}, {:.4})",
        uv.0,
        uv.1
    )]
    OutOfBounds {
        /// `0` for surface 0, `1` for surface 1.
        surface_index: usize,
        /// Pre-clamp parameter (where the iteration wanted to go).
        uv: (f64, f64),
        /// Post-clamp parameter (the actual returned point).
        clamped_uv: (f64, f64),
    },
}

/// Result of a [`super::pre_check_contact`] call. Indicates whether
/// the caller should proceed with [`super::solve_rolling_ball_contact`]
/// or expect failure.
#[derive(Debug, Clone)]
pub enum ContactFeasibility {
    /// The configuration looks well-conditioned. Proceed.
    LikelyFine,
    /// The configuration may converge but expect tightness near a
    /// face boundary or numerical conditioning. Reason explains.
    LikelyTight {
        /// Human-readable reason this configuration looks tight.
        reason: String,
    },
    /// The configuration is essentially hopeless: no rolling-ball
    /// fillet exists at the requested radius. The caller should
    /// not waste cycles in [`super::solve_rolling_ball_contact`].
    LikelyInfeasible {
        /// Human-readable reason no fillet exists.
        reason: String,
        /// Suggested maximum feasible radius, if computable.
        max_feasible_radius: Option<f64>,
    },
}
