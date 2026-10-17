//! API types for fillet and chamfer specifications.

use serde::{Deserialize, Serialize};

#[path = "types_chamfer.rs"]
mod types_chamfer;
#[path = "types_corner.rs"]
mod types_corner;
#[path = "types_exact_boundary.rs"]
mod types_exact_boundary;
#[path = "types_repair.rs"]
mod types_repair;
#[path = "types_repair_blend_seam.rs"]
mod types_repair_blend_seam;
#[path = "types_repair_exact_boundary.rs"]
mod types_repair_exact_boundary;
#[path = "types_repair_setback_fit.rs"]
mod types_repair_setback_fit;
#[path = "types_repair_terminal_arc.rs"]
mod types_repair_terminal_arc;

pub use types_chamfer::{ChamferKind, ChamferOptions, ChamferSpec};
pub use types_corner::{
    CornerMethod, CornerPatchInfo, CornerSubPatchFamily, CornerSubPatchProvenance,
    CornerSubPatchRole,
};
pub use types_exact_boundary::{
    CapContinuityReport, CapContinuityTarget, ExactBoundaryCandidateOutcome,
    ExactBoundaryCapReport, ExactBoundaryReportDiagnostic, ExactBoundaryReportProvenance,
    ExactBoundaryReportSector, ExactBoundaryReportTransactionStatus,
    ExactBoundaryRouteDecisionOutcome, ExactBoundaryRouteDecisionStep,
    ExactBoundaryToleranceEvidence, INFORMATIONAL_DIAGNOSTIC_CODES,
    diagnostic_code_is_informational,
};
pub use types_repair::{
    ApexConstructionDiagnostic, ApexConstructionRoute, ApexConstructionRouteDecision,
    BlendRepairReport, CapSupportContinuityDiagnostic, CoherenceValidationDiagnostic,
    CylinderTrimFarSideStrategy, CylinderTrimMaterializationDeclineReason,
    CylinderTrimMaterializationDiagnostic, ExactBoundaryCandidatePerformanceReport,
    SERVING_SHAPE_DECLINE_PREFIX, ServingSelfIntersectionRetry, SideSupportContinuationDiagnostic,
    SupportCellSplitAxisDiagnostic, SupportCellSplitAxisStatus, SupportCellSplitDiagnostic,
    SupportCellSplitFallbackReason, SupportSeamDecision, SupportSeamDeclineReason,
    SupportSeamRelocationReport, SupportSeamRollbackReason, build_blend_repair_action_artifact,
    support_seam_decline_reason_label, support_seam_rollback_reason_label,
};
pub use types_repair_blend_seam::*;
pub use types_repair_exact_boundary::ExactBoundaryCapDiagnostic;
pub use types_repair_setback_fit::*;
pub use types_repair_terminal_arc::{
    TerminalArcTrimDecline, TerminalArcTrimDeclineReason, TerminalArcTrimReport,
};

use popkern_core::topology::EdgeHandle;

use crate::error::BlendWarning;
pub use crate::types_blend_options::{
    BlendProfile, CapKind, CornerCapPreference, CornerConstruction, CornerSetback, EdgeEndSetback,
    InterferencePolicy,
};

/// Result of a successful blend operation.
#[derive(Debug, Clone)]
pub struct BlendResult {
    /// The modified BRep with blends applied.
    pub brep: popkern_core::topology::BRep,
    /// Per-edge warnings emitted during blending.
    pub edge_warnings: Vec<(EdgeHandle, BlendWarning)>,
    /// §5.3 per-corner patch reports — one entry per vertex where
    /// three or more filleted edges met and a corner patch was
    /// inserted. Empty when the operation produced no vertex
    /// blends. See [`CornerPatchInfo`] for the schema.
    pub corner_patches: Vec<CornerPatchInfo>,
    /// Structured diagnostics for repair and hardening actions performed
    /// during blend construction.
    pub repair_report: BlendRepairReport,
}

/// Specification for a single fillet edge.
#[derive(Debug, Clone)]
pub struct FilletSpec {
    /// The edge to fillet.
    pub edge: EdgeHandle,
    /// The fillet radius profile.
    pub radius: FilletRadius,
}

impl FilletSpec {
    /// Copies of `specs` with every radius law canonicalized (see
    /// [`FilletRadius::canonical`]): a constant-valued
    /// Variable/ByPoints law is a constant law and must take the
    /// analytic cylinder pipeline — the swept-NURBS spelling of the
    /// same geometry shipped open pair corners because every corner
    /// settlement rightly recognizes cylinder strips only.
    pub fn canonicalized_all(specs: &[FilletSpec]) -> Vec<FilletSpec> {
        specs
            .iter()
            .map(|spec| FilletSpec {
                edge: spec.edge,
                radius: spec.radius.canonical(),
            })
            .collect()
    }
}

/// Fillet radius profile along an edge.
#[derive(Debug, Clone)]
pub enum FilletRadius {
    /// Constant radius along the entire edge.
    Constant(f64),
    /// Linearly varying radius from start to end of the edge.
    Variable {
        /// Radius at the start vertex.
        start: f64,
        /// Radius at the end vertex.
        end: f64,
    },
    /// Radius defined at specific parameter positions along the edge.
    ByPoints(Vec<(f64, f64)>),
}

impl FilletRadius {
    /// Canonical form of the law: a profile whose value is constant
    /// along the edge IS a constant law. `Variable` with equal ends
    /// and `ByPoints` with a single value collapse to `Constant`, so
    /// identical geometry always takes the identical pipeline (the
    /// analytic cylinder strip and its corner settlements) regardless
    /// of how the law was spelled. Exact value equality only — a law
    /// with any real variation is preserved bit-for-bit, so this
    /// carries no tolerance semantics.
    pub fn canonical(&self) -> FilletRadius {
        match self {
            FilletRadius::Constant(r) => FilletRadius::Constant(*r),
            FilletRadius::Variable { start, end } if start == end => FilletRadius::Constant(*start),
            FilletRadius::ByPoints(points)
                if !points.is_empty() && points.iter().all(|(_, v)| *v == points[0].1) =>
            {
                FilletRadius::Constant(points[0].1)
            }
            other => other.clone(),
        }
    }
}

impl From<f64> for FilletRadius {
    /// Convenience: a bare `f64` becomes a constant-radius profile.
    /// Lets `modify_fillet_radius(brep, face, 0.2, tol)` and
    /// `modify_fillet_radius(brep, face, FilletRadius::Variable { ... }, tol)`
    /// share the same entry point.
    fn from(r: f64) -> Self {
        FilletRadius::Constant(r)
    }
}

/// Options controlling fillet behaviour.
#[derive(Debug, Clone)]
pub struct FilletOptions {
    /// Whether to automatically propagate fillets along tangent-continuous edge chains.
    pub propagate_tangent: bool,
    /// How to handle vertex blends where three or more fillets meet.
    pub vertex_mode: VertexMode,
    /// Which construction closes setback-eligible corners, and whether
    /// it may fall back (PK-603). Default: the kernel's ladder.
    pub corner_cap: CornerCapPreference,
    /// User-set setback distances at the corners the setback
    /// construction serves (PK-623). Default: the kernel's automatic
    /// setbacks.
    pub corner_setback: CornerSetback,
    /// Requested blend-boundary continuity level. The fillet engine
    /// tries to build a surface that meets this continuity with the
    /// two support faces at every contact point; falls back to a
    /// lower level and emits [`BlendWarning::ContinuityDowngraded`]
    /// if the requested level isn't achievable for the surface pair.
    ///
    /// - [`ContinuityLevel::G0`] — positional only. Debug / special
    ///   cases. Currently produces the same output as `G1` since the
    ///   pipeline has no strict-G0 fast path.
    /// - [`ContinuityLevel::G1`] (default) — tangent continuity. The
    ///   circular-arc rolling-ball sweep; sufficient for
    ///   machining and visual smoothness on non-Class-A surfaces.
    /// - [`ContinuityLevel::G2`] — curvature continuity. Routes
    ///   through the spine walker and skins degree-5 Bezier
    ///   cross-sections. Required for §15 Class-A surfaces (no
    ///   highlight-line break at the blend boundary).
    ///
    /// **Cost.** G2 is materially more expensive than G1 — the
    /// walker performs ~20 contact-solver iterations per spine
    /// sample where G1 does a single analytic offset intersect.
    /// Leave this `G1` for interactive preview; enable `G2` for
    /// final / Class-A output.
    ///
    /// **Fallback.** If the walker fails or the G2 sweep errors,
    /// the pipeline emits a [`BlendWarning::ContinuityDowngraded`]
    /// with the reason and produces the G1 result. Consult
    /// [`BlendQuality::achieved_continuity`](crate::quality::BlendQuality::achieved_continuity)
    /// to see what was actually delivered.
    ///
    /// See §2.4 / §7 G2 in `docs/requirements/fillets_blends.md`.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use popkern_blend::{ContinuityLevel, FilletOptions, FilletSpec, FilletRadius, fillet_edges};
    /// # use popkern_core::topology::BRep;
    /// # use popkern_core::tolerance::Tolerance;
    /// # let brep: BRep = unimplemented!();
    /// # let edge = unimplemented!();
    /// let specs = vec![FilletSpec {
    ///     edge,
    ///     radius: FilletRadius::Constant(0.1),
    /// }];
    /// let opts = FilletOptions {
    ///     continuity: ContinuityLevel::G2,
    ///     ..FilletOptions::default()
    /// };
    /// let result = fillet_edges(&brep, &specs, &opts, Tolerance::default(), None);
    /// ```
    pub continuity: ContinuityLevel,
    /// **§13.4 soft-supports mode** for low-precision input
    /// geometry (e.g. lossy STEP imports with vertex drift).
    /// When `true`, the fillet pipeline runs against a relaxed
    /// tolerance derived from the caller's input by scaling the
    /// linear threshold up by `soft_support_factor` (default 10×).
    /// The relaxed tolerance propagates to:
    ///
    /// - The walker's contact solver (loosened Newton residual)
    /// - The §13.1 degeneracy detectors (loosened collapsed-
    ///   edge / near-tangent / zero-area thresholds)
    /// - The §13.1(2) sliver absorber (loosened sliver
    ///   classification thresholds)
    /// - The §13.2 post-construction self-intersection check
    ///   (loosened coincidence threshold)
    ///
    /// **When to use.** Enable for operations on BReps that have
    /// been round-tripped through a lossy format (STEP, IGES,
    /// STL re-meshing) or produced by a boolean op that left
    /// sub-tolerance drift. Leave disabled for authored geometry
    /// that is known to be tight.
    ///
    /// **Always produces a warning** — every relaxed result
    /// carries one
    /// [`crate::error::BlendWarning::SoftSupportTolerancesApplied`]
    /// entry so downstream consumers can distinguish soft-mode
    /// from strict. See §13.4 in `docs/requirements/fillets_blends.md`.
    pub soft_supports: bool,
    /// Linear tolerance multiplication factor applied when
    /// [`Self::soft_supports`] is `true`. Default `10.0` — matches
    /// the common loose-precision convention of "one order
    /// of magnitude looser than the calling tolerance."
    ///
    /// Ignored when `soft_supports = false`.
    pub soft_support_factor: f64,
    /// §8.3 cap type applied to every blend endpoint that
    /// terminates at a free boundary (i.e. doesn't meet another
    /// blend at a shared vertex). Default [`CapKind::Open`] — the
    /// historical behaviour where the blend simply ends at the
    /// face trim with no explicit cap surface.
    ///
    /// The four non-default cap kinds each produce a dedicated
    /// cap surface inserted between the blend surface and the
    /// support faces at the endpoint:
    ///
    /// - [`CapKind::Flat`] — a planar cap perpendicular to the
    ///   spine tangent at the endpoint. G0 with the blend and
    ///   supports. Cheapest and most common for chamfer-style
    ///   closures.
    /// - [`CapKind::Rounded`] — a spherical cap of radius equal
    ///   to the blend's local minor radius, centred at the
    ///   rolling-ball centre-curve endpoint. G1 with the blend
    ///   surface (same radius of curvature at the join).
    /// - [`CapKind::TangentMatch`] — a NURBS cap of revolution
    ///   that continues the strip's ruling with zero curvature
    ///   (Class-A: no curvature jump across the section).
    /// - [`CapKind::UserDefined`] — [BEYOND] the caller provides
    ///   the cap surface. Not yet implemented.
    ///
    /// When a blend endpoint DOES meet another blend (shared
    /// vertex), the `vertex_mode` field drives the join instead
    /// and this field is ignored.
    pub cap_kind: CapKind,
    /// **§5.3 R8 C7 opt-in** — route N=3 trihedral vertex blends
    /// through the rational-Gregory three-quad decomposition
    /// (Chiyokura–Kimura). When `true`, a qualifying corner
    /// (planar supports + manifold hex output) is split into 3
    /// rectangular Gregory sub-patches meeting at a central point;
    /// otherwise the polar-NURBS fallback runs unchanged.
    /// Default `false` — ships wired-but-dormant pending the R8
    /// C8 downstream scorecard. See ADR
    /// `2026-04-17-fillet-rational-gregory-rounds-6-7-8.md`.
    pub use_rational_gregory: bool,
    /// §5.3 item #1 — use quintic (bidegree-5) rectangular Gregory
    /// sub-patches instead of cubic. Default `true` after the
    /// cubic-vs-quintic cross-fixture scorecard showed quintic
    /// reduces §2.4 by 12-17pp across all non-uniform radii
    /// configurations with no regressions. Requires
    /// `use_rational_gregory: true`.
    pub use_deg5_gregory: bool,
    /// §5.3 item #3 extended — Karciauskas-Peters biquintic-full
    /// with spoke-G1 solver. Default `false`; requires
    /// `use_rational_gregory` and `use_deg5_gregory`.
    pub use_biquintic_full_gregory: bool,
    /// §5.3 Phase-9 — biseptic (degree-7) Gregory with spoke-G1
    /// solver; 4 free CPs per spoke (vs 2 at biquintic-full)
    /// gives 4–6× residual reduction. Requires
    /// `use_biquintic_full_gregory`. Default `false`.
    pub use_biseptic_gregory: bool,
    /// §5.3 item #3 Phase-10a — Peters two-stripe construction.
    /// Drops the rational blend at corner-adjacent positions
    /// giving 6 DOF per spoke side (full `{B_1^7..B_6^7}`
    /// coverage) — closes the structural G1 gap. Default
    /// `false`; requires `use_biseptic_gregory`.
    pub use_two_stripe_gregory: bool,
    /// Karciauskas-Peters bidegree-(9, 9) guided
    /// construction replacing two-stripe sub-patches; G1 and
    /// G2 across every shared spoke by construction. Default
    /// `false`; requires `use_two_stripe_gregory`. ADR
    /// `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`.
    pub use_kp_binonic_gregory: bool,
    /// K-P 2015 bi-6 polar construction for N ≥ 5-valent
    /// vertex blends. Default `false`. See the `kp_polar_bi6_promote`
    /// module for usage, scope, diagnostics, and G1 gate plumbing.
    pub use_kp_polar_gregory: bool,
    /// PK-250 serving-boundary self-intersection retry: after the
    /// full pipeline, measure each served direct-apex corner region
    /// and re-run once with offending vertices excluded (vertex-blend
    /// serves them; recorded on the repair report). Default `false`:
    /// default arming is the PK-245 campaign — always-on it exposes
    /// pre-existing hidden serving debt on 9 direct-apex fixtures
    /// (measured 2026-08-07) whose pins accept today's serve.
    pub serving_self_intersection_retry: bool,
    /// **Direct-apex dispatch** — route apex-region
    /// construction (vertices with ≥ 2 fillet edges in Smooth mode)
    /// through the direct-apex planner + builder instead of the vertex-blend ladder
    /// `insert_vertex_blend`. Default `false` as a GLOBAL flag: a
    /// blanket flip mis-routes classes the machinery does not own
    /// (saddles, N ≥ 5 hubs, conic profiles, variable radii —
    /// 2026-08-09, nine fixture failures); serving classes graduate
    /// per-class through `direct_apex_promotion`. The flag remains
    /// the A/B override routing EVERY apex direct (PK-255
    /// measurements: journal 2026-08-09). On plan or build error the
    /// dispatch falls back per-apex to `insert_vertex_blend` and
    /// emits `BlendWarning::DirectApexFellBackToVertexBlend` — never a
    /// hard failure; the worst case is "behaves like vertex-blend."
    pub use_direct_apex_construction: bool,
    /// PK-185 — resolution for mid-face strip interference that is
    /// not seam-completable (the PK-159 analytic completion only
    /// covers parallel-cylinder pairs). Default
    /// [`InterferencePolicy::Decline`] (transactional
    /// `FilletInterference`); see [`InterferencePolicy`] for the
    /// partial-service alternative.
    pub interference_policy: InterferencePolicy,

    /// Cross-section profile (PK-189): circular rolling-ball arc
    /// (default), conic-by-rho, or the curvature-continuous G2
    /// section. See [`BlendProfile`] for the rho semantics and the
    /// circular-subsumption contract.
    pub profile: BlendProfile,
}

impl Default for FilletOptions {
    fn default() -> Self {
        Self {
            propagate_tangent: true,
            vertex_mode: VertexMode::Smooth,
            continuity: ContinuityLevel::G1,
            soft_supports: false,
            soft_support_factor: 10.0,
            cap_kind: CapKind::Open,
            corner_cap: CornerCapPreference::default(),
            corner_setback: CornerSetback::default(),
            // §5.3 R8 C9 — default-ON per the C4b.2c scorecard
            // (Gregory wins 30-70pp vs polar-NURBS on the §2.4
            // metric; uniform-radii falls through to sphere).
            use_rational_gregory: true,
            use_deg5_gregory: true,
            use_biquintic_full_gregory: false,
            use_biseptic_gregory: false,
            use_two_stripe_gregory: false,
            use_kp_binonic_gregory: false,
            use_kp_polar_gregory: false,
            serving_self_intersection_retry: true,
            use_direct_apex_construction: false,
            interference_policy: InterferencePolicy::default(),
            profile: BlendProfile::default(),
        }
    }
}

/// §7 blend-boundary continuity level.
///
/// Each level subsumes the lower ones. See `docs/requirements/fillets_blends.md`
/// §7 for the mathematical conditions and acceptance thresholds.
///
/// - **G0**: positional. `||P_blend - P_face|| < tol_geom` at the boundary.
///   For debugging / special cases.
/// - **G1**: tangent. `|N_blend × N_face| = 0` at the boundary;
///   angular tolerance validates. Default; sufficient for machining
///   and visual smoothness on non-Class-A surfaces.
/// - **G2**: curvature. Second fundamental forms match at the
///   boundary; `kappa_n` continuous in all tangent directions.
///   Required for §15 Class-A surfaces — eliminates highlight-line
///   discontinuities. Acceptance: `|Δκ_n| / max(|κ_n|, 1/r) < 0.01`.
///
/// `G3` (curvature-derivative continuity) is reserved as [BEYOND]
/// and not currently representable here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum ContinuityLevel {
    /// Positional continuity only. Debug / special cases.
    G0,
    /// Tangent continuity — the historical rolling-ball arc sweep.
    G1,
    /// Curvature continuity — degree-5 Bezier sweep with Weingarten-
    /// map normal curvature matching at both contact points.
    G2,
}

impl ContinuityLevel {
    /// commercial-style display string (`"G0"`, `"G1"`, `"G2"`).
    /// Useful for diagnostic messages and UI display.
    pub fn as_str(&self) -> &'static str {
        match self {
            ContinuityLevel::G0 => "G0",
            ContinuityLevel::G1 => "G1",
            ContinuityLevel::G2 => "G2",
        }
    }
}

impl std::fmt::Display for ContinuityLevel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FilletOptions {
    /// §13.4 tolerance relaxation. Returns `(relaxed_tolerance,
    /// optional_warning)` — when soft-supports mode is active and
    /// the factor is > 1, the returned tolerance has `linear` and
    /// `angular` scaled up by the factor, and the warning side
    /// carries a `BlendWarning::SoftSupportTolerancesApplied` that
    /// the caller should attach to the result. Otherwise returns
    /// the caller's tolerance unchanged and `None`.
    pub(crate) fn relax_tolerance(
        &self,
        tolerance: popkern_core::tolerance::Tolerance,
    ) -> (
        popkern_core::tolerance::Tolerance,
        Option<crate::error::BlendWarning>,
    ) {
        if self.soft_supports && self.soft_support_factor > 1.0 {
            let relaxed = tolerance
                .with_linear(tolerance.linear() * self.soft_support_factor)
                .with_angular(tolerance.angular() * self.soft_support_factor);
            let warn = crate::error::BlendWarning::SoftSupportTolerancesApplied {
                factor: self.soft_support_factor,
            };
            (relaxed, Some(warn))
        } else {
            (tolerance, None)
        }
    }
}

/// Vertex blend strategy at multi-fillet junctions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VertexMode {
    /// Smooth G1-continuous vertex blend.
    Smooth,
    /// Mitred (flat) vertex blend.
    Mitred,
}
