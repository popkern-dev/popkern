//! API types for fillet and chamfer specifications.

#[path = "types_chamfer.rs"]
mod types_chamfer;
#[path = "types_corner.rs"]
mod types_corner;
#[path = "types_repair.rs"]
mod types_repair;

pub use types_chamfer::{ChamferKind, ChamferOptions, ChamferSpec};
pub use types_corner::{
    CornerMethod, CornerPatchInfo, CornerSubPatchFamily, CornerSubPatchProvenance,
    CornerSubPatchRole,
};
pub use types_repair::{
    ApexConstructionDiagnostic, ApexConstructionRoute, ApexConstructionRouteDecision,
    BlendRepairReport, CapSupportContinuityDiagnostic, CylinderTrimFarSideStrategy,
    CylinderTrimMaterializationDeclineReason, CylinderTrimMaterializationDiagnostic,
    ExactBoundaryCapDiagnostic, SupportSeamDecision, SupportSeamDeclineReason,
    SupportSeamRelocationReport, SupportSeamRollbackReason, build_blend_repair_action_artifact,
    support_seam_decline_reason_label, support_seam_rollback_reason_label,
};

use popkern_core::topology::EdgeHandle;

use crate::error::BlendWarning;

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
    ///   legacy circular-arc rolling-ball sweep; sufficient for
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
    /// use popkern_fillet::{ContinuityLevel, FilletOptions, FilletSpec, FilletRadius, fillet_edges};
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
    /// - [`CapKind::TangentMatch`] — a NURBS cap whose boundary
    ///   tangents match both support-face tangents at the blend
    ///   endpoint contact points. G1 with the supports and the
    ///   blend. Most expensive, highest quality.
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
    /// otherwise the legacy polar-NURBS fallback runs unchanged.
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
    /// Phase H — Karciauskas-Peters bidegree-(9, 9) guided
    /// construction replacing two-stripe sub-patches; G1 and
    /// G2 across every shared spoke by construction. Default
    /// `false`; requires `use_two_stripe_gregory`. ADR
    /// `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`.
    pub use_kp_binonic_gregory: bool,
    /// Phase I — K-P 2015 bi-6 polar construction for N ≥ 5-valent
    /// vertex blends. Default `false`. See the `kp_polar_bi6_promote`
    /// module for usage, scope, diagnostics, and G1 gate plumbing.
    pub use_kp_polar_gregory: bool,
    /// **Phase 6c-direct-apex dispatch** — route apex-region construction (vertices
    /// with ≥ 2 incident fillet edges in `VertexMode::Smooth`)
    /// through the direct-apex planner + builder (`plan_apex_resolution`
    /// then `construct_apex` in `crate::apex_construction`) instead
    /// of the legacy `insert_vertex_blend` path.  Default `false`
    /// — the legacy path runs unchanged until direct-apex default flips the
    /// default after the validation matrix is green.
    ///
    /// On direct-apex plan or build error, the dispatch falls back per-apex
    /// to `insert_vertex_blend` and emits
    /// [`BlendWarning::DirectApexFellBackToLegacy`](crate::error::BlendWarning::DirectApexFellBackToLegacy)
    /// (added in direct-apex dispatch.c).  Never causes a hard failure; the worst
    /// case is "behaves like legacy."
    pub use_direct_apex_construction: bool,
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
            use_direct_apex_construction: false,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
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

/// §8.3 blend cap type for free-boundary endpoints.
///
/// When a filleted edge terminates without meeting another blend
/// (i.e. at a face boundary, a cut, or an open sheet edge), the
/// engine applies one of five strategies to close the blend end.
///
/// The default is [`CapKind::Open`] — no explicit cap surface —
/// which preserves the historical behaviour for callers that do
/// not set [`FilletOptions::cap_kind`]. The other variants each
/// produce a dedicated cap face matching a specific continuity
/// requirement (G0 for flat, G1 for rounded and tangent-match).
///
/// See [`FilletOptions::cap_kind`] for the motivation and
/// continuity guarantees of each kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapKind {
    /// No cap surface — the blend ends at its natural trim boundary
    /// and the free edge is reported to the caller as-is. This is
    /// the default and matches the pre-§8.3 behaviour.
    Open,
    /// Planar cap perpendicular to the spine tangent at the blend
    /// endpoint. G0 with the blend and supports. The cheapest cap
    /// kind and the most common for chamfer-style end closures.
    Flat,
    /// Hemispherical cap of radius equal to the blend's local
    /// minor radius, centred at the rolling-ball centre-curve
    /// endpoint. G1 with the blend surface (matching radii of
    /// curvature at the join).
    Rounded,
    /// NURBS cap surface whose boundary tangents match both
    /// support-face tangents at the blend endpoint contact points.
    /// G1 with the supports and the blend. The most expensive
    /// cap kind and the highest quality for Class-A surfacing.
    TangentMatch,
    /// [BEYOND] The caller supplies the cap surface via a callback
    /// or resource handle. Placeholder — not yet implemented;
    /// requesting this kind currently falls back to [`Self::Open`]
    /// with a warning.
    UserDefined,
}
