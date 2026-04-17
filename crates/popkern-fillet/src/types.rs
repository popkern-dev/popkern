//! API types for fillet and chamfer specifications.

use popkern_core::topology::{EdgeHandle, FaceHandle, VertexHandle};

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
}

/// §5.3 corner-resolution method used to build a single vertex-
/// blend patch. Reported on [`CornerPatchInfo::method`] so a UI
/// can tell at a glance which algorithm fired and colour the
/// corners accordingly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CornerMethod {
    /// Rolling-ball sphere. Used when the incident blend strips
    /// all share one radius and the convexity is uniform (e.g. a
    /// cube corner with equal-radius fillets). G1 by construction.
    RollingBallSphere,
    /// Polar NURBS surface fit through the boundary curves —
    /// fallback when no sphere fits cleanly. Currently G0 with
    /// the incoming blends.
    PolarNurbs,
    /// Flat (mitred) plane fit through the boundary points. Final
    /// legacy fallback; G0, not suitable for Class-A surfaces.
    Mitred,
    /// §5.3 Coons patch for N=4 corners. Bilinear boundary
    /// interpolation with G1 at each boundary. [Pending — C3.]
    Coons,
    /// §5.3 Gregory rational patch for general N. G1 reliable,
    /// G2 implemented for N=3 and N=4 only. [Pending — C3/C4.]
    Gregory,
    /// §5.3 Plate-minimisation (thin-plate energy). [BEYOND, not
    /// currently implemented.]
    PlateMin,
}

impl CornerMethod {
    /// Display string for diagnostics and UI.
    pub fn as_str(&self) -> &'static str {
        match self {
            CornerMethod::RollingBallSphere => "RollingBallSphere",
            CornerMethod::PolarNurbs => "PolarNurbs",
            CornerMethod::Mitred => "Mitred",
            CornerMethod::Coons => "Coons",
            CornerMethod::Gregory => "Gregory",
            CornerMethod::PlateMin => "PlateMin",
        }
    }
}

impl std::fmt::Display for CornerMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// §5.3 per-corner patch report emitted on
/// [`BlendResult::corner_patches`]. One entry per vertex where a
/// blend patch was inserted. A UI uses this to annotate each corner
/// with its method, N, and achieved continuity.
#[derive(Debug, Clone)]
pub struct CornerPatchInfo {
    /// The original BRep vertex the corner was resolved at.
    pub vertex: VertexHandle,
    /// Number of incoming blend boundaries (the N in "N-sided
    /// patch"). 2 for a two-blend corner; 3 for a cube-style
    /// trihedral vertex; 4+ for higher-valence junctions.
    pub n_sides: usize,
    /// Which construction method produced the patch — see
    /// [`CornerMethod`]. A UI should display this alongside the
    /// achieved continuity so callers can distinguish a G1
    /// Rolling-Ball-Sphere (clean analytic output) from a G0
    /// Mitred fallback (flat plane, Class-A failure).
    pub method: CornerMethod,
    /// The [`ContinuityLevel`] actually achieved at the corner-vs-
    /// edge boundaries. Lower than the requested
    /// [`FilletOptions::continuity`] when the method couldn't
    /// deliver — see [`crate::error::BlendWarning::CornerContinuityDowngraded`].
    pub achieved_continuity: ContinuityLevel,
    /// §5.3 C6 — the vertex-blend classification that drove
    /// method dispatch. `SameSideConvex` / `SameSideConcave` are
    /// the comfort-zone cases where a rolling-ball sphere is
    /// tried first; `OppositeSide` / `Mixed` skip the sphere and
    /// fall straight to Gregory (rolling ball can't fit a saddle);
    /// `Tangent` / `FreeBoundary` flag quality caveats a UI
    /// should surface alongside the corner annotation.
    pub classification: crate::vertex_classifier::VertexBlendClassification,
    /// The corner face handle in the output BRep.
    pub face: FaceHandle,
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
    /// **Always produces a warning.** When this flag causes the
    /// tolerance to be relaxed, the result's `edge_warnings`
    /// vector contains exactly one
    /// [`crate::error::BlendWarning::SoftSupportTolerancesApplied`]
    /// entry with the actual factor used, so downstream
    /// consumers can distinguish soft-mode results from strict
    /// ones for logging or auditing.
    ///
    /// See §13.4 "Soft Support Faces" in
    /// `docs/requirements/fillets_blends.md`.
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

/// Specification for a single chamfer edge.
#[derive(Debug, Clone)]
pub struct ChamferSpec {
    /// The edge to chamfer.
    pub edge: EdgeHandle,
    /// The chamfer geometry.
    pub kind: ChamferKind,
}

/// Chamfer geometry definition.
#[derive(Debug, Clone)]
pub enum ChamferKind {
    /// Equal-distance chamfer (45-degree bevel).
    Equal(f64),
    /// Two-distance chamfer with specified face orientation.
    TwoDistance {
        /// Reference face for the first distance.
        face_1: FaceHandle,
        /// Distance on the first face.
        dist_1: f64,
        /// Distance on the second (opposite) face.
        dist_2: f64,
    },
    /// Distance-and-angle chamfer.
    DistanceAngle {
        /// Reference face for the distance measurement.
        face: FaceHandle,
        /// Distance from the edge on the reference face.
        distance: f64,
        /// Angle in radians from the reference face.
        angle: f64,
    },
}

/// Options controlling chamfer behaviour.
#[derive(Debug, Clone)]
pub struct ChamferOptions {
    /// Whether to automatically propagate chamfers along tangent-continuous edge chains.
    pub propagate_tangent: bool,
}

impl Default for ChamferOptions {
    fn default() -> Self {
        Self {
            propagate_tangent: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fillet_options_default() {
        let opts = FilletOptions::default();
        assert!(opts.propagate_tangent);
        assert_eq!(opts.vertex_mode, VertexMode::Smooth);
    }

    #[test]
    fn chamfer_options_default() {
        let opts = ChamferOptions::default();
        assert!(opts.propagate_tangent);
    }

    #[test]
    fn fillet_radius_constant() {
        let r = FilletRadius::Constant(2.5);
        match r {
            FilletRadius::Constant(v) => assert!((v - 2.5).abs() < f64::EPSILON),
            _ => panic!("expected Constant"),
        }
    }
}
