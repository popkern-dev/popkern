//! Non-fatal per-edge warnings for fillet/chamfer operations.
//!
//! Extracted from `error.rs` to keep each module under the 500-LoC
//! structural limit and to give Phase I (polar K-P) headroom for
//! its fallback-ladder warning variants.
//!
//! Re-exported as `crate::error::BlendWarning` for API stability;
//! callers never see the module boundary.

use popkern_core::topology::{EdgeHandle, VertexHandle};

/// Non-fatal per-edge warning.
#[derive(Debug, Clone)]
pub enum BlendWarning {
    /// Edge was auto-propagated from a chain.
    Propagated,
    /// Radius was clamped to fit available space.
    RadiusClamped {
        /// Originally requested radius.
        requested: f64,
        /// Actual radius used after clamping.
        actual: f64,
    },
    /// Vertex mode was downgraded because the requested mode is not yet
    /// supported for this configuration.
    VertexModeFallback {
        /// The mode that was requested.
        requested: crate::types::VertexMode,
        /// The mode that was actually used.
        actual: crate::types::VertexMode,
        /// The vertex where the fallback occurred, if identified.
        vertex: Option<VertexHandle>,
    },
    /// Two fillet trim curves overlap on a shared face.
    FilletOverlap {
        /// Index of the first fillet edge in the specs list.
        fillet_a: usize,
        /// Index of the second fillet edge in the specs list.
        fillet_b: usize,
    },
    /// The blend output surface may have self-intersection risk at
    /// the requested radius — the caller should consider reducing
    /// the radius or running a post-operation validation.
    SelfIntersectionRisk {
        /// The radius that triggered the warning.
        radius: f64,
        /// The estimated safe maximum radius.
        safe_radius: f64,
    },
    /// A fillet spec was soft-skipped from the batch by the
    /// §13.1(5) collapsed-edge pre-pass because the target edge's
    /// two endpoint vertices were separated by less than the
    /// geometric tolerance. The rest of the batch proceeds
    /// normally; the caller should check `warnings` in addition
    /// to `failures` to understand the full set of spec outcomes.
    CollapsedEdgeSkipped {
        /// The edge that was dropped.
        edge: EdgeHandle,
        /// Measured 3D length between the edge's endpoint vertices.
        length: f64,
    },
    /// A sliver face was detected in the fillet output topology but
    /// could not be absorbed into a neighbour (no suitable target
    /// or the topology surgery would be invalid). The output BRep
    /// contains the sliver unchanged; downstream consumers may
    /// choose to ignore it, run a separate healing pass, or fail.
    ///
    /// §13.1(2) "remove, stitch neighbors" — commercial kernels
    /// ship this warning for the residual cases where absorption
    /// is impossible. See `construct::sliver_absorb` for the
    /// classifier details.
    SliverRetained {
        /// The residual sliver face.
        face: popkern_core::topology::FaceHandle,
        /// Short description of the sliver category that was
        /// classified but not absorbed.
        kind: &'static str,
    },
    /// §13.4 soft-supports mode was active for this fillet
    /// operation: the pipeline's internal tolerance was scaled up
    /// from the caller's value to accommodate low-precision input
    /// geometry. Downstream consumers can distinguish soft-mode
    /// results from strict ones by checking for this warning.
    ///
    /// Emitted exactly once per fillet call when
    /// [`crate::types::FilletOptions::soft_supports`] is `true`.
    SoftSupportTolerancesApplied {
        /// The multiplicative factor applied to `tolerance.linear()`
        /// when producing the relaxed tolerance. Defaults to 10,
        /// set via
        /// [`crate::types::FilletOptions::soft_support_factor`].
        factor: f64,
    },
    /// §8.3 cap was produced at a blend endpoint. Carries the
    /// requested kind and a flag describing whether topology was
    /// modified to insert the cap face (on an open-shell or
    /// free-boundary endpoint) or whether the cap is represented
    /// implicitly by an existing adjacent solid face (the common
    /// case on closed bodies, where the adjacent third face at a
    /// corner already serves as the natural "flat" cap).
    ///
    /// Emitted once per `(edge, endpoint)` pair whenever
    /// [`crate::types::FilletOptions::cap_kind`] is non-`Open`.
    /// Callers can inspect the per-edge `EdgeGeometry::start_cap` /
    /// `end_cap` fields for the raw cap geometry (surface +
    /// boundary polyline) regardless of the warning.
    CapApplied {
        /// The filleted edge the cap terminates.
        edge: EdgeHandle,
        /// `true` if the cap is at the spine's start endpoint
        /// (`t = 0`); `false` if at the end (`t = 1`).
        at_start: bool,
        /// The cap kind the caller requested.
        kind: crate::types::CapKind,
        /// `true` if cap topology was inserted into the output
        /// BRep; `false` if the endpoint was already closed by an
        /// adjacent face and no insertion was needed.
        topology_inserted: bool,
    },
    /// §7 continuity fallback: the caller requested a higher
    /// continuity level than the engine was able to deliver for
    /// this edge. The output is still a valid fillet — it just
    /// sits at a lower level than requested.
    ///
    /// **When this fires.** The primary trigger is a G2 request
    /// where either (a) the spine walker couldn't produce contact
    /// solutions for the surface pair (e.g. degenerate normals,
    /// analytic plane-plane path that bypasses the walker), or
    /// (b) `g2_sweep_from_contact_solutions` errored on at least
    /// one spine sample (`G2Error::ZeroChord`,
    /// `DegenerateNormal`, or `TangentDegenerate`).
    ///
    /// **UI guidance.** Display this alongside
    /// [`crate::quality::BlendQuality::achieved_continuity`] so the
    /// caller sees both "you asked for G2" and "you got G1 because
    /// the cross-section at sample 3 had a zero chord." Consumer
    /// code that requires strict G2 should treat this as a failure
    /// and reject the output.
    ContinuityDowngraded {
        /// The filleted edge whose blend surface was downgraded.
        edge: EdgeHandle,
        /// What the caller asked for (from
        /// [`crate::types::FilletOptions::continuity`]).
        requested: crate::types::ContinuityLevel,
        /// What the engine actually delivered on this edge's blend
        /// face. Always lower than `requested`.
        achieved: crate::types::ContinuityLevel,
        /// Human-readable reason for the downgrade. Comes from the
        /// underlying walker / G2-sweep error so UI consumers can
        /// surface the specific geometric condition that forced
        /// the fallback.
        reason: String,
    },
    /// §5.3 corner-patch continuity fallback at a high-valence
    /// vertex. Fires when either (a) the construction reported an
    /// achieved continuity lower than `requested`, or (b) the
    /// corner has five or more incident fillet strips
    /// (`n_sides >= 5`), where even the curvature-matched polar
    /// NURBS patch degrades in interior quality because a single
    /// degree-2 radial cross-section can't capture the full
    /// curvature variation around N boundaries.
    ///
    /// **UI guidance.** Display alongside the
    /// [`crate::CornerPatchInfo`] for the same vertex so the caller
    /// sees `(method = Gregory, n_sides = 6, achieved = G1)`
    /// paired with "corner continuity downgraded because N ≥ 5."
    /// Consumer code that requires strict Class-A at N≥5 should
    /// treat this warning as a failure and reject the corner
    /// (Plate-minimisation / true rational Gregory for N≥5 is
    /// planned for §5.3 C6+).
    ///
    /// Keyed on a representative incident edge on the outer
    /// `Vec<(EdgeHandle, BlendWarning)>` to keep the warning
    /// plumbing uniform with other vertex-scoped warnings
    /// (`VertexModeFallback`); `vertex` inside the body carries
    /// the authoritative handle.
    CornerContinuityDowngraded {
        /// The vertex whose corner patch was downgraded.
        vertex: VertexHandle,
        /// Valence of the corner — number of incident filleted
        /// edges meeting at this vertex.
        n_sides: usize,
        /// The continuity level the caller asked for (from
        /// [`crate::types::FilletOptions::continuity`]).
        requested: crate::types::ContinuityLevel,
        /// The continuity level the corner construction actually
        /// delivered.
        achieved: crate::types::ContinuityLevel,
        /// Human-readable reason (e.g. "N >= 5, polar NURBS
        /// interior may degrade").
        reason: String,
    },
    /// §5.3 item #6 — a face-edge bridge in the vertex-blend
    /// boundary construction fell back to a straight-line cubic
    /// Bezier because the generalised surface-hugging bridge
    /// couldn't be built on the face's supporting surface (e.g.,
    /// degenerate UV projection at a singular NURBS point).
    ///
    /// The straight-line bridge still interpolates the arc
    /// endpoints correctly and the downstream Gregory fit
    /// re-samples the bridge — so the blend surface is valid
    /// but may sit slightly *off* the supporting face near the
    /// bridge, rather than hugging it. Surfaces downstream of
    /// the blend (mesh tessellation, boolean consumers) may
    /// show a visible seam where the bridge leaves the face.
    ///
    /// Callers can surface this warning in a UI as "curved-
    /// support blend fell back to planar approximation" — in
    /// the common case (planar supports), this variant is never
    /// emitted.
    SurfaceBridgeFallback {
        /// The vertex whose blend boundary triggered the fallback.
        vertex: VertexHandle,
        /// Short description of which surface type the fallback
        /// occurred on (e.g., `"Cylinder"`, `"NURBS"`) — for
        /// logging / UI display.
        surface_kind: &'static str,
    },
    /// §5.3 item #7 — a tangent corner (dihedral ≈ π at one or
    /// more incident edges) was resolved by a planar cap with
    /// G1-tangent-match to the fillet arcs, rather than a curved
    /// blend. Informational, not a failure: at a tangent
    /// configuration the rolling-ball contact degenerates so the
    /// planar cap is the geometrically-correct choice. A UI can
    /// surface this as a neutral "tangent join" badge on the
    /// affected corner.
    TangentCornerCap {
        /// The tangent-classified corner vertex.
        vertex: VertexHandle,
        /// Valence of the corner — number of incident filleted edges.
        n_sides: usize,
    },
    /// §5.3 item #10 — the rational-Gregory three-quad dispatch
    /// produced a sub-patch that self-intersects (fold-over or
    /// curl-around detected by
    /// [`GregoryPatch::self_intersects_with_tolerance`](popkern_core::geometry::gregory::GregoryPatch::self_intersects_with_tolerance)).
    /// The corner was remediated by falling back to the single-
    /// patch polar-NURBS construction, which always produces a
    /// manifold (but lower-quality) blend. UIs can surface this as
    /// an "overlap remediation" badge so the user knows the corner
    /// used the safe fallback instead of the preferred Gregory
    /// geometry.
    CornerRemediationOverlap {
        /// The vertex where the Gregory construction self-
        /// intersected and the remediation was applied.
        vertex: VertexHandle,
        /// Valence of the corner.
        n_sides: usize,
    },
    /// Phase H.8 step-D — the K-P bidegree-9 construction shipped
    /// the corner but could not reach the commercial-grade G1
    /// tolerance. commercial-style graceful degradation: ship the
    /// best-effort patch rather than failing. Common under extreme
    /// radius asymmetry. UIs surface as "reduced quality" badge.
    CornerContinuityRelaxed {
        /// The polar vertex.
        vertex: VertexHandle,
        /// Valence.
        n_sides: usize,
        /// Worst-case G1 coplanarity on the interior
        /// `t ∈ [0.05, 0.95]` dense sampling.
        worst_coplanarity: f64,
    },
    /// Phase I — polar K-P bi-nonic unavailable; bicubic polar rung
    /// took over. Ships G1 but not the intended bi-nonic quality.
    PolarKpBinonicToBicubic {
        /// The polar vertex (valence ≥ 5).
        vertex: VertexHandle,
        /// Valence.
        n_sides: usize,
        /// Worst-case G1 coplanarity achieved.
        worst_coplanarity: f64,
    },
    /// Phase I — polar K-P failed; sphere-fit shortcut took over.
    /// Common on uniform-radii convex N-valence vertices.
    PolarKpToSphereFit {
        /// The polar vertex.
        vertex: VertexHandle,
        /// Valence.
        n_sides: usize,
    },
    /// Phase I — polar K-P failed; single-patch polar-NURBS took
    /// over. Significant quality drop vs K-P sector construction.
    PolarKpToPolarNurbs {
        /// The polar vertex.
        vertex: VertexHandle,
        /// Valence.
        n_sides: usize,
        /// Worst-case G1 coplanarity on the polar-NURBS surface.
        worst_coplanarity: f64,
    },
    /// Phase I — all upstream options failed; mitred plane G0-only
    /// last resort. Indicates severe geometry / numerical issues.
    PolarKpToMitredPlane {
        /// The polar vertex that fell all the way to mitre.
        vertex: VertexHandle,
        /// Valence.
        n_sides: usize,
    },
    /// **Phase 6c-direct-apex dispatch** — apex-region construction was requested
    /// to use the direct-apex builder
    /// ([`crate::types::FilletOptions::use_direct_apex_construction`]
    /// = `true`) but the direct-apex path could not finish for this apex,
    /// so the dispatch fell back to the legacy
    /// `insert_vertex_blend` per-vertex.  The output BRep is the
    /// same as the legacy path would produce; this warning is the
    /// only signal that direct-apex was tried first.
    ///
    /// Common reasons during the direct-apex build.b/c → direct-apex build.d/e/f transition
    /// (when the direct-apex builder errors on Sphere/Gregory paths that
    /// haven't been implemented yet) — the warning's `reason`
    /// field carries the direct-apex error message verbatim.
    DirectApexFellBackToLegacy {
        /// The apex vertex where direct-apex was attempted.
        vertex: VertexHandle,
        /// The direct-apex plan or build error that triggered fallback,
        /// rendered through its `Display` impl for end-user
        /// consumption.  Non-empty.
        reason: String,
    },
}

/// Trace a smooth-vertex-blend fallback to Mitred under the direct-apex
/// trace lever — the insertion error is otherwise swallowed by the
/// warning-only fallback arm.
pub(crate) fn trace_smooth_blend_mitred_fallback(
    vertex: popkern_core::topology::VertexHandle,
    err: &str,
) {
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] smooth vertex blend at {vertex:?} fell back to Mitred: {err}"
        );
    }
}
