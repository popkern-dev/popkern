//! Blend-option enums split out of `types.rs` (500-line cap):
//! the cross-section profile (PK-189), the strip-interference
//! resolution policy (PK-185), the end-cap kind (§8.3), the corner
//! construction preference (PK-603) and the user corner setback
//! (PK-623).

/// Cross-section profile of a blend (requirements 2.2, PK-189).
///
/// The rolling-ball contact curves and end tangents are identical
/// for every profile — only the section INTERIOR differs — so
/// profiles compose with variable radius and with the corner-cap
/// dispatch exactly like the circular sweep does.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum BlendProfile {
    /// Circular rolling-ball section (default; the historical
    /// behaviour, analytic torus/cylinder fast paths eligible).
    #[default]
    Circular,
    /// Conic (rational quadratic) section with an explicit shoulder
    /// factor. `rho` is the classic conic fullness parameter in
    /// (0, 1): the section's shoulder point sits at
    /// `M + rho·(A − M)` where `M` is the contact-chord midpoint
    /// and `A` the tangent-intersection apex. `rho = 0.5` is
    /// parabolic; smaller is flatter, larger is fuller. At each
    /// station the CIRCULAR section corresponds to
    /// `rho = cos(θ/2) / (1 + cos(θ/2))` for section opening angle
    /// `θ`, and supplying that value reproduces the circular blend
    /// to roundoff — the conic path subsumes the arc, it does not
    /// fork it.
    ///
    /// Conic sections are exact rational-quadratic rows — there is
    /// no approximation or fit. Stations whose opening angle
    /// exceeds the single-segment rational range decline with a
    /// structured error (no silent downgrade to circular).
    Conic {
        /// Shoulder fullness in (0, 1); see the variant docs.
        rho: f64,
    },
    /// Curvature-continuous degree-5 section — the existing G2
    /// sweep. Selecting this is equivalent to requesting
    /// [`crate::types::ContinuityLevel::G2`].
    CurvatureContinuous,
}

impl BlendProfile {
    /// Rational-quadratic middle weight for a conic profile:
    /// `w = rho / (1 − rho)`.
    pub fn conic_weight(rho: f64) -> f64 {
        rho / (1.0 - rho)
    }

    /// The rho value at which a conic section reproduces the
    /// circular arc for a station with opening angle
    /// `theta` (radians): `cos(θ/2) / (1 + cos(θ/2))`.
    pub fn circular_rho(theta: f64) -> f64 {
        let w = (0.5 * theta).cos();
        w / (1.0 + w)
    }
}

/// PK-185 — resolution policy for non-completable mid-face strip
/// interference. See [`crate::types::FilletOptions::interference_policy`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum InterferencePolicy {
    /// Decline the whole operation transactionally (default — the
    /// historical behaviour).
    #[default]
    Decline,
    /// Drop both members of each interfering pair and serve the
    /// remaining fillets; every dropped pair is reported with a
    /// [`BlendWarning::FilletPairDropped`](crate::error::BlendWarning::FilletPairDropped)
    /// carrying the reason. When nothing remains after dropping,
    /// the operation still declines (there is no partial result to
    /// serve).
    ServeNonInterfering,
}

/// §8.3 blend cap type for free-boundary endpoints.
///
/// When a filleted edge terminates without meeting another blend
/// (i.e. at a face boundary, a cut, or an open sheet edge), the
/// engine applies one of five strategies to close the blend end.
///
/// The default is [`CapKind::Open`] — no explicit cap surface —
/// which preserves the historical behaviour for callers that do
/// not set [`crate::FilletOptions::cap_kind`]. The other variants each
/// produce a dedicated cap face matching a specific continuity
/// requirement (G0 for flat, G1 for rounded and tangent-match).
///
/// See [`crate::FilletOptions::cap_kind`] for the motivation and
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
    /// NURBS cap of revolution about the end tangent, its meridian
    /// continuing the strip's ruling with zero curvature to a flat
    /// pole: G2 with straight constant-radius strips, G1 with the
    /// supports. See [`crate::cap::build_tangent_match_cap`].
    TangentMatch,
    /// [BEYOND] The caller supplies the cap surface via a callback
    /// or resource handle. Placeholder — not yet implemented;
    /// requesting this kind currently falls back to [`Self::Open`]
    /// with a warning.
    UserDefined,
}

/// Which construction closes a setback-eligible fillet corner (three or
/// more rounds meeting, smooth corner mode) — PK-603. Corners the
/// setback rung never sees are unaffected.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CornerConstruction {
    /// The kernel's ladder: the setback G2 cap, its measured G1
    /// downgrade, the fitted cap, then the classic vertex-blend paths.
    #[default]
    Auto,
    /// The setback G2 cap first (as `Auto`, stated): direct-apex
    /// promotion does not pre-empt it.
    Setback,
    /// The fitted cap (the variational corner fit, PK-574) first, then
    /// the setback cap, then the classic paths.
    Fitted,
    /// Skip the setback rung: the classic vertex-blend paths (direct
    /// apex, boundary ring, exact-boundary ring, Gregory, mitre).
    Classic,
}

/// The corner construction preference of a fillet (PK-603).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CornerCapPreference {
    /// The preferred construction.
    pub construction: CornerConstruction,
    /// Serve the preferred construction or fail the fillet with
    /// [`crate::BlendFailure::CornerCapUnavailable`], never falling
    /// back. Only `Setback` and `Fitted` can decline; `Auto` and
    /// `Classic` have nothing to refuse.
    pub strict: bool,
}

impl CornerConstruction {
    /// The scripting spelling: `auto`, `setback`, `fitted`, `classic`.
    pub fn spelling(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Setback => "setback",
            Self::Fitted => "fitted",
            Self::Classic => "classic",
        }
    }
}

impl CornerCapPreference {
    /// Whether the preference names a setback-rung construction, so
    /// direct-apex promotion must not claim a rung-eligible corner
    /// before the rung runs.
    pub fn names_setback_rung(&self) -> bool {
        matches!(
            self.construction,
            CornerConstruction::Setback | CornerConstruction::Fitted
        )
    }

    /// Whether a decline of the preferred construction fails the
    /// fillet.
    pub fn declines_when_unavailable(&self) -> bool {
        self.strict && self.names_setback_rung()
    }
}

/// User-set setback distances at the corners the setback construction
/// serves (PK-623, `docs/architecture/fillet-corner-setback.md`).
///
/// A distance is the arclength along a round's two rails from where the
/// sharp edge's vertex projects onto them to the round's end section;
/// on a straight edge, the distance along the edge from the vertex. A
/// corner with any user setback is served at exactly those stations or
/// the fillet fails ([`crate::BlendFailure::SetbackOutOfRange`],
/// [`crate::BlendFailure::CornerCapUnavailable`]); it never falls back
/// to an automatic setback.
///
/// ```
/// use popkern_blend::CornerSetback;
/// // Every round set back 0.15 at every setback corner.
/// let uniform = CornerSetback::uniform(0.15);
/// assert_eq!(uniform.distance, Some(0.15));
/// assert!(CornerSetback::default().is_automatic());
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CornerSetback {
    /// Every round's setback at every corner the setback construction
    /// serves; `None` keeps the automatic setbacks.
    pub distance: Option<f64>,
    /// Per round end, overriding `distance` there.
    pub ends: Vec<EdgeEndSetback>,
}

/// One round end's user setback (PK-623).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeEndSetback {
    /// The corner vertex.
    pub vertex: popkern_core::topology::VertexHandle,
    /// The filleted edge ending at `vertex`.
    pub edge: popkern_core::topology::EdgeHandle,
    /// The setback distance (model units).
    pub distance: f64,
}

impl CornerSetback {
    /// Every round set back `distance` at every setback corner.
    pub fn uniform(distance: f64) -> Self {
        CornerSetback {
            distance: Some(distance),
            ends: Vec::new(),
        }
    }

    /// This setback with `edge`'s end at `vertex` set back `distance`.
    pub fn with_end(
        mut self,
        vertex: popkern_core::topology::VertexHandle,
        edge: popkern_core::topology::EdgeHandle,
        distance: f64,
    ) -> Self {
        self.ends.push(EdgeEndSetback {
            vertex,
            edge,
            distance,
        });
        self
    }

    /// Whether every setback is the kernel's own.
    pub fn is_automatic(&self) -> bool {
        self.distance.is_none() && self.ends.is_empty()
    }

    /// The user setback of `edge`'s end at `vertex`, if any: its own
    /// override, else the uniform distance.
    pub fn at(
        &self,
        vertex: popkern_core::topology::VertexHandle,
        edge: popkern_core::topology::EdgeHandle,
    ) -> Option<f64> {
        self.ends
            .iter()
            .rev()
            .find(|e| e.vertex == vertex && e.edge == edge)
            .map(|e| e.distance)
            .or(self.distance)
    }

    /// Whether an override names `vertex`.
    pub fn names_vertex(&self, vertex: popkern_core::topology::VertexHandle) -> bool {
        self.ends.iter().any(|e| e.vertex == vertex)
    }
}

impl CornerSetback {
    /// Check the setbacks before any geometry work: finite positive
    /// distances; each override on a filleted edge (`filleted`) at one of
    /// its own vertices; a smooth, setback-capable corner mode.
    ///
    /// # Errors
    ///
    /// [`crate::BlendError::InvalidParameter`] naming the offending value.
    pub fn validate(
        &self,
        brep: &popkern_core::topology::BRep,
        filleted: &[popkern_core::topology::EdgeHandle],
        vertex_mode: &crate::types::VertexMode,
        corner_cap: CornerCapPreference,
    ) -> Result<(), crate::BlendError> {
        if self.is_automatic() {
            return Ok(());
        }
        let invalid = |value: String, reason: &str| crate::BlendError::InvalidParameter {
            parameter: "options.corner_setback".to_string(),
            value,
            reason: reason.to_string(),
        };
        if *vertex_mode != crate::types::VertexMode::Smooth {
            return Err(invalid(
                format!("{vertex_mode:?}"),
                "a mitred corner has no setback; use the smooth corner mode",
            ));
        }
        if corner_cap.construction == CornerConstruction::Classic {
            return Err(invalid(
                "classic".to_string(),
                "the classic corner paths do not set rounds back; drop the classic corner cap",
            ));
        }
        let distances = self
            .distance
            .iter()
            .chain(self.ends.iter().map(|e| &e.distance));
        if let Some(d) = distances.copied().find(|d| !(d.is_finite() && *d > 0.0)) {
            return Err(invalid(
                format!("{d}"),
                "a setback distance must be finite and positive",
            ));
        }
        for end in &self.ends {
            if !filleted.contains(&end.edge) {
                return Err(invalid(
                    format!("{:?}", end.edge),
                    "the setback's edge is not one of the filleted edges",
                ));
            }
            let at_vertex = brep
                .half_edges
                .values()
                .any(|h| h.edge == end.edge && h.origin == end.vertex);
            if !at_vertex {
                return Err(invalid(
                    format!("{:?} at {:?}", end.edge, end.vertex),
                    "the setback's vertex is not an end of its edge",
                ));
            }
        }
        Ok(())
    }
}
