//! Typed outcomes for the setback surgery (PK-574): which stage
//! declined, and which cap construction the surgery commits.
//!
//! The rung needs the stage to choose its next step: a carve decline
//! leaves no region for any cap, while a cap decline (KP20 assembly,
//! acceptance or postcondition) is exactly where the variational fit
//! can serve the same carved region.

/// The crossing carve's decline when a setback arc does not cross its
/// strip's trimmed boundary exactly twice.
pub(crate) const ARC_BOUNDARY_CROSSINGS: &str = "arc/boundary crossings";

/// The marker on [`ARC_BOUNDARY_CROSSINGS`] declines with more than two
/// crossings: the arc runs across its strip's trimmed boundary mid-arc,
/// the corners border overrunning a narrow strip (PK-594), which the
/// bridged border retries (PK-611). Fewer than two crossings is a
/// different class (PK-615: pyramid-3 h0.4 wide at G1 found none) and
/// does not retry.
pub(crate) const ARC_OVERRUN: &str = "a mid-arc crossing overruns the strip";

/// Where the setback surgery declined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SetbackStage {
    /// Border extraction or the topology carve: no carved region.
    Carve,
    /// The cap construction on the carved region.
    Cap,
}

/// A structured setback decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SetbackDecline {
    pub stage: SetbackStage,
    pub reason: String,
    /// The fit's certified bounds when a variational-fit cap reached
    /// the certifier and was declined on them (boxed: the decline is an
    /// error value).
    pub fit_certified: Option<Box<crate::types::SetbackFitCertified>>,
    /// The border the declining attempt was carved on.
    pub border: crate::types::SetbackFitBorder,
    /// The fit's interior reading (in `1 / r_min`) when it was declined
    /// uncertified for being unable to serve (PK-620).
    pub interior: Option<f64>,
}

impl SetbackDecline {
    /// A cap-stage decline.
    pub(crate) fn cap(error: CapCommitError) -> Self {
        SetbackDecline {
            stage: SetbackStage::Cap,
            reason: error.reason,
            fit_certified: error.fit_certified.map(Box::new),
            border: crate::types::SetbackFitBorder::Corners,
            interior: error.interior,
        }
    }
}

/// A cap commit's failure: the reason, and the fit's certified bounds
/// when the failure was a certificate below the fit mode's level.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CapCommitError {
    pub reason: String,
    pub fit_certified: Option<crate::types::SetbackFitCertified>,
    /// The fit's interior reading when it was declined uncertified (see
    /// [`SetbackDecline::interior`]).
    pub interior: Option<f64>,
}

impl From<String> for CapCommitError {
    fn from(reason: String) -> Self {
        CapCommitError {
            reason,
            fit_certified: None,
            interior: None,
        }
    }
}

impl From<&str> for CapCommitError {
    fn from(reason: &str) -> Self {
        CapCommitError::from(reason.to_string())
    }
}

/// Every untyped failure inside the surgery is a carve failure unless
/// the cap commit marks it [`SetbackDecline::cap`].
impl From<String> for SetbackDecline {
    fn from(reason: String) -> Self {
        SetbackDecline {
            stage: SetbackStage::Carve,
            reason,
            fit_certified: None,
            border: crate::types::SetbackFitBorder::Corners,
            interior: None,
        }
    }
}

impl From<&str> for SetbackDecline {
    fn from(reason: &str) -> Self {
        SetbackDecline::from(reason.to_string())
    }
}

impl std::fmt::Display for SetbackDecline {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let stage = match self.stage {
            SetbackStage::Carve => "carve",
            SetbackStage::Cap => "cap",
        };
        write!(f, "{stage}: {}", self.reason)
    }
}

/// Which cap the surgery commits on the carved region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum SetbackCapStrategy {
    /// The KP20 fan. Its seam-G2 verdict is read by the rung after the
    /// commit (`setback_fan_commit::fan_seam_g2_failure`).
    Kp20,
    /// The exact rolling ball (PK-864): the carved region, verified to be
    /// the ball's spherical triangle, takes the sphere directly — no KP20
    /// cap is fitted only to be replaced.
    RollingBall,
    /// The variational fit in `mode`; `r_min` (the corner's smallest
    /// fillet radius) scales the fairness bound its layout refit keys on.
    VariationalFit {
        mode: crate::types::SetbackFitMode,
        r_min: f64,
        /// The border to carve first; a corners border that overruns a
        /// narrow strip retries on the bridged border (PK-611).
        border: crate::types::SetbackFitBorder,
        /// The interior reading (in `1 / r_min`) of an unfair fan the rung
        /// holds: a thin fitted sector serves while it stays under the
        /// larger of this and the fairness bound (PK-621).
        rival: Option<f64>,
        /// The reading a fit must beat to serve (the fit the ladder
        /// holds, or an unfair fan's): a fit no fairer is declined
        /// before certification (PK-620).
        held: Option<f64>,
        /// Whether a fit past the fairness bound defers its proof to
        /// the ladder's end (PK-620).
        defer: bool,
    },
}
