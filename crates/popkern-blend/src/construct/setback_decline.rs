//! Typed outcomes for the setback surgery (PK-574): which stage
//! declined, and which cap construction the surgery commits.
//!
//! The rung needs the stage to choose its next step: a carve decline
//! leaves no region for any cap, while a cap decline (KP20 assembly,
//! acceptance or postcondition) is exactly where the variational fit
//! can serve the same carved region.

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
    /// the certifier and was declined on them.
    pub fit_certified: Option<crate::types::SetbackFitCertified>,
}

impl SetbackDecline {
    /// A cap-stage decline.
    pub(crate) fn cap(error: CapCommitError) -> Self {
        SetbackDecline {
            stage: SetbackStage::Cap,
            reason: error.reason,
            fit_certified: error.fit_certified,
        }
    }
}

/// A cap commit's failure: the reason, and the fit's certified bounds
/// when the failure was a certificate below the fit mode's level.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CapCommitError {
    pub reason: String,
    pub fit_certified: Option<crate::types::SetbackFitCertified>,
}

impl From<String> for CapCommitError {
    fn from(reason: String) -> Self {
        CapCommitError {
            reason,
            fit_certified: None,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SetbackCapStrategy {
    /// The KP20 fan. Its seam-G2 verdict is read by the rung after the
    /// commit (`setback_fan_commit::fan_seam_g2_failure`).
    Kp20,
    /// The variational fit.
    VariationalFit(crate::types::SetbackFitMode),
}
