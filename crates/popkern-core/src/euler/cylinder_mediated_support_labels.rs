//! Diagnostic labels for cylinder-mediated support planning.

use super::cylinder_mediated_support::CylinderMediatedSupportContinuationAggregateDelta;

impl CylinderMediatedSupportContinuationAggregateDelta {
    /// Stable diagnostic label for operator reports and plan logs.
    pub fn diagnostic_label(&self) -> String {
        match self {
            Self::Known {
                vertices,
                edges,
                faces,
                loops,
            } => format!("dV={vertices} dE={edges} dF={faces} dL={loops}"),
            Self::Unknown { reason } => format!("unknown reason={reason}"),
        }
    }
}
