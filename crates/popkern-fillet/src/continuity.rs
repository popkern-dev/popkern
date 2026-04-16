//! §7 continuity-downgrade plumbing.
//!
//! Collected in one module so the warning-emission and fallback-
//! tagging logic don't bloat `lib.rs` and `geometry_compute.rs`
//! past the 500-line structural cap. See
//! [`crate::types::ContinuityLevel`] and
//! [`crate::error::BlendWarning::ContinuityDowngraded`] for the
//! public surface.

use popkern_core::topology::EdgeHandle;

use crate::construct::EdgeGeometry;
use crate::error::BlendWarning;
use crate::types::{ContinuityLevel, FilletOptions};

/// Append one `ContinuityDowngraded` warning per edge whose
/// `EdgeGeometry::continuity_downgrade` is populated.
pub(crate) fn emit_downgrade_warnings(
    edge_geometries: &[EdgeGeometry],
    options: &FilletOptions,
    warnings: &mut Vec<(EdgeHandle, BlendWarning)>,
) {
    for geom in edge_geometries {
        if let Some((achieved, reason)) = &geom.continuity_downgrade {
            let edge = geom.analysis.edge;
            warnings.push((
                edge,
                BlendWarning::ContinuityDowngraded {
                    edge,
                    requested: options.continuity,
                    achieved: *achieved,
                    reason: reason.clone(),
                },
            ));
        }
    }
}

/// Promote a walker/sweep `downgrade_reason` string into the typed
/// `(achieved, reason)` pair stored on `EdgeGeometry`. Returns
/// `None` when no downgrade occurred. The achieved level is always
/// `G1` in the current pipeline — the fallback paths land on arc
/// sweeps / analytic torus output.
pub(crate) fn to_field(downgrade_reason: Option<String>) -> Option<(ContinuityLevel, String)> {
    downgrade_reason.map(|r| (ContinuityLevel::G1, r))
}
