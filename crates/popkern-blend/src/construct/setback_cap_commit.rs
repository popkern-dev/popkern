//! The cap commit on a carved setback region: the KP20 fan or the
//! variational fit (PK-574), by strategy.

use super::*;
use crate::construct::kp20_paper_border_extract::Kp20SetbackArc;
use crate::construct::setback_decline::{CapCommitError, SetbackCapStrategy};

/// Commit the cap `strategy` names over `blend_face`. Returns the cap
/// faces and, for the fit, its certified serve.
#[allow(clippy::too_many_arguments)]
pub(super) fn commit_cap(
    brep: &mut BRep,
    blend_face: FaceHandle,
    side_edges: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    arcs: &[Kp20SetbackArc],
    vertex_position: Point3,
    presolved: Option<&crate::construct::setback_cap_presolve::SetbackCapSolve>,
    strategy: SetbackCapStrategy,
    tolerance: &Tolerance,
) -> Result<(Vec<FaceHandle>, Option<setback_fit_fan::FitServe>), CapCommitError> {
    match strategy {
        SetbackCapStrategy::Kp20 { gate_seam_g2 } => {
            if crate::dispatch_overrides::kp20_cap_decline_forced() {
                return Err("KP20 cap forced to decline (dispatch override)".into());
            }
            super::setback_fan_commit::commit_cap_fan(
                brep,
                blend_face,
                side_edges,
                fillet_surfaces,
                vertex_position,
                presolved,
                gate_seam_g2,
            )
            .map(|faces| (faces, None))
            .map_err(CapCommitError::from)
        }
        SetbackCapStrategy::VariationalFit(mode) => super::setback_fit_fan::commit_fit_fan(
            brep,
            blend_face,
            fillet_surfaces,
            arcs,
            mode,
            tolerance,
        )
        .map(|(faces, serve)| (faces, Some(serve))),
    }
}
