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
        SetbackCapStrategy::Kp20 => {
            if crate::dispatch_overrides::kp20_cap_decline_forced() {
                return Err("KP20 cap forced to decline (dispatch override)".into());
            }
            // The ranked assemblies (PK-613): the fairest first, then each
            // alternate on a fan postcondition failure, from the same carve.
            let alternates = presolved.map_or(&[][..], |solve| &solve.alternates[..]);
            let snapshot = (!alternates.is_empty()).then(|| brep.clone());
            let commit = |brep: &mut BRep,
                          solve: Option<
                &crate::construct::setback_cap_presolve::SetbackCapSolve,
            >| {
                super::setback_fan_commit::commit_cap_fan(
                    brep,
                    blend_face,
                    side_edges,
                    fillet_surfaces,
                    vertex_position,
                    solve,
                )
            };
            let first = match commit(brep, presolved) {
                Ok(faces) => return Ok((faces, None)),
                Err(e) => e,
            };
            if let (Some(solve), Some(snapshot)) = (presolved, snapshot) {
                for alternate in alternates {
                    *brep = snapshot.clone();
                    let retry = crate::construct::setback_cap_presolve::SetbackCapSolve {
                        border: solve.border.clone(),
                        extract: solve.extract.clone(),
                        assembly: Ok(alternate.clone()),
                        alternates: Vec::new(),
                    };
                    if let Ok(faces) = commit(brep, Some(&retry)) {
                        return Ok((faces, None));
                    }
                }
                *brep = snapshot;
            }
            Err(CapCommitError::from(first))
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
